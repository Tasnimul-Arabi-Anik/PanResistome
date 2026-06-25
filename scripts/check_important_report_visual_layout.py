#!/usr/bin/env python3
"""Visual-layout QA for PanResistome important/results.html.

This checker is intentionally report-focused rather than biology-focused. It
verifies that the generated HTML has the expected dashboard structure, figure
cards, captions, download links, accessibility basics, responsive CSS, and
optional browser-rendered layout checks.
"""

from __future__ import annotations

import argparse
import html.parser
import re
import sys
from pathlib import Path
from urllib.parse import urlparse

from check_important_report_outputs import REQUIRED_SECTION_IDS, _linked_files, _resolve_paths


REQUIRED_UI_CLASSES = [
    "report-header",
    "sidebar",
    "sidebar-links",
    "section",
    "summary-card",
    "figure-card",
    "figure-guidance",
    "section-focus",
    "table-card",
    "analysis-card",
    "warning-box",
    "download-card",
    "back-to-top",
]

REQUIRED_CSS_TOKENS = [
    ":root",
    "--primary: #0f766e",
    "--warning: #f97316",
    "overflow-x: hidden",
    "max-width: 100%",
    "@media (max-width: 920px)",
    "@media (max-width: 520px)",
    ".table-scroll",
    ".figure-card",
    ".figure-guidance",
    ".section-focus",
    ".download-bar",
    ".download-card-grid",
    ".details-block",
]

GENERIC_CAPTION_PHRASES = [
    "Report-facing figure with PNG",
    "Report-facing figure with companion downloads",
    "Figure preview with companion PNG",
    "Supporting report figure",
]

VIEWPORTS = {
    "desktop": (1440, 1000),
    "laptop": (1180, 820),
    "mobile": (390, 844),
}


def _attrs_to_dict(attrs: list[tuple[str, str | None]]) -> dict[str, str]:
    return {key: value or "" for key, value in attrs}


def _classes(attrs: dict[str, str]) -> set[str]:
    return {token.strip() for token in attrs.get("class", "").split() if token.strip()}


class ReportHTMLScanner(html.parser.HTMLParser):
    """Collect report UI structure while keeping figure-card context."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.ids: set[str] = set()
        self.class_counts: dict[str, int] = {}
        self.tags: list[tuple[str, dict[str, str]]] = []
        self.figures: list[dict[str, object]] = []
        self._active_figures: list[dict[str, object]] = []
        self._stack: list[tuple[str, dict[str, object] | None]] = []
        self.details_count = 0
        self.summary_count = 0
        self.section_count = 0
        self.h1_count = 0
        self.h2_count = 0
        self.img_without_alt: list[str] = []
        self.iframe_without_title: list[str] = []

    def handle_starttag(self, tag: str, attrs_raw: list[tuple[str, str | None]]) -> None:
        attrs = _attrs_to_dict(attrs_raw)
        classes = _classes(attrs)
        self.tags.append((tag, attrs))
        if "id" in attrs:
            self.ids.add(attrs["id"])
        for class_name in classes:
            self.class_counts[class_name] = self.class_counts.get(class_name, 0) + 1
        if tag == "section":
            self.section_count += 1
        if tag == "h1":
            self.h1_count += 1
        if tag == "h2":
            self.h2_count += 1
        if tag == "details":
            self.details_count += 1
        if tag == "summary":
            self.summary_count += 1
        if tag == "img" and not attrs.get("alt", "").strip():
            self.img_without_alt.append(attrs.get("src", ""))
        if tag == "iframe" and not attrs.get("title", "").strip():
            self.iframe_without_title.append(attrs.get("src", ""))

        figure_obj: dict[str, object] | None = None
        if tag == "div" and "figure-card" in classes:
            figure_obj = {
                "h3": 0,
                "captions": 0,
                "images": 0,
                "downloads": 0,
                "png": 0,
                "svg": 0,
                "pdf": 0,
                "data": 0,
                "missing_alt": [],
            }
            self.figures.append(figure_obj)
            self._active_figures.append(figure_obj)
        if self._active_figures:
            figure = self._active_figures[-1]
            if tag == "h3":
                figure["h3"] = int(figure["h3"]) + 1
            if tag == "p" and "figure-caption" in classes:
                figure["captions"] = int(figure["captions"]) + 1
            if tag == "img":
                figure["images"] = int(figure["images"]) + 1
                if not attrs.get("alt", "").strip():
                    missing_alt = figure["missing_alt"]
                    assert isinstance(missing_alt, list)
                    missing_alt.append(attrs.get("src", ""))
            if tag == "div" and (("figure-downloads" in classes) or ("download-bar" in classes)):
                figure["downloads"] = int(figure["downloads"]) + 1
            if tag == "a":
                href = attrs.get("href", "")
                suffix = Path(urlparse(href).path).suffix.lower()
                if suffix == ".png":
                    figure["png"] = int(figure["png"]) + 1
                elif suffix == ".svg":
                    figure["svg"] = int(figure["svg"]) + 1
                elif suffix == ".pdf":
                    figure["pdf"] = int(figure["pdf"]) + 1
                elif suffix in {".tsv", ".csv"}:
                    figure["data"] = int(figure["data"]) + 1
        self._stack.append((tag, figure_obj))

    def handle_endtag(self, tag: str) -> None:
        while self._stack:
            open_tag, figure_obj = self._stack.pop()
            if figure_obj is not None and figure_obj in self._active_figures:
                self._active_figures.remove(figure_obj)
            if open_tag == tag:
                break


def _scan_html(html_text: str) -> ReportHTMLScanner:
    scanner = ReportHTMLScanner()
    scanner.feed(html_text)
    return scanner


def _check_static_layout(html_text: str, scanner: ReportHTMLScanner, errors: list[str]) -> None:
    if "file://" in html_text:
        errors.append("results.html contains file:// links")
    for phrase in GENERIC_CAPTION_PHRASES:
        if phrase in html_text:
            errors.append(f"results.html contains generic caption phrase: {phrase}")

    missing_sections = [section_id for section_id in REQUIRED_SECTION_IDS if section_id not in scanner.ids]
    if missing_sections:
        errors.append("Missing section anchors: " + ", ".join(missing_sections))
    if scanner.h1_count < 1:
        errors.append("results.html is missing a semantic h1")
    if scanner.h2_count < len(REQUIRED_SECTION_IDS):
        errors.append(
            f"results.html has {scanner.h2_count} h2 headings; expected at least {len(REQUIRED_SECTION_IDS)}"
        )

    missing_classes = [class_name for class_name in REQUIRED_UI_CLASSES if scanner.class_counts.get(class_name, 0) == 0]
    if missing_classes:
        errors.append("Missing UI classes: " + ", ".join(missing_classes))
    missing_css = [token for token in REQUIRED_CSS_TOKENS if token not in html_text]
    if missing_css:
        errors.append("Missing responsive/style tokens: " + ", ".join(missing_css))

    if scanner.img_without_alt:
        errors.append(f"{len(scanner.img_without_alt)} image(s) are missing alt text")
    if scanner.iframe_without_title:
        errors.append(f"{len(scanner.iframe_without_title)} iframe(s) are missing title attributes")
    if scanner.details_count == 0 or scanner.summary_count == 0:
        errors.append("results.html has no collapsible details/summary blocks for dense sections")

    if not scanner.figures:
        errors.append("results.html has no figure cards")
    for index, figure in enumerate(scanner.figures, 1):
        missing_parts: list[str] = []
        for field, label in [
            ("h3", "title"),
            ("captions", "caption"),
            ("images", "image"),
            ("downloads", "download block"),
            ("png", "PNG link"),
            ("svg", "SVG link"),
            ("data", "data link"),
        ]:
            if int(figure[field]) < 1:
                missing_parts.append(label)
        if missing_parts:
            errors.append(f"Figure card {index} is missing: {', '.join(missing_parts)}")
        missing_alt = figure["missing_alt"]
        if isinstance(missing_alt, list) and missing_alt:
            errors.append(f"Figure card {index} has image(s) without alt text")

    table_cards = scanner.class_counts.get("table-card", 0)
    table_scrolls = scanner.class_counts.get("table-scroll", 0)
    if table_cards and table_scrolls < max(1, table_cards // 2):
        errors.append(
            f"Only {table_scrolls} table-scroll containers for {table_cards} table cards; dense tables may spill into the page"
        )


def _run_browser_checks(
    important_dir: Path,
    screenshots_dir: Path | None,
    section_screenshots: bool,
    max_overflow_px: int,
    require_browser: bool,
    errors: list[str],
    notes: list[str],
) -> None:
    try:
        from playwright.sync_api import Error as PlaywrightError
        from playwright.sync_api import sync_playwright
    except Exception as exc:  # pragma: no cover - depends on optional local tooling
        message = f"Browser layout QA skipped because Playwright is unavailable: {exc}"
        if require_browser:
            errors.append(message)
        else:
            notes.append(message)
        return

    screenshots_dir = screenshots_dir.resolve() if screenshots_dir else None
    if screenshots_dir:
        screenshots_dir.mkdir(parents=True, exist_ok=True)

    js_metrics = """
    () => {
      const allowedOverflowSelectors = [
        '.table-scroll',
        '.figure-box',
        '.heatmap-scroll',
        'pre',
        'code',
        'iframe'
      ];
      const allowed = (el) => allowedOverflowSelectors.some(sel => el.closest(sel));
      const doc = document.documentElement;
      const body = document.body;
      const overflow = Math.max(doc.scrollWidth, body.scrollWidth) - window.innerWidth;
      const offenders = [];
      for (const el of Array.from(document.body.querySelectorAll('*'))) {
        const rect = el.getBoundingClientRect();
        if (rect.width <= 0 || rect.height <= 0 || allowed(el)) continue;
        if (rect.right > window.innerWidth + 4 || rect.left < -4) {
          offenders.push({
            tag: el.tagName.toLowerCase(),
            id: el.id || '',
            cls: el.className ? String(el.className).slice(0, 120) : '',
            left: Math.round(rect.left),
            right: Math.round(rect.right),
            width: Math.round(rect.width)
          });
        }
        if (offenders.length >= 8) break;
      }
      const blankFigures = Array.from(document.querySelectorAll('.figure-card img')).filter(img => {
        const rect = img.getBoundingClientRect();
        return rect.width < 40 || rect.height < 30;
      }).map(img => img.getAttribute('src')).slice(0, 8);
      const header = document.querySelector('.report-header');
      const sidebar = document.querySelector('.sidebar');
      return {
        overflow,
        offenders,
        blankFigures,
        figureCards: document.querySelectorAll('.figure-card').length,
        tableCards: document.querySelectorAll('.table-card').length,
        downloadCards: document.querySelectorAll('.download-card').length,
        headerVisible: !!header && header.getBoundingClientRect().height > 40,
        sidebarVisible: !!sidebar && sidebar.getBoundingClientRect().height > 40
      };
    }
    """

    try:  # pragma: no cover - exercised only when local browser tooling exists
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch()
            for name, (width, height) in VIEWPORTS.items():
                page = browser.new_page(viewport={"width": width, "height": height}, device_scale_factor=1)
                page.goto((important_dir / "results.html").as_uri(), wait_until="load")
                page.wait_for_timeout(500)
                metrics = page.evaluate(js_metrics)
                if metrics["overflow"] > max_overflow_px:
                    errors.append(
                        f"{name} viewport has {metrics['overflow']:.0f}px whole-page horizontal overflow"
                    )
                if metrics["offenders"]:
                    errors.append(f"{name} viewport has overflowing elements: {metrics['offenders']}")
                if metrics["blankFigures"]:
                    errors.append(f"{name} viewport has tiny/blank figure previews: {metrics['blankFigures']}")
                if not metrics["headerVisible"]:
                    errors.append(f"{name} viewport does not render a visible report header")
                if not metrics["sidebarVisible"]:
                    errors.append(f"{name} viewport does not render a visible sidebar/navigation")
                if screenshots_dir:
                    page.screenshot(path=str(screenshots_dir / f"results_{name}.png"), full_page=False)
                    if section_screenshots:
                        for section_id in ["featured", "geography", "metadata-associations", "lineage", "downloads"]:
                            locator = page.locator(f"section#{section_id}")
                            if locator.count():
                                locator.screenshot(path=str(screenshots_dir / f"{name}_{section_id}.png"))
                page.close()
            browser.close()
    except PlaywrightError as exc:  # pragma: no cover - depends on optional local tooling
        message = f"Browser layout QA could not run: {exc}"
        if require_browser:
            errors.append(message)
        else:
            notes.append(message)


def validate(
    sample_dir: Path,
    browser: str = "auto",
    screenshots_dir: Path | None = None,
    section_screenshots: bool = False,
    max_overflow_px: int = 8,
) -> tuple[list[str], list[str]]:
    _root_dir, important_dir = _resolve_paths(sample_dir)
    html_text = (important_dir / "results.html").read_text(encoding="utf-8", errors="ignore")
    scanner = _scan_html(html_text)
    errors: list[str] = []
    notes: list[str] = []

    _check_static_layout(html_text, scanner, errors)

    for link in sorted(_linked_files(html_text)):
        target = (important_dir / link).resolve()
        if not target.exists():
            errors.append(f"Broken relative link: {link}")

    if browser != "skip":
        _run_browser_checks(
            important_dir,
            screenshots_dir,
            section_screenshots,
            max_overflow_px,
            browser == "require",
            errors,
            notes,
        )
    return errors, notes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample_dir", type=Path, help="Sample output directory, important directory, or results.html path")
    parser.add_argument(
        "--browser",
        choices=["auto", "skip", "require"],
        default="auto",
        help="Run optional Playwright browser checks when available. Use 'require' to fail if browser QA cannot run.",
    )
    parser.add_argument("--screenshots-dir", type=Path, help="Optional directory for viewport screenshots")
    parser.add_argument("--section-screenshots", action="store_true", help="Also capture key section screenshots")
    parser.add_argument("--max-overflow-px", type=int, default=8, help="Allowed whole-page horizontal overflow in browser checks")
    args = parser.parse_args(argv)

    errors, notes = validate(
        args.sample_dir,
        browser=args.browser,
        screenshots_dir=args.screenshots_dir,
        section_screenshots=args.section_screenshots,
        max_overflow_px=args.max_overflow_px,
    )
    for note in notes:
        print(f"NOTE: {note}", file=sys.stderr)
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print("important report visual-layout QA passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
