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
import json
import re
import shutil
import subprocess
import sys
import tempfile
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
    "best-figure-note",
    "interactive-explorer-card",
    "explorer-frame",
    "report-storyline",
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
    ".best-figure-note",
    ".interactive-explorer-card",
    ".explorer-frame",
    ".report-storyline",
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

SECTION_SCREENSHOT_IDS = ["featured", "geography", "metadata-associations", "lineage", "downloads"]
CHROME_CANDIDATES = ("google-chrome-stable", "google-chrome", "chromium", "chromium-browser")
CHROME_BASE_FLAGS = [
    "--headless",
    "--disable-gpu",
    "--disable-dev-shm-usage",
    "--no-sandbox",
    "--run-all-compositor-stages-before-draw",
    "--virtual-time-budget=3000",
]

QA_METRICS_JS = """
(() => {
  const sectionIds = __SECTION_IDS__;
  const allowedOverflowSelectors = [
    '.table-scroll',
    '.figure-box',
    '.heatmap-scroll',
    'pre',
    'code',
    'iframe',
    '.sidebar-links'
  ];
  const allowed = (el) => allowedOverflowSelectors.some(sel => el.closest(sel));
  const textOf = (el) => (el.textContent || '').trim().replace(/\\s+/g, ' ').slice(0, 120);
  const measure = () => {
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
          text: textOf(el),
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
    const clippedHeaderDownloads = Array.from(document.querySelectorAll('.report-header .download-button, .report-header .downloads a')).filter(el => {
      const rect = el.getBoundingClientRect();
      return rect.width <= 0 || rect.right > window.innerWidth + 4 || rect.left < -4;
    }).map(textOf).slice(0, 8);
    const sections = {};
    for (const id of sectionIds) {
      const section = document.getElementById(id);
      const rect = section ? section.getBoundingClientRect() : null;
      sections[id] = {
        present: !!section,
        height: rect ? Math.round(rect.height) : 0,
        textLength: section ? textOf(section).length : 0
      };
    }
    return {
      overflow,
      offenders,
      blankFigures,
      clippedHeaderDownloads,
      sections,
      figureCards: document.querySelectorAll('.figure-card').length,
      tableCards: document.querySelectorAll('.table-card').length,
      downloadCards: document.querySelectorAll('.download-card').length,
      headerVisible: !!header && header.getBoundingClientRect().height > 40,
      sidebarVisible: !!sidebar && sidebar.getBoundingClientRect().height > 40
    };
  };
  const finish = () => {
    const metrics = measure();
    document.body.innerHTML = '<pre id="qa-metrics"></pre>';
    document.getElementById('qa-metrics').textContent = JSON.stringify(metrics);
  };
  if (document.readyState === 'complete') {
    setTimeout(finish, 500);
  } else {
    window.addEventListener('load', () => setTimeout(finish, 500));
  }
})();
"""

SCROLL_TO_ANCHOR_JS = """
(() => {
  const scrollToHash = () => {
    const id = decodeURIComponent((window.location.hash || '').replace(/^#/, ''));
    if (!id) return;
    const target = document.getElementById(id);
    if (target) target.scrollIntoView({block: 'start'});
  };
  if (document.readyState === 'complete') {
    setTimeout(scrollToHash, 500);
  } else {
    window.addEventListener('load', () => setTimeout(scrollToHash, 500));
  }
})();
"""


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
        self.iframes: list[dict[str, str]] = []

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
        if tag == "iframe":
            self.iframes.append(attrs)
            if not attrs.get("title", "").strip():
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
    iframe_without_lazy = [
        attrs.get("src", "")
        for attrs in scanner.iframes
        if attrs.get("loading", "").strip().lower() != "lazy"
    ]
    if iframe_without_lazy:
        errors.append(f"{len(iframe_without_lazy)} iframe(s) are not lazy loaded")
    if "Load embedded explorer" not in html_text:
        errors.append("results.html does not expose collapsed embedded explorer controls")
    if "Best figure to start with" not in html_text:
        errors.append("results.html does not include best-figure guidance notes")
    if "Report storyline" not in html_text:
        errors.append("results.html does not include a report storyline card")
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


def _check_browser_metrics(metrics: dict[str, object], name: str, max_overflow_px: int, errors: list[str]) -> None:
    overflow = float(metrics.get("overflow", 0) or 0)
    if overflow > max_overflow_px:
        errors.append(f"{name} viewport has {overflow:.0f}px whole-page horizontal overflow")
    offenders = metrics.get("offenders") or []
    if offenders:
        errors.append(f"{name} viewport has overflowing elements: {offenders}")
    blank_figures = metrics.get("blankFigures") or []
    if blank_figures:
        errors.append(f"{name} viewport has tiny/blank figure previews: {blank_figures}")
    clipped_downloads = metrics.get("clippedHeaderDownloads") or []
    if clipped_downloads:
        errors.append(f"{name} viewport has clipped header download buttons: {clipped_downloads}")
    if not metrics.get("headerVisible"):
        errors.append(f"{name} viewport does not render a visible report header")
    if not metrics.get("sidebarVisible"):
        errors.append(f"{name} viewport does not render a visible sidebar/navigation")
    sections = metrics.get("sections") or {}
    if isinstance(sections, dict):
        for section_id in SECTION_SCREENSHOT_IDS:
            section = sections.get(section_id, {})
            if not isinstance(section, dict) or not section.get("present"):
                errors.append(f"{name} viewport cannot find section #{section_id}")
            elif int(section.get("height", 0) or 0) < 20 or int(section.get("textLength", 0) or 0) < 10:
                errors.append(f"{name} viewport has blank/tiny section #{section_id}")


def _playwright_metrics_js() -> str:
    return """
    () => {
      const sectionIds = __SECTION_IDS__;
      const allowedOverflowSelectors = [
        '.table-scroll',
        '.figure-box',
        '.heatmap-scroll',
        'pre',
        'code',
        'iframe',
        '.sidebar-links'
      ];
      const allowed = (el) => allowedOverflowSelectors.some(sel => el.closest(sel));
      const textOf = (el) => (el.textContent || '').trim().replace(/\\s+/g, ' ').slice(0, 120);
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
            text: textOf(el),
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
      const clippedHeaderDownloads = Array.from(document.querySelectorAll('.report-header .download-button, .report-header .downloads a')).filter(el => {
        const rect = el.getBoundingClientRect();
        return rect.width <= 0 || rect.right > window.innerWidth + 4 || rect.left < -4;
      }).map(textOf).slice(0, 8);
      const sections = {};
      for (const id of sectionIds) {
        const section = document.getElementById(id);
        const rect = section ? section.getBoundingClientRect() : null;
        sections[id] = {
          present: !!section,
          height: rect ? Math.round(rect.height) : 0,
          textLength: section ? textOf(section).length : 0
        };
      }
      return {
        overflow,
        offenders,
        blankFigures,
        clippedHeaderDownloads,
        sections,
        figureCards: document.querySelectorAll('.figure-card').length,
        tableCards: document.querySelectorAll('.table-card').length,
        downloadCards: document.querySelectorAll('.download-card').length,
        headerVisible: !!header && header.getBoundingClientRect().height > 40,
        sidebarVisible: !!sidebar && sidebar.getBoundingClientRect().height > 40
      };
    }
    """.replace("__SECTION_IDS__", json.dumps(SECTION_SCREENSHOT_IDS))


def _run_playwright_browser_checks(
    important_dir: Path,
    screenshots_dir: Path | None,
    section_screenshots: bool,
    max_overflow_px: int,
    errors: list[str],
    notes: list[str],
) -> bool:
    try:
        from playwright.sync_api import Error as PlaywrightError
        from playwright.sync_api import sync_playwright
    except Exception as exc:  # pragma: no cover - depends on optional local tooling
        notes.append(f"Playwright unavailable: {exc}")
        return False

    try:  # pragma: no cover - exercised only when local browser tooling exists
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch()
            for name, (width, height) in VIEWPORTS.items():
                page = browser.new_page(viewport={"width": width, "height": height}, device_scale_factor=1)
                page.goto((important_dir / "results.html").as_uri(), wait_until="load")
                page.wait_for_timeout(500)
                metrics = page.evaluate(_playwright_metrics_js())
                _check_browser_metrics(metrics, name, max_overflow_px, errors)
                if screenshots_dir:
                    page.screenshot(path=str(screenshots_dir / f"results_{name}.png"), full_page=False)
                    if section_screenshots:
                        for section_id in SECTION_SCREENSHOT_IDS:
                            locator = page.locator(f"section#{section_id}")
                            if locator.count():
                                locator.screenshot(path=str(screenshots_dir / f"{name}_{section_id}.png"))
                page.close()
            browser.close()
        notes.append("Browser layout QA used playwright")
        return True
    except PlaywrightError as exc:  # pragma: no cover - depends on optional local tooling
        notes.append(f"Playwright browser layout QA could not run: {exc}")
        return False


def _find_chrome_executable(candidates: tuple[str, ...] = CHROME_CANDIDATES) -> str | None:
    for candidate in candidates:
        executable = shutil.which(candidate)
        if executable:
            return executable
    return None


def _chrome_screenshot_command(executable: str, url: str, output_path: Path, width: int, height: int) -> list[str]:
    return [
        executable,
        *CHROME_BASE_FLAGS,
        f"--window-size={width},{height}",
        f"--screenshot={output_path}",
        url,
    ]


def _chrome_dump_dom_command(executable: str, url: str, width: int, height: int) -> list[str]:
    return [
        executable,
        *CHROME_BASE_FLAGS,
        f"--window-size={width},{height}",
        "--dump-dom",
        url,
    ]


def _html_with_base_and_script(html_text: str, base_uri: str, script: str) -> str:
    base_tag = f'<base href="{base_uri}">\n'
    if re.search(r"<head[^>]*>", html_text, flags=re.IGNORECASE):
        html_text = re.sub(r"(<head[^>]*>)", r"\1\n" + base_tag, html_text, count=1, flags=re.IGNORECASE)
    else:
        html_text = base_tag + html_text
    script_tag = f"\n<script>{script}</script>\n"
    if re.search(r"</body>", html_text, flags=re.IGNORECASE):
        return re.sub(r"</body>", lambda _match: script_tag + "</body>", html_text, count=1, flags=re.IGNORECASE)
    return html_text + script_tag


def _extract_chrome_metrics(dom_text: str) -> dict[str, object] | None:
    match = re.search(r"<pre[^>]*id=[\"']qa-metrics[\"'][^>]*>(.*?)</pre>", dom_text, flags=re.DOTALL)
    if not match:
        return None
    text = html.unescape(match.group(1))
    try:
        value = json.loads(text)
    except json.JSONDecodeError:
        return None
    return value if isinstance(value, dict) else None


def _run_chrome_browser_checks(
    important_dir: Path,
    html_text: str,
    screenshots_dir: Path | None,
    section_screenshots: bool,
    max_overflow_px: int,
    errors: list[str],
    notes: list[str],
) -> bool:
    executable = _find_chrome_executable()
    if not executable:
        notes.append("Chrome/Chromium unavailable")
        return False

    base_uri = important_dir.resolve().as_uri().rstrip("/") + "/"
    metrics_js = QA_METRICS_JS.replace("__SECTION_IDS__", json.dumps(SECTION_SCREENSHOT_IDS))
    screenshot_js = SCROLL_TO_ANCHOR_JS
    with tempfile.TemporaryDirectory(prefix="panr2_report_browser_qa_") as tmpdir:
        tmp = Path(tmpdir)
        metrics_html = tmp / "metrics.html"
        metrics_html.write_text(_html_with_base_and_script(html_text, base_uri, metrics_js), encoding="utf-8")
        screenshot_html = tmp / "screenshot.html"
        screenshot_html.write_text(_html_with_base_and_script(html_text, base_uri, screenshot_js), encoding="utf-8")

        for name, (width, height) in VIEWPORTS.items():
            command = _chrome_dump_dom_command(executable, metrics_html.resolve().as_uri(), width, height)
            result = subprocess.run(command, check=False, capture_output=True, text=True, timeout=60)
            if result.returncode != 0:
                errors.append(f"Chrome browser QA failed for {name}: {result.stderr.strip()[:500]}")
                continue
            metrics = _extract_chrome_metrics(result.stdout)
            if metrics is None:
                errors.append(f"Chrome browser QA could not parse metrics for {name}")
                continue
            _check_browser_metrics(metrics, name, max_overflow_px, errors)
            if screenshots_dir:
                screenshot_command = _chrome_screenshot_command(
                    executable,
                    screenshot_html.resolve().as_uri(),
                    screenshots_dir / f"results_{name}.png",
                    width,
                    height,
                )
                screenshot_result = subprocess.run(screenshot_command, check=False, capture_output=True, text=True, timeout=60)
                if screenshot_result.returncode != 0:
                    errors.append(f"Chrome screenshot failed for {name}: {screenshot_result.stderr.strip()[:500]}")

        if screenshots_dir and section_screenshots:
            for section_id in SECTION_SCREENSHOT_IDS:
                command = _chrome_screenshot_command(
                    executable,
                    screenshot_html.resolve().as_uri() + f"#{section_id}",
                    screenshots_dir / f"section_{section_id}.png",
                    1440,
                    1000,
                )
                result = subprocess.run(command, check=False, capture_output=True, text=True, timeout=60)
                if result.returncode != 0:
                    errors.append(f"Chrome section screenshot failed for #{section_id}: {result.stderr.strip()[:500]}")
    notes.append(f"Browser layout QA used chrome: {executable}")
    return True


def _run_browser_checks(
    important_dir: Path,
    html_text: str,
    screenshots_dir: Path | None,
    section_screenshots: bool,
    max_overflow_px: int,
    require_browser: bool,
    errors: list[str],
    notes: list[str],
) -> None:
    screenshots_dir = screenshots_dir.resolve() if screenshots_dir else None
    if screenshots_dir:
        screenshots_dir.mkdir(parents=True, exist_ok=True)

    if _run_playwright_browser_checks(
        important_dir,
        screenshots_dir,
        section_screenshots,
        max_overflow_px,
        errors,
        notes,
    ):
        return
    if _run_chrome_browser_checks(
        important_dir,
        html_text,
        screenshots_dir,
        section_screenshots,
        max_overflow_px,
        errors,
        notes,
    ):
        return

    message = "Browser layout QA skipped because no Playwright or Chrome/Chromium engine was available"
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
            html_text,
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
