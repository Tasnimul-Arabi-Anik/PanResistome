import sys
import unittest
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import check_important_report_visual_layout as checker


class ImportantReportVisualLayoutCheckerTests(unittest.TestCase):
    def test_find_chrome_executable_uses_candidate_order(self):
        def fake_which(name):
            return f"/usr/bin/{name}" if name == "chromium" else None

        with mock.patch("shutil.which", side_effect=fake_which):
            self.assertEqual(
                checker._find_chrome_executable(("missing-browser", "chromium", "google-chrome")),
                "/usr/bin/chromium",
            )

    def test_chrome_screenshot_command_uses_safe_headless_flags(self):
        command = checker._chrome_screenshot_command(
            "/usr/bin/google-chrome",
            "file:///tmp/report.html",
            Path("/tmp/report.png"),
            390,
            844,
        )

        self.assertEqual(command[0], "/usr/bin/google-chrome")
        for flag in ["--headless", "--disable-gpu", "--disable-dev-shm-usage", "--no-sandbox"]:
            self.assertIn(flag, command)
        self.assertIn("--window-size=390,844", command)
        self.assertIn("--screenshot=/tmp/report.png", command)
        self.assertEqual(command[-1], "file:///tmp/report.html")

    def test_extract_chrome_metrics_unescapes_dumped_dom_json(self):
        metrics = checker._extract_chrome_metrics(
            '<html><body><pre id="qa-metrics">{&quot;overflow&quot;:0,&quot;headerVisible&quot;:true}</pre></body></html>'
        )

        self.assertEqual(metrics["overflow"], 0)
        self.assertTrue(metrics["headerVisible"])

    def test_section_scroll_script_targets_requested_anchor(self):
        script = checker._scroll_to_section_script("geography")

        self.assertIn('"geography"', script)
        self.assertIn("scrollIntoView", script)
        self.assertIn("window.scrollBy", script)

    def test_browser_require_fails_when_no_engine_is_available(self):
        errors = []
        notes = []
        with (
            mock.patch("check_important_report_visual_layout._run_playwright_browser_checks", return_value=False),
            mock.patch("check_important_report_visual_layout._find_chrome_executable", return_value=None),
        ):
            checker._run_browser_checks(
                Path("/tmp"),
                "<html><body></body></html>",
                None,
                False,
                8,
                True,
                errors,
                notes,
            )

        self.assertTrue(errors)
        self.assertIn("no Playwright or Chrome/Chromium", errors[0])

    def test_browser_auto_records_note_when_no_engine_is_available(self):
        errors = []
        notes = []
        with (
            mock.patch("check_important_report_visual_layout._run_playwright_browser_checks", return_value=False),
            mock.patch("check_important_report_visual_layout._find_chrome_executable", return_value=None),
        ):
            checker._run_browser_checks(
                Path("/tmp"),
                "<html><body></body></html>",
                None,
                False,
                8,
                False,
                errors,
                notes,
            )

        self.assertFalse(errors)
        self.assertTrue(any("no Playwright or Chrome/Chromium" in note for note in notes))


if __name__ == "__main__":
    unittest.main()
