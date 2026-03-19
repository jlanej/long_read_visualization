#!/usr/bin/env python3
"""Capture screenshots of the visualization server for documentation.

This script uses Selenium with a headless Chromium browser to capture
screenshots of the IGV.js multi-panel visualization server.  It is
designed to be run from the ``generate-docs`` GitHub Actions workflow
but can also be invoked locally.

Usage
-----
    python3 scripts/capture_screenshots.py \\
        --url http://localhost:8080 \\
        --output-dir docs/screenshots
"""

import argparse
import os
import sys
import time


def main():
    parser = argparse.ArgumentParser(
        description="Capture screenshots of the visualization server")
    parser.add_argument("--url", default="http://localhost:8080",
                        help="Server URL (default: http://localhost:8080)")
    parser.add_argument("--output-dir", default="docs/screenshots",
                        help="Output directory for screenshots")
    parser.add_argument("--wait", type=int, default=8,
                        help="Seconds to wait for page load (default: 8)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    try:
        from selenium import webdriver
        from selenium.webdriver.chrome.options import Options
        from selenium.webdriver.chrome.service import Service
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
    except ImportError:
        print("selenium is required: pip install selenium", file=sys.stderr)
        # Generate placeholder documentation even without screenshots
        _generate_placeholder(args.output_dir)
        return

    chrome_options = Options()
    chrome_options.add_argument("--headless")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--window-size=1920,1080")
    chrome_options.add_argument("--disable-gpu")

    driver = None
    try:
        driver = webdriver.Chrome(options=chrome_options)
        driver.set_window_size(1920, 1080)

        # 1. Capture main view (initial load)
        print("Capturing main view...")
        driver.get(args.url)
        time.sleep(args.wait)
        driver.save_screenshot(
            os.path.join(args.output_dir, "01_main_view.png"))

        # 2. Capture with first region loaded
        print("Capturing first region view...")
        time.sleep(2)
        driver.save_screenshot(
            os.path.join(args.output_dir, "02_region_view.png"))

        # 3. Navigate to next region and capture
        print("Capturing second region...")
        try:
            next_btn = driver.find_element(By.ID, "nextRegion")
            next_btn.click()
            time.sleep(3)
            driver.save_screenshot(
                os.path.join(args.output_dir, "03_next_region.png"))
        except Exception as e:
            print(f"  Warning: Could not navigate to next region: {e}")

        print(f"Screenshots saved to {args.output_dir}/")

    except Exception as e:
        print(f"Error capturing screenshots: {e}", file=sys.stderr)
        _generate_placeholder(args.output_dir)
    finally:
        if driver:
            driver.quit()


def _generate_placeholder(output_dir):
    """Create a placeholder text file when screenshots cannot be captured."""
    path = os.path.join(output_dir, "README.md")
    with open(path, "w") as f:
        f.write("# Screenshots\n\n")
        f.write("Screenshots will be generated automatically when the\n")
        f.write("`generate-docs` workflow runs on the `main` branch.\n")
    print(f"Wrote placeholder to {path}")


if __name__ == "__main__":
    main()
