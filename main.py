"""
Protein Hunter v5
Main entry point.
"""

from __future__ import annotations

import time

from core.startup import StartupChecker


def main() -> None:
    """Run Protein Hunter."""

    checker = StartupChecker()

    if not checker.run():
        raise SystemExit(1)

    from core.logger import logger

    logger.info("Protein Hunter started")
    logger.success("Startup check passed")

    with logger.section("Logger Test"):
        with logger.timer("Example step"):
            time.sleep(0.5)

    logger.summary()

    logger.success("Program finished")


if __name__ == "__main__":
    main()