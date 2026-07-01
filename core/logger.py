"""
Protein Hunter v5
Production logger module.

This module provides:
- Colored console logging
- Timestamped file logging
- logs/latest.log
- SUCCESS log level
- Section logging
- Runtime timers
- Exception logging with traceback
"""

from __future__ import annotations

import logging
import sys
import time
import traceback
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import Iterator

from colorama import Fore, Style, init

from core.constants import LOG_DIR

init(autoreset=True)

SUCCESS_LEVEL = 25
logging.addLevelName(SUCCESS_LEVEL, "SUCCESS")


class ConsoleFormatter(logging.Formatter):
    """Colored console formatter."""

    COLORS = {
        "DEBUG": Fore.CYAN,
        "INFO": Fore.BLUE,
        "SUCCESS": Fore.GREEN,
        "WARNING": Fore.YELLOW,
        "ERROR": Fore.RED,
        "CRITICAL": Fore.RED + Style.BRIGHT,
    }

    ICONS = {
        "DEBUG": "●",
        "INFO": "ℹ",
        "SUCCESS": "✓",
        "WARNING": "⚠",
        "ERROR": "✗",
        "CRITICAL": "✗",
    }

    def format(self, record: logging.LogRecord) -> str:
        color = self.COLORS.get(record.levelname, "")
        icon = self.ICONS.get(record.levelname, "•")
        time_text = datetime.fromtimestamp(record.created).strftime("%H:%M:%S")
        message = record.getMessage()

        return f"{color}[{time_text}] {icon} {message}{Style.RESET_ALL}"


class FileFormatter(logging.Formatter):
    """Plain text file formatter."""

    def __init__(self) -> None:
        super().__init__(
            fmt="%(asctime)s | %(levelname)-8s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )


class LoggerManager:
    """Main logger manager for Protein Hunter."""

    def __init__(
        self,
        name: str = "ProteinHunter",
        log_dir: Path = LOG_DIR,
        level: int = logging.INFO,
    ) -> None:
        self.name = name
        self.log_dir = Path(log_dir)
        self.level = level
        self.timers: dict[str, float] = {}

        self.log_dir.mkdir(parents=True, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_file = self.log_dir / f"{timestamp}.log"
        self.latest_log = self.log_dir / "latest.log"

        self._logger = logging.getLogger(name)
        self._logger.setLevel(level)
        self._logger.handlers.clear()
        self._logger.propagate = False

        self._setup_handlers()

    def _setup_handlers(self) -> None:
        """Set up console and file handlers."""

        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(self.level)
        console_handler.setFormatter(ConsoleFormatter())
        self._logger.addHandler(console_handler)

        file_formatter = FileFormatter()

        timestamp_handler = logging.FileHandler(
            self.log_file,
            encoding="utf-8",
        )
        timestamp_handler.setLevel(logging.DEBUG)
        timestamp_handler.setFormatter(file_formatter)
        self._logger.addHandler(timestamp_handler)

        latest_handler = logging.FileHandler(
            self.latest_log,
            mode="w",
            encoding="utf-8",
        )
        latest_handler.setLevel(logging.DEBUG)
        latest_handler.setFormatter(file_formatter)
        self._logger.addHandler(latest_handler)

    # ------------------------------------------------------
    # Basic log methods
    # ------------------------------------------------------

    def debug(self, message: str) -> None:
        self._logger.debug(message)

    def info(self, message: str) -> None:
        self._logger.info(message)

    def success(self, message: str) -> None:
        self._logger.log(SUCCESS_LEVEL, message)

    def warning(self, message: str) -> None:
        self._logger.warning(message)

    def error(self, message: str) -> None:
        self._logger.error(message)

    def critical(self, message: str) -> None:
        self._logger.critical(message)

    def exception(self, message: str | BaseException) -> None:
        """Log an exception with traceback.

        Can be called inside an except block as:
            logger.exception(e)
        """

        if isinstance(message, BaseException):
            self._logger.error(str(message))
        else:
            self._logger.error(message)

        tb = traceback.format_exc()
        if tb and tb.strip() != "NoneType: None":
            self._logger.debug(tb)

    # ------------------------------------------------------
    # Section
    # ------------------------------------------------------

    @contextmanager
    def section(self, title: str) -> Iterator[None]:
        """Log a visible section block."""

        line = "=" * 60
        self.info("")
        self.info(line)
        self.info(title)
        self.info(line)

        try:
            yield
        except Exception as exc:
            self.exception(exc)
            raise
        finally:
            self.info("")

    # ------------------------------------------------------
    # Timer
    # ------------------------------------------------------

    @contextmanager
    def timer(self, name: str) -> Iterator[None]:
        """Measure and log runtime for a block."""

        start = time.perf_counter()
        self.info(f"{name} started")

        try:
            yield
        except Exception as exc:
            elapsed = time.perf_counter() - start
            self.timers[name] = elapsed
            self.error(f"{name} failed ({elapsed:.2f} sec)")
            self.exception(exc)
            raise
        else:
            elapsed = time.perf_counter() - start
            self.timers[name] = elapsed
            self.success(f"{name} finished ({elapsed:.2f} sec)")

    # ------------------------------------------------------
    # Summary
    # ------------------------------------------------------

    def summary(self) -> None:
        """Print runtime summary."""

        if not self.timers:
            self.info("No timed steps recorded.")
            return

        self.info("")
        self.info("=" * 60)
        self.info("Runtime Summary")
        self.info("=" * 60)

        total = 0.0

        for name, seconds in self.timers.items():
            self.info(f"{name:<30} {seconds:8.2f} sec")
            total += seconds

        self.success(f"Total runtime {total:.2f} sec")

    # ------------------------------------------------------
    # Utility
    # ------------------------------------------------------

    def get_log_file(self) -> Path:
        """Return timestamped log file path."""

        return self.log_file

    def get_latest_log(self) -> Path:
        """Return latest.log path."""

        return self.latest_log


logger = LoggerManager()