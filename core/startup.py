"""
Protein Hunter v5
Startup Checker

Checks environment before analysis starts.
"""

from __future__ import annotations

import importlib
import shutil
import socket
import sys
from pathlib import Path

from colorama import Fore, Style, init

init(autoreset=True)


class StartupChecker:

    def __init__(self):

        self.failed = False

    # ==========================================================
    # Internal
    # ==========================================================

    def _ok(self, message):

        print(Fore.GREEN + "✓ " + message)

    def _warn(self, message):

        print(Fore.YELLOW + "⚠ " + message)

    def _error(self, message):

        self.failed = True
        print(Fore.RED + "✗ " + message)

    # ==========================================================
    # Banner
    # ==========================================================

    def banner(self):

        print("=" * 60)
        print(" Protein Hunter v5")
        print("=" * 60)
        print()

    # ==========================================================
    # Python
    # ==========================================================

    def check_python(self):

        version = sys.version_info

        if version >= (3, 12):
            self._ok(f"Python {version.major}.{version.minor}.{version.micro}")

        else:
            self._error(
                f"Python 3.12+ required "
                f"(current {version.major}.{version.minor}.{version.micro})"
            )

    # ==========================================================
    # Packages
    # ==========================================================

    def check_packages(self):

        required = [
            "Bio",
            "pandas",
            "requests",
            "yaml",
            "openpyxl",
            "numpy",
            "tqdm",
        ]

        print()

        for pkg in required:

            try:

                importlib.import_module(pkg)

                self._ok(pkg)

            except Exception:

                self._error(pkg)

    # ==========================================================
    # Config
    # ==========================================================

    def check_config(self):

        print()

        cfg = Path("config.yaml")

        if cfg.exists():
            self._ok("config.yaml")

        else:
            self._error("config.yaml not found")

    # ==========================================================
    # FASTA
    # ==========================================================

    def check_file(self, path: Path, name: str):

        if path.exists():

            self._ok(name)

        else:

            self._error(f"{name} not found")

    # ==========================================================
    # BLAST
    # ==========================================================

    def check_blast(self):

        print()

        for exe in ("blastp", "makeblastdb"):

            if shutil.which(exe):

                self._ok(exe)

            else:

                self._error(exe)

    # ==========================================================
    # Internet
    # ==========================================================

    def check_internet(self):

        print()

        try:

            socket.create_connection(("8.8.8.8", 53), timeout=3)

            self._ok("Internet")

        except OSError:

            self._warn("Internet unavailable")

    # ==========================================================
    # Summary
    # ==========================================================

    def summary(self):

        print()
        print("-" * 60)

        if self.failed:

            print(Fore.RED + "Startup check FAILED")

        else:

            print(Fore.GREEN + "Startup check PASSED")

        print("-" * 60)
        print()

    # ==========================================================
    # Main
    # ==========================================================

    def run(self):

        self.banner()

        self.check_python()

        self.check_packages()

        self.check_config()

        self.check_blast()

        self.check_internet()

        self.summary()

        return not self.failed