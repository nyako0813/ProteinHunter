"""
Protein Hunter v5
Global constants

All URLs, filenames and application constants are defined here.
"""

from pathlib import Path

# ==========================================================
# Application
# ==========================================================

APP_NAME = "Protein Hunter"
APP_VERSION = "5.0"

# ==========================================================
# Python
# ==========================================================

MINIMUM_PYTHON = (3, 12)

# ==========================================================
# Directories
# ==========================================================

PROJECT_ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = PROJECT_ROOT / "data"

INPUT_DIR = DATA_DIR / "input"

OUTPUT_DIR = DATA_DIR / "output"

DATABASE_DIR = DATA_DIR / "databases"

CACHE_DIR = PROJECT_ROOT / ".cache"

LOG_DIR = PROJECT_ROOT / "logs"

# ==========================================================
# Cache
# ==========================================================

CDD_CACHE = CACHE_DIR / "cdd"

PFAM_CACHE = CACHE_DIR / "pfam"

UNIPROT_CACHE = CACHE_DIR / "uniprot"

ALPHAFOLD_CACHE = CACHE_DIR / "alphafold"

BLAST_CACHE = CACHE_DIR / "blast"

CACHE_DIRS = (
    CACHE_DIR,
    CDD_CACHE,
    PFAM_CACHE,
    UNIPROT_CACHE,
    ALPHAFOLD_CACHE,
    BLAST_CACHE,
)

# ==========================================================
# API
# ==========================================================

CDD_URL = (
    "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"
)

UNIPROT_URL = (
    "https://rest.uniprot.org/uniprotkb"
)

ALPHAFOLD_URL = (
    "https://alphafold.ebi.ac.uk/api/prediction"
)

# Pfam/HMMER
HMMER_SUBMIT_URL = (
    "https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan"
)

HMMER_RESULT_URL = (
    "https://www.ebi.ac.uk/Tools/hmmer/results"
)

# ==========================================================
# BLAST
# ==========================================================

BLASTP = "blastp"

MAKEBLASTDB = "makeblastdb"

BLAST_DB_EXTENSIONS = (
    ".pin",
    ".phr",
    ".psq",
)

# ==========================================================
# Network
# ==========================================================

DEFAULT_TIMEOUT = 60

DOWNLOAD_TIMEOUT = 180

API_RETRY = 3

API_SLEEP = 5

# ==========================================================
# Console
# ==========================================================

OK = "✓"

WARNING = "⚠"

ERROR = "✗"

INFO = "ℹ"

# ==========================================================
# Excel
# ==========================================================

RESULT1_SHEET = "result1"

RESULT2_SHEET = "result2"

SUMMARY_SHEET = "summary"

# ==========================================================
# Score
# ==========================================================

DEFAULT_SCORE = 0

# ==========================================================
# Logging
# ==========================================================

LOG_FORMAT = (
    "%(asctime)s | %(levelname)-8s | %(message)s"
)

DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

# ==========================================================
# Progress
# ==========================================================

TQDM_WIDTH = 100