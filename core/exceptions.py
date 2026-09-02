"""Custom exception classes used by ProteinHunter.

The exceptions in this module provide clear, beginner-friendly error messages
for configuration, startup checks, BLAST work, annotation, caching, scoring,
and output generation.
"""

from __future__ import annotations


class ProteinHunterError(Exception):
    """Base class for all custom ProteinHunter exceptions."""

    default_message: str = "ProteinHunter encountered an unexpected problem."

    def __init__(self, message: str | None = None) -> None:
        """Create an exception with a helpful message."""
        super().__init__(message or self.default_message)


class ConfigError(ProteinHunterError):
    """Raised when a configuration file or setting is invalid."""

    default_message = "Please check your configuration settings."


class StartupCheckError(ProteinHunterError):
    """Raised when an environment or startup check fails."""

    default_message = "ProteinHunter could not complete its startup checks."


class FileValidationError(ProteinHunterError):
    """Raised when an input file is missing, unreadable, or invalid."""

    default_message = "Please check that the input file exists and is valid."


class BlastError(ProteinHunterError):
    """Base class for BLAST-related errors."""

    default_message = "A BLAST step failed. Please check the BLAST setup and inputs."


class BlastDatabaseError(BlastError):
    """Raised when a BLAST database cannot be found or used."""

    default_message = "The BLAST database could not be found or opened."


class BlastExecutionError(BlastError):
    """Raised when a BLAST command fails while running."""

    default_message = "The BLAST command did not finish successfully."


class BlastParseError(BlastError):
    """Raised when BLAST results cannot be parsed."""

    default_message = "ProteinHunter could not read the BLAST results."


class AnnotationError(ProteinHunterError):
    """Base class for annotation-related errors."""

    default_message = "An annotation step failed. Please check the annotation inputs."


class CDDAnnotationError(AnnotationError):
    """Raised when CDD annotation fails."""

    default_message = "CDD annotation did not finish successfully."


class PfamAnnotationError(AnnotationError):
    """Raised when Pfam annotation fails."""

    default_message = "Pfam annotation did not finish successfully."


class UniProtAnnotationError(AnnotationError):
    """Raised when UniProt annotation fails."""

    default_message = "UniProt annotation did not finish successfully."


class AlphaFoldAnnotationError(AnnotationError):
    """Raised when AlphaFold annotation fails."""

    default_message = "AlphaFold annotation did not finish successfully."


class GeneContextError(AnnotationError):
    """Raised when gene neighborhood or context analysis fails."""

    default_message = "ProteinHunter could not analyze the gene context."


class StringPpiAnnotationError(AnnotationError):
    """Raised when STRING PPI evidence cannot be downloaded/parsed and no cache exists."""

    default_message = "STRING PPI evidence could not be retrieved."


class CacheError(ProteinHunterError):
    """Raised when cached data cannot be read or written."""

    default_message = "ProteinHunter could not use the cache."


class ScoringError(ProteinHunterError):
    """Raised when candidate scoring fails."""

    default_message = "ProteinHunter could not calculate candidate scores."


class ExcelOutputError(ProteinHunterError):
    """Raised when Excel output cannot be created."""

    default_message = "ProteinHunter could not create the Excel output file."


__all__: tuple[str, ...] = (
    "ProteinHunterError",
    "ConfigError",
    "StartupCheckError",
    "FileValidationError",
    "BlastError",
    "BlastDatabaseError",
    "BlastExecutionError",
    "BlastParseError",
    "AnnotationError",
    "CDDAnnotationError",
    "PfamAnnotationError",
    "UniProtAnnotationError",
    "AlphaFoldAnnotationError",
    "GeneContextError",
    "StringPpiAnnotationError",
    "CacheError",
    "ScoringError",
    "ExcelOutputError",
)
