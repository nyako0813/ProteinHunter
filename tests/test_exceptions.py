"""Tests for ProteinHunter custom exceptions."""

from __future__ import annotations

from core.exceptions import (
    __all__,
    AlphaFoldAnnotationError,
    AnnotationError,
    BlastDatabaseError,
    BlastError,
    BlastExecutionError,
    BlastParseError,
    CDDAnnotationError,
    CacheError,
    ConfigError,
    ExcelOutputError,
    FileValidationError,
    GeneContextError,
    PfamAnnotationError,
    ProteinHunterError,
    ScoringError,
    StartupCheckError,
    UniProtAnnotationError,
)


CUSTOM_EXCEPTIONS: tuple[type[ProteinHunterError], ...] = (
    ConfigError,
    StartupCheckError,
    FileValidationError,
    BlastError,
    BlastDatabaseError,
    BlastExecutionError,
    BlastParseError,
    AnnotationError,
    CDDAnnotationError,
    PfamAnnotationError,
    UniProtAnnotationError,
    AlphaFoldAnnotationError,
    GeneContextError,
    CacheError,
    ScoringError,
    ExcelOutputError,
)

BLAST_EXCEPTIONS: tuple[type[BlastError], ...] = (
    BlastDatabaseError,
    BlastExecutionError,
    BlastParseError,
)

ANNOTATION_EXCEPTIONS: tuple[type[AnnotationError], ...] = (
    CDDAnnotationError,
    PfamAnnotationError,
    UniProtAnnotationError,
    AlphaFoldAnnotationError,
    GeneContextError,
)


def test_all_custom_exceptions_inherit_from_base_error() -> None:
    """Every custom exception should inherit from ProteinHunterError."""
    for exception_class in CUSTOM_EXCEPTIONS:
        assert issubclass(exception_class, ProteinHunterError)


def test_blast_specific_exceptions_inherit_from_blast_error() -> None:
    """BLAST-specific exceptions should inherit from BlastError."""
    for exception_class in BLAST_EXCEPTIONS:
        assert issubclass(exception_class, BlastError)


def test_annotation_specific_exceptions_inherit_from_annotation_error() -> None:
    """Annotation-specific exceptions should inherit from AnnotationError."""
    for exception_class in ANNOTATION_EXCEPTIONS:
        assert issubclass(exception_class, AnnotationError)


def test_custom_messages_are_preserved() -> None:
    """Custom messages should be returned unchanged."""
    message = "The FASTA file path is missing."

    error = FileValidationError(message)

    assert str(error) == message


def test_default_messages_are_usable() -> None:
    """Default messages should be non-empty and beginner-friendly."""
    for exception_class in CUSTOM_EXCEPTIONS:
        message = str(exception_class())
        assert message
        assert "None" not in message


def test_all_exports_include_exception_classes() -> None:
    """The public export list should include all exception class names."""
    expected_names = {
        "ProteinHunterError",
        *(exception_class.__name__ for exception_class in CUSTOM_EXCEPTIONS),
    }

    assert set(__all__) == expected_names
