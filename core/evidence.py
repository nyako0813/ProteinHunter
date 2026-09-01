"""Evidence model shared by the interaction scoring engine.

This module keeps three things separate, on purpose (see
``docs/integrated_scoring_design.md`` and the project design specification):

1. ``raw_value``      -- the untouched measurement (a BLAST bitscore, a
   distance in bp, a set of shared keywords, ...).
2. ``normalized_value`` -- the same signal mapped onto ``0.0``-``1.0`` so
   different evidence types can be combined.
3. ``status``          -- whether the evidence was actually evaluated at
   all. A component that was never evaluated (no GFF coordinates, no
   annotation text, an annotation step that was skipped) must never be
   silently treated as "evaluated, score 0". Mixing those two cases was the
   main accuracy problem in the old additive scorer.

Nothing here imports pandas, openpyxl, or any I/O code, so it can be unit
tested in isolation and reused by any future scoring model.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class EvidenceStatus(str, Enum):
    """Lifecycle state of one evidence component.

    Only :attr:`AVAILABLE` contributes to a score. Every other status is
    excluded from the scoring denominator -- it is never converted into a
    zero or a negative value. This mirrors the "missing is not negative"
    rule used by both ProteinHunter's design goals and
    ProteinInteractionHunter's evidence model.
    """

    AVAILABLE = "AVAILABLE"
    """The component was evaluated and produced a usable value (which may
    legitimately be 0.0, e.g. "evaluated, no shared keyword found")."""

    MISSING = "MISSING"
    """The required input data was absent (no annotation, no sequence)."""

    NOT_RUN = "NOT_RUN"
    """The evidence engine was disabled in configuration."""

    NOT_APPLICABLE = "NOT_APPLICABLE"
    """The comparison does not make sense for this pair (e.g. different
    contigs for a genomic-distance measurement)."""

    FAILED = "FAILED"
    """The evidence engine raised an error while evaluating this pair."""

    MALFORMED = "MALFORMED"
    """The input needed for this component was present but invalid."""

    EXCLUDED = "EXCLUDED"
    """The component was intentionally left out of this run (e.g. by
    configuration or a deduplication rule)."""


#: Statuses that never contribute weight to a score denominator.
NON_CONTRIBUTING_STATUSES: frozenset[EvidenceStatus] = frozenset(
    {
        EvidenceStatus.MISSING,
        EvidenceStatus.NOT_RUN,
        EvidenceStatus.NOT_APPLICABLE,
        EvidenceStatus.FAILED,
        EvidenceStatus.MALFORMED,
        EvidenceStatus.EXCLUDED,
    }
)


def clamp01(value: float) -> float:
    """Clamp a float into the closed ``[0.0, 1.0]`` range."""
    if value < 0.0:
        return 0.0
    if value > 1.0:
        return 1.0
    return value


def linear_normalize(value: float, low: float, high: float) -> float:
    """Map ``value`` linearly onto ``[0.0, 1.0]`` given ``[low, high]``.

    ``value <= low`` maps to ``0.0`` and ``value >= high`` maps to ``1.0``.
    ``low`` may be greater than ``high`` to invert the scale (for example,
    turning "smaller distance is better" into a 0-1 support score).
    """
    if low == high:
        raise ValueError("linear_normalize requires low != high")
    fraction = (value - low) / (high - low)
    return clamp01(fraction)


@dataclass(slots=True)
class EvidenceComponent:
    """One normalized, auditable evidence signal for a candidate pair.

    ``category`` groups components that measure correlated signals (see
    ``analysis/scoring_engine.py`` for how categories are capped so that
    correlated evidence cannot be double counted).
    """

    name: str
    category: str
    status: EvidenceStatus
    raw_value: object = None
    normalized_value: float | None = None
    weight: float = 0.0
    is_negative: bool = False
    source: str = ""
    explanation: str = ""

    def __post_init__(self) -> None:
        if self.status is EvidenceStatus.AVAILABLE:
            if self.normalized_value is None:
                raise ValueError(
                    f"evidence component '{self.name}' is AVAILABLE but has "
                    "no normalized_value"
                )
            if not 0.0 <= self.normalized_value <= 1.0:
                raise ValueError(
                    f"evidence component '{self.name}' normalized_value "
                    f"must be within [0.0, 1.0], got {self.normalized_value!r}"
                )
            if self.weight < 0.0:
                raise ValueError(
                    f"evidence component '{self.name}' weight must be >= 0, "
                    f"got {self.weight!r}"
                )

    @property
    def effective_weight(self) -> float:
        """Weight actually usable in a denominator (0 unless AVAILABLE)."""
        if self.status is not EvidenceStatus.AVAILABLE:
            return 0.0
        return self.weight

    @property
    def contribution(self) -> float:
        """``effective_weight * normalized_value``, 0 when unavailable."""
        if self.status is not EvidenceStatus.AVAILABLE or self.normalized_value is None:
            return 0.0
        return self.effective_weight * self.normalized_value

    @classmethod
    def available(
        cls,
        name: str,
        category: str,
        normalized_value: float,
        weight: float,
        *,
        raw_value: object = None,
        is_negative: bool = False,
        source: str = "",
        explanation: str = "",
    ) -> "EvidenceComponent":
        """Convenience constructor for an evaluated component."""
        return cls(
            name=name,
            category=category,
            status=EvidenceStatus.AVAILABLE,
            raw_value=raw_value,
            normalized_value=clamp01(normalized_value),
            weight=weight,
            is_negative=is_negative,
            source=source,
            explanation=explanation,
        )

    @classmethod
    def unavailable(
        cls,
        name: str,
        category: str,
        status: EvidenceStatus,
        *,
        raw_value: object = None,
        source: str = "",
        explanation: str = "",
    ) -> "EvidenceComponent":
        """Convenience constructor for a component that was not evaluated."""
        if status is EvidenceStatus.AVAILABLE:
            raise ValueError("use EvidenceComponent.available() for AVAILABLE status")
        return cls(
            name=name,
            category=category,
            status=status,
            raw_value=raw_value,
            normalized_value=None,
            weight=0.0,
            source=source,
            explanation=explanation,
        )


__all__: tuple[str, ...] = (
    "EvidenceComponent",
    "EvidenceStatus",
    "NON_CONTRIBUTING_STATUSES",
    "clamp01",
    "linear_normalize",
)
