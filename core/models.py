"""Simple data models used across ProteinHunter.

These dataclasses keep BLAST hits, annotation results, candidate scores, and
protein records in predictable Python structures that are easy to pass between
pipeline steps.
"""

from __future__ import annotations

from dataclasses import dataclass, field


MetadataValue = str | int | float | bool | None


@dataclass(slots=True)
class BlastHit:
    """A single BLAST match for a query protein."""

    query_id: str
    subject_id: str
    percent_identity: float
    alignment_length: int
    evalue: float
    bitscore: float
    source: str = "blast"


@dataclass(slots=True)
class DomainHit:
    """A domain or motif database hit found on a protein sequence."""

    source: str
    accession: str
    name: str
    description: str = ""
    evalue: float | None = None
    bitscore: float | None = None
    start: int | None = None
    end: int | None = None


@dataclass(slots=True)
class AnnotationResult:
    """Result from one annotation source for one protein."""

    protein_id: str
    source: str
    success: bool
    cached: bool = False
    error: str | None = None
    domains: list[DomainHit] = field(default_factory=list)
    motifs: list[str] = field(default_factory=list)
    metadata: dict[str, MetadataValue] = field(default_factory=dict)


@dataclass(slots=True)
class CandidateScore:
    """Score summary for a candidate protein."""

    protein_id: str
    total_score: float = 0.0
    components: dict[str, float] = field(default_factory=dict)
    reasons: list[str] = field(default_factory=list)

    def add_component(
        self,
        name: str,
        score: float,
        reason: str | None = None,
    ) -> None:
        """Add or replace a score component and refresh the total score."""
        self.components[name] = score
        self.total_score = sum(self.components.values())

        if reason:
            self.reasons.append(reason)


@dataclass(slots=True)
class ProteinRecord:
    """A protein sequence with its hits, annotations, score, and notes."""

    protein_id: str
    description: str = ""
    sequence: str = ""
    positive_hits: list[BlastHit] = field(default_factory=list)
    negative_hits: list[BlastHit] = field(default_factory=list)
    domains: list[DomainHit] = field(default_factory=list)
    motifs: list[str] = field(default_factory=list)
    annotations: dict[str, AnnotationResult] = field(default_factory=dict)
    score: CandidateScore | None = None
    alphafold_url: str | None = None
    uniprot_accession: str | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        """Return the number of amino acids in the protein sequence."""
        return len(self.sequence)

    @property
    def has_negative_hit(self) -> bool:
        """Return True when at least one negative BLAST hit is present."""
        return bool(self.negative_hits)


__all__: tuple[str, ...] = (
    "AnnotationResult",
    "BlastHit",
    "CandidateScore",
    "DomainHit",
    "MetadataValue",
    "ProteinRecord",
)
