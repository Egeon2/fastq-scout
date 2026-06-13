from dataclasses import dataclass


@dataclass(frozen=True)
class RunContext:
    layout: str
    library_type: str
    expected_species: str | None = None

    @property
    def is_paired(self) -> bool:
        return self.layout == "paired"

    @property
    def is_transcriptome(self) -> bool:
        return self.library_type == "transcriptome"


def build_context(args) -> RunContext:
    return RunContext(
        layout=args.layout,
        library_type=args.library_type,
        expected_species=getattr(args, "expected_species", None),
    )
