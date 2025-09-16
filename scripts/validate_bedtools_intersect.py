#!/usr/bin/env python
"""Validation utilities for 'bedtools intersect -wo -s' output.

The input must be the result of:
* -a GFF3/GTF (annotated features, 1-based)
* -b BAM (alignments, 0-based)
* -wo (report only overlaps and include overlap length)
* -s (strand-aware intersections)

Expected column order (16 columns):
1-9 Feature (GFF3/GTF): chr, source, type, start, end, score, strand,
                        phase/frame, attributes
10-15 Read (BAM-derived): chr, start, end, name, score, strand
16 Overlap length (bp)

Exposes:
- FileFormatError: exception raised for malformed lines/fields.
- Record: a validated, parsed view of a single output line.
- validate_first_n: fail-fast validation of the first N lines of a file.
- parse_all: a streaming generator of '(line_number, Record)' tuples.
"""


from pathlib import Path
from typing import Dict, Literal, Iterator


class FileFormatError(BaseException):
    """Raised when a line or field does not match the expected format."""


# pylint: disable=too-many-locals
# pylint: disable=too-many-arguments
# pylint: disable=too-many-instance-attributes
# pylint: disable=too-many-positional-arguments
class Record:
    """A single validated 'bedtools intersect -wo -s' record.

    The expected column order is listed in the module docstring.

    Attributes:
        feat_chr: Feature chromosome of scaffold (with or without the 'chr'
            prefix).
        feat_source: Program or data source that produced the feature
            data source.
        feat_type: Feature type name.
        feat_start: Feature start (1-based).
        feat_end: Feature end (1-based).
        feat_score: Floating-point value or '.' if missing.
        feat_strand: Feature's strand defined as + (forward) or - (reverse).
        feat_phase_frame: '0', '1', or '2' indicating the feature's first base
            position within a codon, or '.' if missing.
        feat_attrs: Dict of parsed GFF3/GTF-style attributes (keys lowercased).

        read_chr: Read chromosome of scaffold (with or without the 'chr'
            prefix).
        read_start: Read start (0-based).
        read_end: Read end (0-based).
        read_name: Read name.
        read_score: Floating-point value or '.' if missing.
        read_strand: Read's strand defined as + (forward) or - (reverse).

        overlap_len: Number of overlapping base pairs between the feature and
            the read.
    """

    FIELDS = [
        "feat_chr",
        "feat_source",
        "feat_type",
        "feat_start",
        "feat_end",
        "feat_score",
        "feat_strand",
        "feat_phase_frame",
        "feat_attrs",
        "read_chr",
        "read_start",
        "read_end",
        "read_name",
        "read_score",
        "read_strand",
        "overlap_len",
    ]

    def __init__(
        self,
        feat_chr: str,
        feat_source: str,
        feat_type: str,
        feat_start: int,
        feat_end: int,
        feat_score: float | Literal["."],
        feat_strand: Literal["+", "-"],
        feat_phase_frame: int | Literal["."],
        feat_attrs: Dict[str, str],
        read_chr: str,
        read_start: int,
        read_end: int,
        read_name: str,
        read_score: float | Literal["."],
        read_strand: Literal["+", "-"],
        overlap_len: int,
    ) -> None:
        """Initialize a validated record.

        All arguments are expected to be already coerced to their target types.
        """
        self.feat_chr = feat_chr
        self.feat_source = feat_source
        self.feat_type = feat_type
        self.feat_start = feat_start
        self.feat_end = feat_end
        self.feat_score = feat_score
        self.feat_strand = feat_strand
        self.feat_phase_frame = feat_phase_frame
        self.feat_attrs = feat_attrs
        self.read_chr = read_chr
        self.read_start = read_start
        self.read_end = read_end
        self.read_name = read_name
        self.read_score = read_score
        self.read_strand = read_strand
        self.overlap_len = overlap_len

        self.cross_checks()

    @classmethod
    def from_line(cls, line: str, sep: str = "\t") -> "Record":
        """Parse and validate a single 'bedtools intersect -wo -s' line.

        Args:
            line: Raw line from the output file.
            sep: Field separator (default: TAB).

        Returns:
            A validated 'Record' instance.

        Raises:
            FileFormatError: If the column count is wrong or any field is
                invalid.
        """
        parts = line.strip().split(sep)

        if len(parts) != 16:
            raise FileFormatError(f"expected 16 columns, found {len(parts)}.")

        feat_chr: str = cls._text(parts[0].strip())
        feat_source: str = cls._text(parts[1].strip())
        feat_type: str = cls._text(parts[2].strip())

        feat_start: int = cls._as_int(
            text=parts[3].strip(), field="feature start (col. 4)"
        )

        feat_end: int = cls._as_int(
            text=parts[4].strip(), field="feature end (col. 5)"
        )

        feat_score: float | Literal["."] = cls._as_float_or_dot(
            text=parts[5].strip(), field="feature score (col. 6)"
        )

        feat_strand: Literal["+", "-"] = cls._as_strand(
            text=parts[6].strip(), field="feature strand (col. 7)"
        )

        feat_phase: int | Literal["."] = cls._as_phase_frame(
            text=parts[7].strip(), field="feature phase/frame (col. 8)"
        )

        feat_attrs: Dict[str, str] = cls._parse_feat_attrs(
            text=parts[8].strip(), field="feature attributes (col. 9)"
        )

        read_chr: str = cls._text(parts[9].strip())

        read_start: int = cls._as_int(
            text=parts[10].strip(), field="read start (col. 11)"
        )

        read_end: int = cls._as_int(
            text=parts[11].strip(), field="read end (col. 12)"
        )

        read_name: str = cls._text(parts[12].strip())

        read_score: float | Literal["."] = cls._as_float_or_dot(
            text=parts[13].strip(), field="read score (col. 14)"
        )

        read_strand: Literal["+", "-"] = cls._as_strand(
            text=parts[14].strip(), field="read strand (col. 15)"
        )

        overlap_len: int = cls._as_int(
            text=parts[15].strip(), field="overlap (col. 16)"
        )

        return cls(
            feat_chr=feat_chr,
            feat_source=feat_source,
            feat_type=feat_type,
            feat_start=feat_start,
            feat_end=feat_end,
            feat_score=feat_score,
            feat_strand=feat_strand,
            feat_phase_frame=feat_phase,
            feat_attrs=feat_attrs,
            read_chr=read_chr,
            read_start=read_start,
            read_end=read_end,
            read_name=read_name,
            read_score=read_score,
            read_strand=read_strand,
            overlap_len=overlap_len,
        )

    def cross_checks(self) -> None:
        """Validate relationships across fields.

        Ensures:
        - Feature and read are on the same chromosome and strand.
        - start <= end for both feature and read.
        - Overlap > 0 (required by `-wo`) and equals the computed overlap.

        Raises:
            FileFormatError: If any rule is violated.
        """
        if self.feat_chr != self.read_chr:
            raise FileFormatError(
                "chromosome mismatch: "
                f"feature (col. 1) is '{self.feat_chr}', read (col. 10) is "
                f"'{self.read_chr}'."
            )

        if self.feat_strand != self.read_strand:
            raise FileFormatError(
                f"strand mismatch: feature (col. 7) is '{self.feat_strand}' "
                f"and read (col. 15) is'{self.read_strand}'. Run bedtools "
                "with '-s' for strand-aware intersections."
            )

        if self.feat_start > self.feat_end:
            raise FileFormatError(
                f"feature range is invalid: start (col. 4) {self.feat_start} "
                f"> end (col. 5) {self.feat_end}."
            )

        if self.read_start > self.read_end:
            raise FileFormatError(
                f"read range is invalid: start (col. 11) {self.read_start} "
                f"> end (col. 12) {self.read_end}."
            )

        overlap = self._compute_overlap(
            f_start=self.feat_start,
            f_end=self.feat_end,
            r_start=self.read_start,
            r_end=self.read_end,
        )

        if overlap == 0:
            raise FileFormatError(
                "computed overlap is 0 between feature (cols. 4–5) and read "
                "(cols. 11–12). bedtools '-wo' outputs only overlaps—check "
                "coordinate conventions and inputs."
            )

        if self.overlap_len != overlap:
            raise FileFormatError(
                f"overlap length mismatch: expected overlap {overlap} but "
                f"column 16 is {self.overlap_len}."
            )

    @staticmethod
    def _text(text: str) -> str:
        """Return the text field unchanged."""
        return text

    @staticmethod
    def _as_int(text: str, field: str) -> int:
        """Parse an integer field.

        Args:
            text: Raw token.
            field: Human-friendly field name for error messages.

        Raises:
            FileFormatError: if 'text' cannot be parsed as an integer.
        """
        try:
            return int(text)
        except ValueError as err:
            raise FileFormatError(
                f"{field}: must be an integer, got {text!r}."
            ) from err

    @staticmethod
    def _as_float_or_dot(text: str, field: str) -> float | Literal["."]:
        """Parse a field that is either '.' or a float.

        Args:
            text: Raw token.
            field: Human-friendly field name for error messages.

        Raises:
            FileFormatError: if 'text' is neither '.' nor a valid float.
        """
        if text == ".":
            return "."
        try:
            return float(text)
        except ValueError as err:
            raise FileFormatError(
                f"{field}: must be a float or the string '.', got {text!r}."
            ) from err

    @staticmethod
    def _as_strand(text: str, field: str) -> Literal["+", "-"]:
        """Parse a strand field: '+' or '-' only.

        Args:
            text: Raw token.
            field: Human-friendly field name for error messages.

        Raises:
            FileFormatError: if 'text' is not '+' or '-'.
        """
        if text == "+":
            return "+"

        if text == "-":
            return "-"

        raise FileFormatError(
            f"{field}: strand must be '+' or '-', got {text!r}."
        )

    @staticmethod
    def _as_phase_frame(text: str, field: str) -> Literal["."] | int:
        """Parse a phase/frame field: 0, 1, 2 or '.'.

        Args:
            text: Raw token.
            field: Human-friendly field name for error messages.

        Raises:
            FileFormatError: if 'text' is not '.', 0, 1, or 2.
        """
        if text == ".":
            return "."

        if text in {"0", "1", "2"}:
            return int(text)

        raise FileFormatError(
            f"{field}: phase/frame must be 0, 1, 2 or '.', got {text!r}."
        )

    @staticmethod
    def _parse_feat_attrs(text: str, field: str) -> Dict[str, str]:
        """Parse a GFF3/GTF-style attributes field into a dict.

        Supports:
            - GFF3: key1=value1;key2=value2
            - GTF: key1 "value1"; key2 "value2";

        Keys are lowercased; surrounding quotes and spaces are stripped.
        Duplicated keys are no duplicated; later entries overwrite earlier
        ones.

        Args:
            text: Raw token.
            field: Human-friendly field name for error messages.

        Returns:
            Dict[str, str] of parsed attributes.

        Raises:
            FileFormatError: If the string is not one if the supported formats.
        """
        pairs = text.split(";")

        if len(pairs[0].split("=")) == 2:
            return {p.split("=")[0].lower(): p.split("=")[1] for p in pairs}

        if len(pairs[0].split('"')) == 3:
            return {
                p.split('"')[0].strip().lower(): p.split('"')[1]
                for p in filter(None, pairs)
            }

        raise FileFormatError(
            f'{field}: must be GFF3 (key=value;..) or GTF (key "value";...); '
            f"got {text!r}."
        )

    @staticmethod
    def _compute_overlap(
        f_start: int, f_end: int, r_start: int, r_end: int
    ) -> int:
        """Compute the number of overlapping bases between feature and read.

        Coordinates:
            - Feature: 1-based, closed [f_start, f_end]
            - Read: 0.based, half-open [r_start, r_end)

        The feature is converted to 0-based half-open before computing.

        Args:
            f_start: Feature start (1-based).
            f_end: Feature end (1-based).
            r_start: Read start (0-based).
            r_end: Read end (0-based).

        Returns:
            Overlap length in base pairs (0 if none).
        """
        low_coord = max(f_start - 1, r_start)
        high_coord = min(f_end, r_end)

        return max(0, high_coord - low_coord)


def validate_first_n(
    intersect_file: Path, n: int = 10, sep: str = "\t"
) -> None:
    """Validate the first 'n' lines of a 'bedtools intersect -wo -s' out file.

    Fails fast if any of the first 'n' lines is invalid.

    Args:
        intersect_file: Path to the intersect output file.
        n: Number of leading lines to validate.
        sep: Field separator (default: TAB).

    Raises:
        FileFormatError: If any of the first 'n' lines has an invalid format.
    """
    with open(intersect_file, "r", encoding="utf-8") as intersect:
        for i in range(n):
            line = intersect.readline()

            if not line:
                break

            try:
                Record.from_line(line=line.strip(), sep=sep)

            except FileFormatError as err:
                raise FileFormatError(
                    f"Invalid format in line {i + 1}: {err}"
                ) from err


def parse_all(
    intersect_file: Path, sep: str = "\t"
) -> Iterator[tuple[int, Record]]:
    """Stream a file, yielding '(line_number, Record)' for each data line.

    Args:
        intersect_file: Path to the intersect output file.
        sep: Field separator (default: TAB).

    Yields:
        Tuples of '(line_number, Record)'.

    Raises:
        FileFormatError: If any line has an invalid format (iteration stops at
            the first error).
    """
    with open(intersect_file, "r", encoding="utf-8") as intersect:
        for line_num, line in enumerate(intersect, 1):
            try:
                yield line_num, Record.from_line(line=line.strip(), sep=sep)

            except FileFormatError as err:
                raise FileFormatError(
                    f"Invalid format in line {line_num}: {err}"
                ) from err
