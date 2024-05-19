#!/usr/bin/env python
"""Extend miRNA start and end coordinates."""

import argparse
from pathlib import Path
from typing import Optional

import gffutils  # type: ignore


class MirnaExtension:
    """Class for extending miRNA start and end coordinates.

    Attributes:
        db: In-memory database of the input GFF3 records.
        db: In-memory database of the updated GFF3 records.
        seq_lengths: Dictionary mapping reference sequence IDs to their lengths.
    """

    def __init__(self) -> None:
        """Initialize class."""
        self.db: Optional[gffutils.FeatureDB] = None
        self.db_out: Optional[gffutils.FeatureDB] = None
        self.seq_lengths: Optional[dict[str, int]] = None

    def set_db(self, path: Path) -> None:
        """Load GFF3 file into `gffutils.FeatureDB`.

        Args:
            gff_file: Path to a GFF3 file.
        """
        try:
            self.db = gffutils.create_db(
                str(path),
                dbfn=":memory:",
                force=True,
                keep_order=True
            )
        except gffutils.exceptions.EmptyInputError:
            pass

    def set_seq_lengths(self, path: Optional[Path] = None) -> None:
        """Set the reference sequence lengths.

        Args:
            path: Path to a tabulated file containing the names of reference sequences
                and the corresponding total length, in bases/basepairs.
        """
        self.seq_lengths = {}
        if path is None:
            for _id in self.db.seqids():  # type: ignore
                self.seq_lengths[_id] = max(
                    record.end for record in self.db.region(_id)  # type: ignore
                )
        else:
            with open(path, encoding="utf-8") as _file:
                for line in _file:
                    ref_seq, length = line.strip().split("\t")
                    try:
                        self.seq_lengths[ref_seq] = int(length)
                    except ValueError as exc:
                        raise ValueError(
                            f'Invalid length: "{length}" for sequence "{ref_seq}";'
                            ' integer value expected'
                        ) from exc

    def process_precursor(
        self,
        precursor: gffutils.Feature,
        n: int = 6
    ) -> list[gffutils.Feature]:
        assert isinstance(self.seq_lengths, dict)

        pre_orig_start = precursor.start
        pre_orig_end = precursor.end

        matures = list(
            self.db.region(  # type: ignore
                strand=precursor.strand,
                seqid=precursor.seqid,
                start=precursor.start,
                end=precursor.end,
                featuretype="miRNA",
                completely_within=True,
            )
        )
        for mir in matures:
            try:
                if self.seq_lengths[mir.seqid] < mir.end:
                    raise ValueError(
                        "The provided GFF3 miRNA annotations and reference sequence"
                        f" lengths are incompatible: end coordinate {mir.end} of miRNA"
                        f' "{mir.id}" exceeds length of reference sequence '
                        f' "{mir.seqid}" ({self.seq_lengths[mir.seqid]} nt)'
                    )
            except KeyError as exc:
                raise KeyError(
                    "The provided GFF3 miRNA annotations and reference sequence lengths"
                    f' are incompatible: reference sequence "{mir.seqid}" of miRNA'
                    f' "{mir.id}" is not available in the provided reference sequence'
                    " lengths table"
                ) from exc
            mir.start = max(mir.start - n, 0)
            mir.end = min(mir.end + n, self.seq_lengths[mir.seqid])

        precursor.start = min(precursor.start, min(mir.start for mir in matures))
        precursor.end = max(precursor.end, max(mir.end for mir in matures))
        precursor.attributes["Name"][0] += f"_-{pre_orig_start - precursor.start}"
        precursor.attributes["Name"][0] += f"_+{precursor.end - pre_orig_end}"

        return [precursor] + matures

    def extend_mirnas(
        self,
        n: int = 6,
    ) -> None:
        """Extend miRNAs start and end coordinates.

        This method elongates the start and end coordinates of mature miRNAs
        by n nucleotides. In the case that this extension makes the start/end
        coordinates to exceed the corresponding primary miRNA boundaries,
        these will be extended as far as the miRNA coordinates.
        If provided, the elongation will take into account the chromosome size.

        Args:
            n: Number of nucleotides to extend miRs start and end coordinates.
        """
        assert isinstance(self.db, gffutils.FeatureDB)

        feats_updated: list[gffutils.Feature] = []

        for pre in self.db.features_of_type(  # type: ignore
            "miRNA_primary_transcript"
        ):
            feats_updated.extend(self.process_precursor(precursor=pre, n=n))

        self.db_out = gffutils.create_db(
            feats_updated,
            dbfn=":memory:",
            force=True,
            keep_order=True,
        )

    def write_gff(self, path: Path, feature_type: Optional[str] = None) -> None:
        """Write features to a GFF3 file.

        Args:
            path: Path to the output file.
            feature_type: Feature type to write. If `None`, all features will be
                written.
        """
        with open(path, "w", encoding="utf-8") as _file:
            if not isinstance(self.db_out, gffutils.FeatureDB):
                _file.write("")
            elif feature_type is None:
                for feature in self.db_out.all_features():  # type: ignore
                    _file.write(str(feature) + "\n")
            else:
                for feature in self.db_out.features_of_type(  # type: ignore
                    feature_type
                ):
                    _file.write(str(feature) + "\n")


def parse_arguments():
    """Parse command-line arguments."""
    description = """Extend miRNA start and end coordinates.

Extends mature miRNA start and end coordinates in a GFF3 file by the indicated number of
nucleotides, but not exceeding chromosome boundaries. If the extension causes the mature
miRNA coordinates to exceed the boundaries of the corresponding precursor (see note in
"Constraints" below), the precursor coordinates will be extended accordingly. In this
case, the precursor name will be appended with "_-n" and/or "_+m", where n and m
represent the extensions on the 5' and 3' end, respectively. The modified mature miRNA
and precursor annotations will be written to separate GFF3 files.

Constraints:
The script was written for GFF3 files obtained from miRBase (http://www.mirbase.org/).
The following assumptions are made:
- The input GFF3 file contains miRNA annotations for a single species.
- The input GFF3 file contains features of type "miRNA_primary_transcript" (referred to
  as "precursors") and "miRNA" (referred to as "mature miRNAs").
- The input GFF3 file contains a Name attribute for each precursor.
- Each precursor contains one or more mature miRNAs.
- Each mature miRNA is a child of exactly one precursor and is completely within the
  boundaries of the precursor.
"""
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 1.1.0",
        help="show program's version number and exit",
    )
    parser.add_argument(
        "input",
        help="path to the GFF3 annotation file",
        type=Path,
    )
    parser.add_argument(
        "-e",
        "--extension",
        help=(
            "number of nucleotides by which to extend the miRNA boundaries; default:"
            " %(default)d"
        ),
        default=6,
        type=int,
        metavar="INT",
    )
    parser.add_argument(
        "--chr",
        help=(
            "path to a tabulated file containing the chromosome names and their"
            " lengths, in basepairs; if not provided, the sequence lengths will be"
            " set to the most distal 3' boundary of the miRNA precursors annotated"
            " for each chromosome in the input GFF3 file; default: %(default)s"
        ),
        default=None,
        type=Path,
        metavar="PATH",
    )
    parser.add_argument(
        "--outdir",
        help="path to the output directory; default: %(default)s",
        type=Path,
        default=Path.cwd(),
        metavar="PATH",
    )
    return parser


def main(args) -> None:
    """Extend miRNAs start/end coordinates."""
    args.outdir.mkdir(parents=True, exist_ok=True)
    mirna_extension = MirnaExtension()
    mirna_extension.set_db(path=args.input)
    if isinstance(mirna_extension.db, gffutils.FeatureDB):
        mirna_extension.set_seq_lengths(path=args.chr)
        mirna_extension.extend_mirnas(n=args.extension)
    mirna_extension.write_gff(
        path=args.outdir / f"extended_primir_annotation_{args.extension}_nt.gff3",
        feature_type="miRNA_primary_transcript",
    )
    mirna_extension.write_gff(
        path=args.outdir / f"extended_mirna_annotation_{args.extension}_nt.gff3",
        feature_type="miRNA",
    )


if __name__ == "__main__":
    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
