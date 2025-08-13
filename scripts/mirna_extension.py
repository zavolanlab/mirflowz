#!/usr/bin/env python
"""Extend miRNA start and end coordinates and ensure name uniqueness."""

import argparse
from pathlib import Path
from typing import Optional

import gffutils  # type: ignore


class AnnotationException(BaseException):
    """A custom exception class for MirnaExtension class."""


class MirnaExtension:
    """Class for updating miRNA annotated coordinates and names.

    Attributes:
        db: In-memory database of the input GFF3 records.
        db_out: In-memory database of the updated GFF3 records.
        seq_lengths: Dictionary mapping reference sequence IDs to their
            lengths.
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
                str(path), dbfn=":memory:", force=True, keep_order=True
            )
        except gffutils.exceptions.EmptyInputError:
            pass
        except AttributeError as err:
            raise AnnotationException(
                "\n\n"
                "Illegal coordinate: The provided GFF3 miRNA annotation file"
                " contains at least one miRNA species starting at position 0."
                " Please, check that all entries start at least at position 1"
                "\n"
            ) from err

    def set_seq_lengths(self, path: Optional[Path] = None) -> None:
        """Set the reference sequence lengths.

        Args:
            path: Path to a tabulated file containing the names of reference
                sequences and the corresponding total length, in base pairs.
        """
        self.seq_lengths = {}
        anno_lengths = {}

        if not isinstance(self.db, gffutils.FeatureDB):
            return

        for _id in self.db.seqids():  # type: ignore
            anno_lengths[_id] = max(
                rec.end for rec in self.db.region(_id)  # type: ignore
            )

        if path is None:
            self.seq_lengths = anno_lengths
        else:
            with open(path, encoding="utf-8") as _file:
                for line in _file:
                    ref_seq, length = line.strip().split("\t")

                    try:
                        max_len = anno_lengths[ref_seq]
                    except KeyError:
                        max_len = 0

                    try:
                        self.seq_lengths[ref_seq] = int(length)
                    except ValueError as exc:
                        raise ValueError(
                            f'Invalid length: "{length}" for sequence'
                            f' "{ref_seq}"; integer value expected'
                        ) from exc

                    if max_len > int(length):
                        raise AnnotationException(
                            "The provided GFF3 miRNA annotations and"
                            " reference sequence lengths are incompatible:"
                            f" end coordinate {max_len} exceeds length of"
                            f' reference sequence "{ref_seq}" ({length} nt)'
                        )

    def adjust_names(
        self, precursor: gffutils.Feature, matures: list[gffutils.Feature]
    ) -> None:
        """Adjust miRNA attributes for uniqueness and consistency.

        This method updates the attributes of precursor and mature miRNA
        entries to ensure consistent naming based on their paralog or sequence
        variant status.

        For precursors:
            - A suffix indicates distinct genomic loci (paralogs) that express
              identical mature sequences. This is typically extracted from the
              'Name' or 'ID' attribute.
            - Format: 'SPECIES-mir-NUMBER[LETTER]-#' (Name) and 'ALIAS_#' (ID)
              where:
                  - 'LETTER' denotes a sequence variant of the mature miRNA
                    (paralogous variant with similar but not identical
                    sequences),
                  - '#' indicates the paralog number (replica/locus index),
                    included when multiple loci express the same or similar
                    miRNAs.

        For mature miRNAs:
            - The replica number is added or replaced as an infix/suffix in
              the name.
            - Formats:
                - 'SPECIES-miR-NUMBER[LETTER]-#-ARM'
                - 'SPECIES-miR-NUMBER[LETTER]-#'
                - 'SPECIES-miR-NUMBER[LETTER]-ARM'

        Cases:
            - If a precursor has multiple genomic instances (paralogs), the
              first occurrence typically lacks a numeric suffix; subsequent
              ones are numbered incrementally.
            - The 'Derives_from' attribute of each mature miRNA is updated to
              match the precursor's 'ID'.
            - If a precursor has a different sequence variant designation
              ('LETTER') than its associated matures, the mature miRNA names
              are updated to match the precursor's designation.

        The 'Alias' attribute remains unchanged.

        Args:
            precursor: 'miRNA primary transcript' feature entry
            matures: list of the corresponding 'mature miRNA' feature(s)
                    entry(s)
        """
        precursor_name_parts = precursor.attributes["Name"][0].split("-")
        precursor_id_parts = precursor.attributes["ID"][0].split("_")

        replica = None
        if len(precursor_name_parts) == 4:
            replica = precursor_name_parts[3]

        elif len(precursor_name_parts) == 3 and len(precursor_id_parts) == 2:
            replica = precursor_id_parts[1]
            precursor_name_parts.append(replica)
            precursor.attributes["Name"][0] = "-".join(precursor_name_parts)

        for mature in matures:
            mature_name_parts = mature.attributes["Name"][0].split("-")
            mature.attributes["Derives_from"][0] = "_".join(precursor_id_parts)

            if replica is not None:
                # Format: SPECIES-miR-NAME%-#-ARM
                if len(mature_name_parts) == 5:
                    mature_name_parts[3] = replica

                # Formats:
                #   - SPECIES-miR-NAME%-#
                #   - SPECIES-miR-NAME%-ARM
                elif len(mature_name_parts) == 4:
                    if mature_name_parts[3].isdigit():
                        mature_name_parts[3] = replica
                    else:
                        mature_name_parts.insert(3, replica)

                # Format: SPECIES-miR-NAME%
                elif len(mature_name_parts) == 3:
                    mature_name_parts.append(replica)
            else:
                if (
                    len(mature_name_parts) >= 3
                    and mature_name_parts[2] != precursor_name_parts[2]
                ):
                    mature_name_parts[2] = precursor_name_parts[2]

            mature.attributes["Name"][0] = "-".join(mature_name_parts)

    def process_precursor(
        self, precursor: gffutils.Feature, n: int = 6
    ) -> list[gffutils.Feature]:
        """Extend miRNAs start and end coordinates and ensure name uniqueness.

        This method elongates the start and end coordinates of mature miRNAs
        by 'n' nucleotides. In the case that this extension makes the start/end
        coordinates to exceed the corresponding primary miRNA boundaries,
        these will be extended as far as the miRNA coordinates.
        If provided, the elongation will take into account the chromosome size.

        In addition, the method `adjust_names` is called to ensure uniqueness
        in the `Name` attribute for both the precursor and its arms.

        Args:
            precursor: 'miRNA primary transcript' feature entry
            n: Number of nucleotides to extend miRs start and end coordinates.
        """
        assert isinstance(self.seq_lengths, dict)

        pre_start = precursor.start
        pre_end = precursor.end

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

        self.adjust_names(precursor=precursor, matures=matures)

        for mir in matures:
            try:
                mir.start = max(mir.start - n, 1)
                mir.end = min(mir.end + n, self.seq_lengths[mir.seqid])
            except KeyError as exc:
                raise KeyError(
                    "The provided GFF3 miRNA annotations and reference"
                    "sequence lengths are incompatible: reference sequence"
                    f' "{mir.seqid}" of miRNA "{mir.id}" is not available in'
                    " the provided reference sequence lengths table"
                ) from exc

        precursor.start = min(
            precursor.start, min(mir.start for mir in matures)
        )
        precursor.end = max(precursor.end, max(mir.end for mir in matures))
        precursor.attributes["Name"][0] += f"_-{pre_start - precursor.start}"
        precursor.attributes["Name"][0] += f"_+{precursor.end - pre_end}"

        return [precursor] + matures

    def update_db(
        self,
        n: int = 6,
    ) -> None:
        """Update miRNA annotations in the local database.

        Using the method `process_precursor` annotated coordinates are
        extended by `n` nucleotides and, if needed, the `Name` attribute is
        modified to contain the sequence replica number.

        Args:
            n: Number of nucleotides to extend miRs start and end coordinates.
        """
        if not isinstance(self.db, gffutils.FeatureDB):
            return

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

    def write_gff(
        self, path: Path, feature_type: Optional[str] = None
    ) -> None:
        """Write features to a GFF3 file.

        Args:
            path: Path to the output file.
            feature_type: Feature type to write. If `None`, all features will
                be written.
        """
        with open(path, "w", encoding="utf-8") as _file:
            if not isinstance(self.db_out, gffutils.FeatureDB):
                _file.write("")
            elif feature_type is None:
                for feature in self.db_out.all_features(  # type: ignore
                    order_by=("featuretype", "strand", "start", "end")
                ):
                    _file.write(str(feature) + "\n")
            else:
                for feature in self.db_out.features_of_type(  # type: ignore
                    feature_type, order_by=("strand", "start", "end")
                ):
                    _file.write(str(feature) + "\n")


def parse_arguments():
    """Parse command-line arguments."""
    description = """Extend miRNA annotated regions and ensure name uniqueness.

This method updates the attributes of precursor and mature miRNA
entries to ensure consistent naming based on their paralog or sequence
variant status.

For precursors:
    - A suffix indicates distinct genomic loci (paralogs) that express
      identical mature sequences. This is typically extracted from the 'Name'
      or 'ID' attribute.
    - Format: 'SPECIES-mir-NUMBER[LETTER]-#' (Name) and 'ALIAS_#' (ID)
      where:
          - 'LETTER' denotes a sequence variant of the mature miRNA (paralogous
            variant with similar but not identical sequences),
          - '#' indicates the paralog number (replica/locus index), included
            when multiple loci express the same or similar miRNAs.

For mature miRNAs:
    - The replica number is added or replaced as an infix/suffix in the name.
    - Formats:
        - 'SPECIES-miR-NUMBER[LETTER]-#-ARM'
        - 'SPECIES-miR-NUMBER[LETTER]-#'
        - 'SPECIES-miR-NUMBER[LETTER]-ARM'

Cases:
    - If a precursor has multiple genomic instances (paralogs), the first
      occurrence typically lacks a numeric suffix; subsequent ones are numbered
      incrementally.
    - The 'Derives_from' attribute of each mature miRNA is updated to match the
      precursor's 'ID'.
    - If a precursor has a different sequence variant designation ('LETTER')
      than its associated matures, the mature miRNA names are updated to match
      the precursor's designation.

Note that the 'Alias' attribute remains unchanged so repeated values may still
be present.

Extend mature miRNA start and end coordinates in a GFF3 file by the indicated
number of nucleotides without exceeding chromosome boundaries. If the
extension causes the mature miRNA coordinates to exceed the boundaries of the
corresponding precursor (see note in "Constraints" below), the precursor
coordinates are extended accordingly. The precursor name will be appended
with '_-n' and/or '_+m', where n and m represent the extensions on the 5' and
3' end, respectively (or 0 otherwise). The modified mature miRNA and precursor
annotations will be written to separate GFF3 files.

Constraints:
The script was written for GFF3 files obtained from miRBase
(http://www.mirbase.org/).

The following assumptions are made:
- The input GFF3 file contains miRNA annotations for a single species.
- The input GFF3 file contains features of type "miRNA_primary_transcript"
  (referred to as "precursors") and "miRNA" (referred to as "mature miRNAs").
- The input GFF3 file contains a 'Name' and 'ID' attribute for each precursor.
- Each precursor contains one or more mature miRNAs.
- Each mature miRNA is a child of exactly one precursor and is completely
  within the boundaries of the precursor.
"""
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 1.2.0",
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
            "number of nucleotides by which to extend the miRNA boundaries;"
            " default: %(default)d"
        ),
        default=6,
        type=int,
        metavar="INT",
    )
    parser.add_argument(
        "--chr",
        help=(
            "path to a tabulated file containing the chromosome names and"
            " their lengths, in basepairs; if not provided, the sequence"
            " lengths will be set to the most distal 3' boundary of the miRNA"
            " precursors annotated for each chromosome in the input GFF3 file;"
            " default: %(default)s"
        ),
        default=None,
        type=Path,
        metavar="PATH",
    )
    parser.add_argument(
        "--outdir",
        help="path to the output directory; default: %(default)s",
        default=Path.cwd(),
        type=Path,
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
        mirna_extension.update_db(n=args.extension)

    mirna_extension.write_gff(
        path=args.outdir
        / f"extended_primir_annotation_{args.extension}_nt.gff3",
        feature_type="miRNA_primary_transcript",
    )
    mirna_extension.write_gff(
        path=args.outdir
        / f"extended_mirna_annotation_{args.extension}_nt.gff3",
        feature_type="miRNA",
    )


if __name__ == "__main__":
    arguments = parse_arguments().parse_args()  # pragma: no cover
    main(arguments)  # pragma: no cover
