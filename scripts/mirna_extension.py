#!/usr/bin/env python

"""Extend miRNAs start and end coordinates by n nucleotides.

This script uses the class MirnaExtension to extend the start and end
coordiantes of the mature miRNAs.This class is build upon two methods. The
first one, turns a GFF3 file into an in-memory database using gffutils. The
second method, iterares over all primary miRNAS in the database to extend
its mature forms' end and start cooridnates by n nucleotides without exceeding
the chromosome boundaries. Thereafter, the primary miRNAs star/end coordinates
will also be extened if and only if, the mature miRNA coordinates exceeds the
primary miRNA ones; this extension will be recorded as _-y and/or _+x being the
the extension on the 5' and 3' respectively. Finally, two different annotation
files will be created; one for the primary miRNAs and one for the mature forms.

Usage:
    mirna_extension.py GFF3 [--outdir Path] [-e int] [--chr Path]
    GFF3 > mirna_extension.py [--outdir Path] [-e int] [--chr Path]
"""

import argparse
from pathlib import Path
import sys
from typing import Dict, Optional
import gffutils  # type: ignore


class MirnaExtension:
    """Class to extend miRNAs start and end coordinates.

    Attributes:
        db:
            in-memory database of the GFF3 file.

    Methods:
        load_gff_file:
            creates an in-memory db from an input file
        extend_mirnas:
            extends the start and end of the mature miRNAs by n nucleotides
    """

    def __init__(self) -> None:
        """Initialize class."""
        self.db = None

    def load_gff_file(self, gff_file: Optional[Path] = None) -> None:
        """Load GFF3 file.

        This method uses the gffutils package to create and in-memory database
        of the GFF3 file.

        Args:
            gff_file:
                path to the input GFF3 file. If None, input is taken from the
                standard input (stdin).
        """
        if gff_file is None:
            self.db = gffutils.create_db(
                sys.stdin, dbfn=":memory:", force=True, keep_order=True
            )
        else:
            self.db = gffutils.create_db(
                str(gff_file), dbfn=":memory:", force=True, keep_order=True
            )

    def adjust_names(
        self,
        primir: gffutils.feature.Feature,
        mirna_list: list[gffutils.feature.Feature],
    ) -> None:
        """Adjust miRNAs' `Name` attribute.

        This method adjusts the miRNAs' `Name` attribute to account for the
        different genomic locations the miRNA sequence is anotated on making
        them unique.

        miRNA primary transcript names have either the format
        `SPECIES-mir-NAME` or `SPECIES-mir-NAME-#` where `#` is the integer
        indicating which replicate is it. If there are several entries with
        the exact same name, the feature `ID` has the format `ID_#`.

        Mature miRNA names present one of the following formats:
            - `SPECIES-miR-NAME-ARM`
            - `SPECIES-miR-NAME-#-ARM`
        Or, if there's a single mature miRNA:
            - `SPECIES-miR-NAME`
            - `SPECIES-miR-NAME-#`

        If one pri-miR is considered to have multiple entries, the suffix of
        either the `Name` or the `ID` is set in its corresponding mature
        miRNA(s) name. The pri-miR `Name` is also updated if the suffix is
        found in the `ID`.

        Args:
            primir:
                miRNA primary transcript feature entry
            mirna_list:
                list with the corresponding mature miRNA feature(s) entry(s)
        """
        suffix = None
        primir_name = primir.attributes["Name"][0].split("-")
        primir_id = primir.attributes["ID"][0].split("_")

        if len(primir_name) == 4:
            suffix = primir_name[-1]

        elif len(primir_name) == 3 and len(primir_id) == 2:
            suffix = primir_id[-1]
            primir_name.append(suffix)

        if suffix:
            for mirna in mirna_list:
                mirna_name = mirna.attributes["Name"][0].split("-")

                if len(mirna_name) == 5:
                    mirna_name[3] = suffix

                elif len(mirna_name) == 4 and mirna_name[-1] != suffix:
                    mirna_name.insert(3, suffix)

                elif len(mirna_name) == 3:
                    mirna_name.append(suffix)

                mirna.attributes["Name"][0] = "-".join(mirna_name)

        primir.attributes["Name"][0] = "-".join(primir_name)

    def extend_mirnas(
        self,
        primir_out: Path,
        mir_out: Path,
        n: int = 6,
        seq_lengths: Optional[dict[str, int]] = None,
    ) -> None:
        """Extend miRNAs start and end coordinates.

        This method elongates the start and end coordinates of mature miRNAs
        by n nucleotides. In the case that this extension makes the start/end
        coordinates to exceed the corresponding primary miRNA boundaries,
        these will be extended as far as the miRNA coordinates.
        If provided, the elongation will take into account the chromosome size.

        Args:
            primir_out:
                path to the output extended precursor miRNA file.
            mir_out:
                path to the output extended mature miRNA file.
            n:
                number of nucleotides to extend miRs start and end coordinates.
            seq_lengths:
                a dictionary that maps reference sequence IDs to their lengths
                (in  nucleotides). If None, sequence lengths are inferred from
                the input GFF3 file.
        """
        assert isinstance(self.db, gffutils.interface.FeatureDB)

        # Set end boundary
        if seq_lengths is None:
            seq_lengths = {}
            for seqid in self.db.seqids():
                seq_lengths[seqid] = max(
                    rec.end for rec in self.db.region(seqid)
                )

        with (
            open(primir_out, "w", encoding="utf-8") as primir,
            open(mir_out, "w", encoding="utf-8") as mirna,
        ):
            for primary_mirna in self.db.features_of_type(
                "miRNA_primary_transcript"
            ):
                seqid = primary_mirna.seqid
                min_start = primary_mirna.start
                max_end = primary_mirna.end

                mature_miRNAs = list(
                    self.db.region(
                        seqid=seqid,
                        start=min_start,
                        end=max_end,
                        strand=primary_mirna.strand,
                        featuretype="miRNA",
                        completely_within=True,
                    )
                )

                self.adjust_names(primary_mirna, mature_miRNAs)

                if mature_miRNAs:
                    for mature_mir in mature_miRNAs:
                        if mature_mir.start - n > 0:
                            mature_mir.start -= n
                        else:
                            mature_mir.start = 0

                        if mature_mir.end + n < seq_lengths[seqid]:
                            mature_mir.end += n
                        else:
                            mature_mir.end = seq_lengths[seqid]

                        min_start = min(min_start, mature_mir.start)
                        max_end = max(max_end, mature_mir.end)

                        mirna.write(str(mature_mir) + "\n")

                    start_diff = primary_mirna.start - min_start
                    end_diff = max_end - primary_mirna.end

                    primary_mirna.start = min_start
                    primary_mirna.end = max_end

                    primary_mirna.attributes["Name"][0] += f"_-{start_diff}"
                    primary_mirna.attributes["Name"][0] += f"_+{end_diff}"
            
                primir.write(str(primary_mirna) + "\n")


def parse_arguments():
    """Command-line arguments parser."""
    parser = argparse.ArgumentParser(
        description="Script to extend miRNAs start and end coordinates."
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s 1.0",
        help="Show program's version number and exit",
    )
    parser.add_argument(
        "input",
        help="Path to the GFF3 annotation file. If not provided, the input \
             will be read from the standard input.",
        type=Path,
    )
    parser.add_argument(
        "--outdir",
        help="Path to the output directory. Default: %(default)s.",
        default=Path.cwd(),
        type=Path,
    )
    parser.add_argument(
        "-e",
        "--extension",
        help="Number of nucleotides to extend the coordinates. Default: \
                %(default)d.",
        default=6,
        type=int,
    )
    parser.add_argument(
        "--chr",
        help="Path to the tabulated file containing the chromosome and its \
            length in basepairs. If not provided, the length will be set to \
            the biggest coordinate of the last miRNA primary transcript.",
        default=None,
        type=Path,
    )

    return parser


def main(arguments) -> None:
    """Extend miRNAs start/end coordinates."""
    outdir = Path(arguments.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    primir_out = outdir / (
        f"extended_primir_annotation_{arguments.extension}_nt.gff3"
    )
    mir_out = (
        outdir / f"extended_mirna_annotation_{arguments.extension}_nt.gff3"
    )

    with open(arguments.input, encoding="utf-8") as in_file:
        if len(in_file.read()) == 0:
            with (
                open(primir_out, "w", encoding="utf-8") as primir,
                open(mir_out, "w", encoding="utf-8") as mir,
            ):
                primir.write("")
                mir.write("")
                return

    # Create dictionary with the ref. sequence length
    seq_lengths: Optional[Dict] = None
    if arguments.chr:
        seq_lengths = {}
        with open(arguments.chr, encoding="utf-8") as f:
            for line in f:
                ref_seq, length = line.strip().split("\t")
                seq_lengths[ref_seq] = int(length)

    m = MirnaExtension()
    m.load_gff_file(arguments.input)
    m.extend_mirnas(
        n=arguments.extension,
        seq_lengths=seq_lengths,
        primir_out=primir_out,
        mir_out=mir_out,
    )


if __name__ == "__main__":
    args = parse_arguments().parse_args()  # pragma: no cover
    main(args)  # pragma: no cover
