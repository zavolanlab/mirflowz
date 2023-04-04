#! usr/bin/env python3
"""Extend primary miRNAs overhang.

This class is build upon two methods. The first one, turns a GFF3 file
into an in-memory database using gffutils. The second method, iterares over
all primary miRNAS in the database and determines the distances to the
nearest mature miRNAs. If either distance is less than n, the function
extends the start/end coordinates of the pirmary miRNA by the necessary
amount, making sure not to extend beyond the borders of the reference
sequence. Finally, it prints to the standard output (stdout) the updated
coordinates.
"""
import gffutils
import sys


class MirnaExtension():
    """Class to extend primary miRNAs overhang.

    Attributes:
        gff_file(str):
            Path to the input GFF3 file. If None, input is taken from the
            standard input (stdin).
        n(int):
            The minimum distance (in nucleotides) from the start/end of the
            primary miRNA to the start/end of the nearest mature miRNA.
        seq_lengths(dict[str, int]):
            A dictionary that maps reference sequences IDs to their lengths
            (in nucleotides). If None, sequence lengths are inferred from
            the input GFF3 file.
        db(FeatureDB):
            in-memory database of the GFF3 file.

    Methods:
        load_gff_file:
            creates an in-memory db from the input file
        extend_primary_mirnas:
            extends the start and/or end of the primary miRNAs if needed
    """

    def __init__(self, gff_file: str = "", n: int = 10, seq_lengths: dict[str, int] = dict()) -> None:
        """Initialize class.

        Args:
            gff_file(str) :
                path to the input GFF3 file. If None, input is taken from the
                standard input (stdin).
            n(int) :
                the minimum distance (in nucleotides) from the start/end of the
                primary miRNA to the start/end of the nearest mature miRNA.
            seq_lengths(dict[str, int]) :
                a dictionary that maps reference sequence IDs to their lengths
                (in  nucleotides). If None, sequence lengths are inferred from
                the input GFF3 file.
        """
        self.gff_file = gff_file
        self.n = n
        self.seq_lengths = seq_lengths
        self.db = None

    def load_gff_file(self) -> None:
        """Load GFF3 file.

        This method uses the gffutils package to create and in-memory database
        of the GFF3 file.
        """
        if self.gff_file == "":
            self.db = gffutils.create_db(sys.stdin, dbfn=':memory:', force=True, keep_order=True)
        else:
            self.db = gffutils.create_db(self.gff_file, dbfn=':memory:', force=True, keep_order=True)
    
    def extend_primary_mirnas(self) -> None:
        """Extend primary miRNAs start/end."""
        # Set end boundary
        if self.seq_lengths is None:
            self.seq_lengths = {}
            for seqid in self.db.seqsids:
                self.seq_lengths[seqid] = max(rec.end for rec in self.db.region(seqid))

        # Iterate over primary miRNAs
        for primary_mirna in self.db.features_of_type('miRNA_primary_transcript'):

            seqid = primary_mirna.seqid
            seq_len = self.seq_lengths[seqid]
            start = int(primary_mirna.start)
            end = int(primary_mirna.end)

            # Determine distances to the nearest mature miRNAs
            left_distance = start
            right_distance = seq_len - end

            mature_miRNAs = list(self.db.region(seqid=primary_mirna.seqid, start=primary_mirna.start, end=primary_mirna.end, featuretype='miRNA', completely_within=True))

            if mature_miRNAs:
                mature_start = min(m.start for m in mature_miRNAs)
                mature_end = max(m.end for m in mature_miRNAs)
                left_distance = mature_start - start
                right_distance = end - mature_end

            # Extend the primary miRNA if needed
            left_shift = self.n - left_distance
            right_shift = self.n - right_distance

            if left_distance < self.n:
                if start - left_shift < 0:
                    left_shift = start
                primary_mirna.start -= left_shift

            if right_distance < self.n:
                if end + right_shift > seq_len:
                    right_shift = seq_len - end
                primary_mirna.end += right_shift

            # Write to standard output
            sys.stdout.write(str(primary_mirna) + '\n')
            for child in mature_miRNAs:
                sys.stdout.write(str(child) + '\n')
