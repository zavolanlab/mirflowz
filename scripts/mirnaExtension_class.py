#! usr/bin/env python3
"""Extend primary miRNAs overhang.

This class is build upon two methods. The first one, turns a GFF3 file
into an in-memory database using gffutils. The second method, iterares over
all primary miRNAS in the database to extend its mature forms' end and start
cooridnates by n nucleotides without exceeding the chromosome boundaries.
Thereafter, the primary miRNAs star/end coordinates will also be extened if
and only if, the mature miRNA coordinates exceeds the primary miRNA ones.
Finally, two different annotation files will be created; on for the primary
miRNAs and one for the mature miRNAs.
"""
import gffutils
import sys


class MirnaExtension():
    """Class to extend miRNAs start and end coordinates.

    Attributes:
        gff_file:
            Path to the input GFF3 file. If None, input is taken from the
            standard input (stdin).
        premir_out:
            Path to the output GFF3 file for pre-miRs
        mir_out:
            Path to the output GFF3 file for miRs
        n:
            number of nucleotides to extend miRs start and end coordinates.
        seq_lengths:
            A dictionary that maps reference sequences IDs to their lengths
            (in nucleotides). If None, sequence lengths are inferred from
            the input GFF3 file.
        db:
            in-memory database of the GFF3 file.

    Methods:
        load_gff_file:
            creates an in-memory db from the input file
        extend_primary_mirnas:
            extends the start and/or end of the primary miRNAs if needed
    """

    def __init__(self, premir_out: str, mir_out: str, gff_file: str = "", 
                 n: int = 6, seq_lengths: dict[str, int] = None) -> None:
        """Initialize class.

        Args:
            gff_file:
                path to the input GFF3 file. If None, input is taken from the
                standard input (stdin).
            premir_out:
                path to the output GF3 file for pre-miRs.
            mir_out:
                path to the output GFF3 file for miRs.
            n:
                number of nucleotides to extend miRs start and end coordinates.
            seq_lengths:
                a dictionary that maps reference sequence IDs to their lengths
                (in  nucleotides). If None, sequence lengths are inferred from
                the input GFF3 file.
        """
        self.gff_file = gff_file
        self.premir_out = premir_out
        self.mir_out = mir_out
        self.n = n
        self.seq_lengths = seq_lengths
        self.db = None

    def load_gff_file(self) -> None:
        """Load GFF3 file.

        This method uses the gffutils package to create and in-memory database
        of the GFF3 file.
        """
        if self.gff_file == "":
            self.db = gffutils.create_db(sys.stdin, dbfn=':memory:', 
                                         force=True, keep_order=True)
        else:
            self.db = gffutils.create_db(self.gff_file, dbfn=':memory:', 
                                         force=True, keep_order=True)
    
    def extend_mirnas(self) -> None:
        """Extend miRNAs start and end coordinates.
        
        This method elongates the start and end coordinates of mature miRNAs
        by n nucleotides. In the case that this extension makes the start/end
        coordinates to exceed the corresponding primary miRNA boundaries,
        these will be extended as far as the miRNA coordinates.
        If provided, the elongation will take into account the chromosome size.
        """
        # Set end boundary
        if self.seq_lengths is None:
            self.seq_lengths = {}
            for seqid in self.db.seqids():
                self.seq_lengths[seqid] = max(rec.end for rec in self.db.region(seqid))

        with open(self.premir_out, 'w') as premir, open(self.mir_out, 'w') as mirna:

            for primary_mirna in self.db.features_of_type('miRNA_primary_transcript'):
                seqid = primary_mirna.seqid
                seq_len = self.seq_lengths[seqid]
                start = int(primary_mirna.start)
                end = int(primary_mirna.end)

                mature_miRNAs = list(self.db.region(seqid=seqid,
                                                    start=start,
                                                    end=end,
                                                    featuretype='miRNA',
                                                    completely_within=True))

                if mature_miRNAs:

                    for mir in mature_miRNAs:
                        if mir.start - self.n > 0:
                            mir.start -= self.n
                        else:
                            mir.start = 0

                        if mir.end + self.n < seq_len:
                            mir.end += self.n
                        else:
                            mir.end = seq_len

                        if mir.start < start:
                            primary_mirna.start = mir.start
                        if mir.end > end:
                            primary_mirna.end = mir.end

                        mirna.write(str(mir) + '\n')

                premir.write(str(primary_mirna) + '\n')