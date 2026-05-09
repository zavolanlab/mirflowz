<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/ascii_alignment_pileup.R#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>script</kbd> `ascii_alignment_pileup.R`

Generates an ASCII-style pileup of read alignments in one or more BAM files
against one or more regions specified in a BED file.

For each BED interval, the script extracts overlapping reads, represents each
read as a fixed-width ASCII string (padded to the region width), optionally
collapses identical alignments into counts, and writes a tab-delimited output
file.

Optionally, a reference sequence (FASTA; `BGZIP`-compressed) can be included as
an additional row, and overlapping annotation features from a GFF/GTF file can
be rendered as arrow-like tracks.

**STRAND HANDLING**

If the BED interval has strand `-`, read strings are complemented to reflect
the opposite strand and are also reversed for sorting and display. If the BED
strand is `+` (or unset/`*`), sequences are displayed as-is.

**ASCII CONVENTIONS**

- Padding character: controlled by `--padding-character` (default `"."`).
- Indel character: controlled by `--indel-character` (default `"-"`), passed to
  the underlying pileup/stacking routine.

**COLLAPSING / COUNTS**

By default, identical alignment strings are collapsed and reported with a count
(frequency). Use `--do-not-collapse-alignments` to output each alignment row
individually (count field empty). When collapsing, alignments with count less
than `--minimum-count` are filtered out.

!!! warning "Alignment insertions"

    Whenever an alignment contains an insertion, it is removed from the read
    sequence. The resulting sequence is counted towards the read with an
    identical sequence, or displayed on its own if the option
    `--do-not-collapse-alignments` is used.

**SORTING**

ASCII-style alignment pileups' final arrangement is controlled by the option
`--sort-by` and the flag `--reverse-sort`.

By default, pileups are sorted by the alignment counts in descending order. If
two alignments have the same amount of counts, these are sorted
from-left-to-right by their first and last nucleotide position.

Use `--sort-by="position"` to arrange the alignments by their first
nucleotide's position (from-left-to-right). If more than one alignment starts
at the same location, these are additionally sorted by counts
(descending order), and by their last nucleotide's position
(from-left-to-right).

Both sorting types can be reversed by using the flag `--reverse-sort`:
from "from-left-to-right" to "from-right-to-left" for sort type "position",
and from "descending" to "ascending" for sort type "counts".

**EXPECTED OUTPUT**

For each BED interval, one TSV is written to `--output-directory`:

- Collapsed (default): `PREFIX.REGION_NAME.min.MINIMUM_COUNT.pileup.tab`
- Uncollapsed: `PREFIX.REGION_NAME.uncollapsed.pileup.tab`

where:

- `PREFIX` is either provided by the option `--prefix`, or the input BAM
  file(s) name otherwise
- `REGION_NAME` is the region name annotated in the BED file the pileup has
  been done for
- `MIN_COUNT`is the minimum amount of collapsed alignments a sequence has to
  have in order to be included. Provided by the option `--minimum-count` with
  a default value of 1


Output has two columns (no header):

1) `seq`: ASCII alignment string (region-width)

2) `count`: integer count (collapsed mode) or label for reference/annotation
   rows or empty/NA (uncollapsed mode)

In addition:

- If `--reference` is provided, a reference row is prepended with `count` set to
  `<seqname>:<start>-<end>:<strand>`
- If `--annotations` is provided, annotation rows are prepended with `count` set
  to the feature label taken from `--annotation-name-field` (default: `"Name"`).


Usage
------

```bash
Rscript ascii_alignment_pileup.R [--help] [--verbose] [OPTIONS] BED BAM [BAM2 ...]
```

Arguments
---------

- <b>`BAM [BAM2 ...]`</b>: BAM file(s) containing the read alignments.
- <b>`BED`</b>:  BED file with the region(s) to create the ASCII-style pileups
  for.


Options
-------

??? information "A note on efficiency"

    For the input queries, consider the `--maximum-region-width` parameter,
    which is provided for safety. While it is possible to increase it, wide
    regions may require excessive memory.

- <b>`--reference=FILE`</b>: Reference genome sequence in FASTA
  format **compressed with `BGZIP`**. If supplied, the reference sequence for
  the query region(s) will be added to the output. Note that on the first run
  with a specific reference genome file, an `FAI` index is generated which will
  take some time.
- <b>`--annotations=FILE`</b>: Annotation file in GFF3/GTF format used to
  annotate sequences. If supplied, features overlapping the query region(s) will
  be visualized in the output. Ensure that the argument to option
  `--annotation-name-field` corresponds to a field in the annotations,
  otherwise the script will fail.
- <b>`--output-directory=DIR`</b>: Output directory. One output file will be
  created for each region in `--bed` and the filenames will be generated from
  the basenames of the supplied BAM file(s) and the name field (4th column) of
  the BED file. (default: working directory)
- <b>`--maximum-region-width=INT`</b>: Maximum input region width. Use with
  care as wide regions will use excessive resources. (default: 200)
- <b>`--do-not-collapse-alignments`</b>: Show alignments of reads with identical
  sequences individually.
- <b>`--minimum-count=INT`</b>: Alignments of reads with less copies than the
  specified number will not be printed. Option is not considered if
  `--do-not-collapse-alignments` is set. (default: 1)
- <b>`--annotation.annotation-name-field=STR`</b>: Annotation field used to
  populate the `name` column in the output. (default: "Name")
- <b>`--padding-character=CHAR`</b>: Character used for padding alignments.
  (default: ".")
- <b>`--indel-character=CHAR`</b>: Character to denote insertions and deletions
  in alignments. (default: "-")
- <b>`--prefix=STRING`</b>: Prefix to be used in the output file name(s). If
  not provided the input BAM file(s) name will be used instead.
- <b>`--sort-by=STRING`</b>: Specify the sort type (either "position" or
  "counts"). (default = "counts")
- <b>`--reverse-sort`</b>: Reverse the sort order; from "from-left-to-right" to
  "from-right-to-left" for sort type "position", and from "descending" to
  "ascending" for sort type "counts".
- <b>`-h`</b> | <b>`--help`</b>: Show this information and die.
- <b>`-v`</b> | <b>`--verbose`</b>: Print log messages to `STDOUT`.

Dependencies
------------

- <b>R version</b>: `>= 3.6.0`
- <b>R packages</b>:
    - <b>`optparse`</b>: `>= 1.6.2`
    - <b>`rtracklayer`</b>: `>= 1.44.0`
    - <b>`GenomicAlignments`</b>: `>=3.22`
    - <b>`tools`</b>: `>=4.5.2`

---
