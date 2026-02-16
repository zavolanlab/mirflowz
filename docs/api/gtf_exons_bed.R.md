<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/gtf_exons_bed.1.1.2.R#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>script</kbd> `gtf_exons_bed.1.1.2.R`

Converts the exon entries of a GTF file to a BED file with one line per exon.

Reads a GTF file, extracts features of type `exon`, and writes them to a BED
file with one record per exon. The BED _name_ (4th column) is taken from a
user-selected GTF attribute (metadata column), and the BED _score_ column is
set to a constant value.

Imports the input GTF as a `GRanges`, filters rows where the feature type
equals `"exon"`, and checks that the requested attribute (from `--name`) exists
in the GTF metadata. Then, assigns:

- `name` = the chosen attribute value (_e.g._, `transcript_id`)
- `score` = the numeric constant provided via `--score`

Finally, exports to BED using `rtracklayer::export(format="bed")`.

**COORDINATE CONVENTIONS**

GTF is conventionally 1-based, while BED uses 0-based start coordinates. When
exporting a `GRanges` to BED, `rtracklayer::export()` handles the BED coordinate
requirements for the output file.

Usage
-----

```bash
Rscript gtf_exons_bed.1.1.2.R [--help] [--verbose] [OPTIONS] --gtf <PATH> --bed <PATH>
```

Arguments
---------

- <b>`-i FILE`</b> | <b>`--gtf=FILE` (required)</b>: Path to the input GTF file.
- <b>`-o FILE`</b> | <b>`--bed=FILE` (required)</b>: Path to the output BED
  file.

Options
---------

- <b>`-n STRING`</b> | <b>`--name=STRING`</b>: Attribute to be used for the
  'name' (fourth) column in the BED output file (check GTF file for available
  options; default: 'transcript_ID').
- <b>`-s NUM`</b> | <b>`--score=NUM`</b>: Score that should be set in each row
  of the output BED file (default: 0).
- <b>`-h`</b> | <b>`--help`</b>: Show this information and die.
- <b>`-u`</b> | <b>`--usage`</b>: Show this information and die.
- <b>`-v`</b> | <b>`--verbose`</b>: Print log messages to `STDOUT`.


Dependencies
------------

- <b>R version</b>: `>= 3.6.0`
- <b>R packages</b>:
    - <b>`optparse`</b>: `>= 1.6.2`
    - <b>`rtracklayer`</b>: `>= 1.44.0`

----
