<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/ascii_pileups_aesthetics_modification.R#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>script</kbd> `ascii_pileups_aesthetics_modification.R`

Enhance ASCII-style alignment pileups.

Filter each ASCII-style alignment pileup in the provided input directory
retaining those alignments with at least  `--min-count` amount of reads that
appear in the first `--max-sequences` lines. If `--canonical` is set, mark the
read(s) matching the canonical sequence(s) inferred from the reference row in
the output.

Optionally, if `--split-arms` is set, split the precursor pileup into
arm-specific pileup(s) with an allowed overhang of &pm; `--overhang`
nucleotides on either side of the mature arm, and the genomic coordinates are
adjusted to the final representation. If no overhang is provided, only reads
fully contained within the exact arm span are kept.

If `--keep-all` is set, the ASCII-style alignment pileup is written even if it
has no aligned sequences.

**EXPECTED INPUT**

Tab-separated pileup files with two columns without header:

1) Alignment string or reference sequence

2) Feature name, genomic coordinates, or read count



The ASCII-style alignment representation is expected to be one of

1) One arm with precursor

2) Two arms with precursor


**EXPECTED OUTPUT**

For each ASCII-style alignment pileup, one TSV is written to `--out-dir` named
`PREFIX.FEATURE_NAME.OVERHANG-shift.tab` where:

- `PREFIX` is provided by the option `--prefix`
- `FEATURE_NAME` is either the precursor being represented, or the mature arm
  if `--split-arms` is set
- `OVERHANG` is either 0, or provided by the option `--overhang` when
  `--split-arms` is set


Usage
-----

```bash
Rscript ascii_pileups_aesthetics_modification.R [--help] [--verbose] [OPTIONS] --in-dir=[DIR] --prefix=[PREFIX]
```


Arguments
---------

- <b>`--in-dir=DIR` (required)</b>: Absolute path from where input files shall
  be read.
- <b>`--prefix=STR` (required)</b>: Prefix to be used in the output file
  name(s).


Options
-------

- <b>`--out-dir DIR`</b>: Absolute path to where output files shall be written.
- <b>`--split-arms`</b>: Split precursor pileups into one mature-arm pileup per
  arm. All subsequent filtering and formatting steps are then applied
  independently to each arm.
- <b>`--min-count INT`</b>: Minimum count for a sequence to be kept
  (default: 1).
- <b>`--max-sequences INT`</b>: Maximum number of top sequences to be
  displayed. It is assumed that the input ASCII-style alignment pileups are
  already sorted in the desired order (default: 30).
- <b>`--overhang INT`</b>: If `--split-arms` is set, number of extra positions
  to retain on each side of the mature arm span. Reads extending beyond this
  interval are removed. If omitted, only reads fully contained within the
  exact arm span are kept.
- <b>`--canonical`</b>: Mark the aligned read(s) corresponding to the canonical
  mature sequence(s).
- <b>`--keep-all`</b>: Write the ASCII-style alignment pileups even if it has
  no aligned sequences.
- <b>`-h`</b> | <b>`--help`</b>: Show this information and die.
- <b>`-u`</b> | <b>`--usage`</b>: Show this information and die.
- <b>`-v`</b> | <b>`--verbose`</b>: Print log messages to `STDOUT`.


Dependencies
------------

- <b>R version</b>: `>= 3.6.0`
- <b>R packages</b>:
    - <b>`optparse`</b>: `>= 1.6.2`
    - <b>`dplyr`</b>: `>= 1.1.4`

---


<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/ascii_pileups_aesthetics_modification.R#L241"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `GetSpan`

```r
GetSpan <- function(str)
```

Get the span of the non-padding portion of a pileup string.

`GetSpan()` returns the first and last positions in a string that are not
`"."`. In the context of ASCII-style pileups, this identifies the coordinates
occupied by the represented element, such as a mature-arm marker, an aligned
read, or a reference sequence.

The function is intended for pileup strings in which `"."` marks padding or
empty positions. It can therefore be used on:

- mature-feature rows represented with `">"` characters
- aligned read strings
- genomic/reference sequence rows

In all cases, the function returns the interval covered by the non-dot
characters. If the string contains only `"."` characters, or no characters at
all, the function returns `NA` for both boundaries.

**Arguments:**

- <b>`str`</b>: A character string from an ASCII-style pileup.


**Returns:**

A named integer vector of length two with elements:

- <b>`start`</b>. The first position in `str` that is not `"."`.
- <b>`end`</b>: The last position in `str` that is not `"."`.

If no such position exists, both values are `NA_integer_`.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/ascii_pileups_aesthetics_modification.R#L278"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `WithinSpan`

```r
WithinSpan <- function(str, ref.span, overhang)
```

Test whether an aligned substring lies within a reference span.

`WithinSpan()` checks whether the aligned portion of a string falls fully
inside a reference span, optionally extended by a user-defined overhang.

**Arguments:**

- <b>`str`</b>: A character string representing an aligned sequence padded with
  `"."`.
- <b>`ref.span`</b>: A named vector with elements `start` and `end` giving the
  reference span.
- <b>`overhang`</b>: A non-negative integer specifying how many positions
  outside the reference span are tolerated on each side.


**Returns:**

`TRUE` if the aligned substring in `str` lies fully within the interval
`[ref.span["start"] - overhang, ref.span["end"] + overhang]`; otherwise `FALSE`.
If `str` contains no aligned substring, the function returns `FALSE`.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/ascii_pileups_aesthetics_modification.R#L307"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `CutToSpan`

```r
CutToSpan <- function(str, span, overhang)
```

Trim a string to a reference span plus optional overhang.

`CutToSpan()` extracts the substring defined by a reference span extended by
an optional overhang. The extracted interval is clipped to the actual string
boundaries.


**Arguments:***

- <b>`str`</b>: A character string to trim.
- <b>`span`</b>: A named vector with elements `start` and `end`.
- <b>`overhang`</b>: A non-negative integer specifying how many extra positions
  to retain on each side of the span.


**Returns:**

A character string corresponding to the requested interval.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/ascii_pileups_aesthetics_modification.R#L346"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `AdjustRefCoords`

```r
AdjustRefCoords <- function(coord, span, overhang, ref.seq)
```

Adjust genomic coordinates after trimming the reference sequence.

`AdjustRefCoords()` recalculates the genomic coordinates of a reference row
after trimming the displayed reference sequence to an arm span plus optional
overhang. The function assumes that character positions in `ref.seq` correspond
directly to increasing genomic coordinates in `coord`.


**Arguments:**

- <b>`coord`</b>: A character string with genomic coordinates in the format
  `"chr:start-end:strand"`.
- <b>`span`</b>: A named vector with elements `start` and `end` defining the
  arm span in the reference sequence.
- <b>`overhang`</b>: A non-negative integer specifying how many extra positions
  to keep on each side of the arm span.
- <b>`ref.seq`</b>: The full reference sequence string before trimming.


**Returns:**

A character string containing the adjusted genomic coordinates in
  the same format as the input.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/ascii_pileups_aesthetics_modification.R#L405"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `SplitArms`

```r
SplitArms <- function(in.pileup, head.lines, overhang)
```

Split a precursor pileup into arm-specific pileups.

`SplitArms()` splits a precursor ASCII-style alignment pileup into one pileup
per mature arm. Each output pileup contains the arm annotation row, the
trimmed reference row, and only those reads whose aligned span falls within
the arm span plus the specified overhang.

For each mature arm:

1. Identify the arm span from the arm annotation row
2. Trim the arm row and reference row to that span plus overhang
3. Adjust the genomic coordinates of the reference row accordingly
4. Retain only reads that lie fully within that interval
5. Trim retained reads to the same interval


**Arguments:**

- <b>`in.pileup`</b>: A data frame representing one input pileup. It must
  contain the columns `seq` and `counts`.
- <b>`overhang`</b>: A non-negative integer specifying the allowed extension
  beyond the arm span on each side. If `NULL`, it is treated as `0`.
- <b>`head.lines`</b>: An integer giving the number of header lines in the
  input pileup. Expected values are: `3` if single arm with precursor, or
  `4` if two arms with precursor.


**Returns:**

A list of data frames, one per mature arm.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/ascii_pileups_aesthetics_modification.R#L483"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `FilterPileup`

```r
FilterPileup <- function(pileup, min.count, max.seq)
```

Filter a pileup by read count and maximum number of displayed sequences.

`FilterPileup()` keeps only read rows with counts greater than or equal to
`min.count`, then keeps at most the first `max.seq` of those rows. The header
rows are preserved, and a separator row labeled `"Counts"` is inserted before
the retained reads. The function assumes that rows above the genomic reference
row are header rows and that rows below it correspond to aligned reads.

**Arguments:**

- <b>`pileup`: A data frame representing one pileup, with columns `seq` and
  `counts`.
- <b>`min.count`</b>: Minimum count required for a read to be retained.
- <b>`max.seq`</b>: Maximum number of read rows to retain.


**Returns:**

A filtered pileup data frame.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/ascii_pileups_aesthetics_modification.R#L535"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `AddCanonical`

```r
AddCanonical <- function(pileup)
```

Mark the canonical aligned read in a pileup.

`AddCanonical()` adds a column named `feat` to a pileup and marks the row or
rows corresponding to the canonical mature sequence for each detected arm.

The canonical sequence is inferred directly from the reference row using the
mature arm span:

1. identify the genomic reference row
2. identify mature arm rows located above the reference row
3. extract the arm span from each mature arm row
4. extract the corresponding substring from the reference row
5. match this sequence against read rows after removing `"."` padding
6. exclude any read containing `"-"` from being considered canonical

The function assumes that the row immediately below the reference row is the
`"Counts"` separator row and therefore does not consider it a read.


**Arguments:**

- <b>`pileup`</b>: A data frame representing one pileup, with columns `seq` and
  `counts`.

**Returns:**

The input pileup with an additional column `feat`. Rows matching the canonical
sequence are marked with `"< <feature_name>"` in that column, while all other
rows contain an empty string.

---
