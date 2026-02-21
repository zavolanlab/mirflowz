<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/sam_uncollapse.pl#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>script</kbd> `sam_uncollapse.pl`

Reverses the collapsing of reads with identical sequences as done with
`fastx_collapser` ([FASTX Toolkit][docs-fastx]) or similar.

Reads and writes files in SAM format. Each line is printed `n` times, where `n`
is the suffix appended to the read/query name via a dash (`-`).


!!! warning "CAUTION"

    Only marginal validation of the input file type/format performed!

Usage
-----

```bash
perl sam_uncollapse.pl [OPTIONS] --in [FILE|SAM] --out [FILE|SAM]
```

Arguments
---------

- <b>`--in [FILE|SAM]` (required)</b>: Path to the input SAM file.
- <b>`--out [FILE|SAM]` (required)</b>: Path to the output SAM file.

Options
-------

- <b>`--suffix`</b>: Add serial number suffix to each `QNAME` during
  uncollapsing (separated by a ".") to allow distinction of multimappers by
  `QNAME`.
- <b>`-h`</b> | <b>`--help`</b>: Show this information and die.
- <b>`-u`</b> | <b>`--usage`</b>: Show this information and die.
- <b>`--quiet`</b>: Shut up!

Requirements
------------

- <b>Perl version</b>: `>= 5.40.2`
- <b>Modules:</b>
    - <b>`Getopt::Long`</b>: `>= 2.58`

---


<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/sam_uncollapse.pl#L73"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>subroutine</kbd> `usage`


Returns usage information for current script

**Accepts**

N/A

**Returns**

String with usage information

**Type**

Specialized


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/sam_uncollapse.pl#L100"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>subroutine</kbd> `sam_uncollapse`

For each line of a SAM file, parses the identifier `QNAME` for the presence of
a number `n` appended to its end via a dash ('-') and re-writes the line `n`
times.

Header lines are reproduced as they are.


**Accepts**

1. Input file [FILE|SAM]
2. Output file [FILE|SAM]
3. Suffix switch: 0 = Do not add serial number suffix, 1 = Add serial number
   suffix

**Type**

Generic

---
