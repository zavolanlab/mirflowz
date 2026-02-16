<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validate_bedtools_intersect.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `validate_bedtools_intersect.py`

Validation utilities for `bedtools intersect -wo -s` output.

The input must be the result of:

- `-a GFF3/GTF` (annotated features, 1-based)
- `-b BAM` (alignments, 0-based)
- `-wo` (report only overlaps and include overlap length)
- `-s` (strand-aware intersections)

Expected column order (16 columns):

- **1-9 Feature (GFF3/GTF)**: chr, source, type, start, end, score, strand,
  phase/frame, attributes
- **10-15 Read (BAM-derived)**: chr, start, end, name, score, strand
- **16** Overlap length (bp)

Exposes:

- `FileFormatError`: exception raised for malformed lines/fields.
- `Record`: a validated, parsed view of a single output line.
- `validate_first_n`: fail-fast validation of the first N lines of a file.
- `parse_all`: a streaming generator of `(line_number, Record)` tuples.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validate_bedtools_intersect.py#L441"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `validate_first_n`

```python
validate_first_n(intersect_file: Path, n: int = 10, sep: str = '\t') → None
```

Validate the first `n` lines of a `bedtools intersect -wo -s` out file.

Fails fast if any of the first `n` lines is invalid.


**Arguments:**

- <b>`intersect_file`</b>:  Path to the intersect output file.
- <b>`n`</b>:  Number of leading lines to validate.
- <b>`sep`</b>:  Field separator (default: `TAB`).



**Raises:**

- <b>`FileFormatError`</b>:  If any of the first `n` lines has an invalid
  format.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validate_bedtools_intersect.py#L472"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_all`

```python
parse_all(intersect_file: Path, sep: str = '\t') → Iterator[tuple[int, Record]]
```

Stream a file, yielding `(line_number, Record)` for each data line.


**Arguments:**

- <b>`intersect_file`</b>:  Path to the intersect output file.
- <b>`sep`</b>:  Field separator (default: `TAB`).


**Yields:**
 Tuples of `(line_number, Record)`.


**Raises:**

- <b>`FileFormatError`</b>:  If any line has an invalid format (iteration
  stops at  the first error).

---

## <kbd>class</kbd> `FileFormatError`

Raised when a line or field does not match the expected format.

---

## <kbd>class</kbd> `Record`

A single validated `bedtools intersect -wo -s` record.

The expected column order is listed in the module docstring.


**Attributes:**

- <b>`feat_chr`</b>:  Feature chromosome of scaffold (with or without the
  `chr` prefix).
- <b>`feat_source`</b>:  Program or data source that produced the feature
  data source.
- <b>`feat_type`</b>:  Feature type name.
- <b>`feat_start`</b>:  Feature start (1-based).
- <b>`feat_end`</b>:  Feature end (1-based).
- <b>`feat_score`</b>:  Floating-point value or `.` if missing.
- <b>`feat_strand`</b>:  Feature's strand defined as `+` (forward) or `-`
  (reverse).
- <b>`feat_phase_frame`</b>:  `0`, `1`, or `2` indicating the feature's first
  base position within a codon, or `.` if missing.
- <b>`feat_attrs`</b>: Dictionary of parsed GFF3/GTF-style attributes (keys
  lower-cased).


- <b>`read_chr`</b>:  Read chromosome of scaffold (with or without the `chr`
  prefix).
- <b>`read_start`</b>:  Read start (0-based).
- <b>`read_end`</b>:  Read end (0-based).
- <b>`read_name`</b>:  Read name.
- <b>`read_score`</b>:  Floating-point value or `.` if missing.
- <b>`read_strand`</b>:  Read's strand defined as `+` (forward) or `-`
  (reverse).


- <b>`overlap_len`</b>:  Number of overlapping base pairs between the feature
  and the read.

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validate_bedtools_intersect.py#L89"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `__init__`

```python
__init__(
    feat_chr: str,
    feat_source: str,
    feat_type: str,
    feat_start: int,
    feat_end: int,
    feat_score: Union[float, Literal['']],
    feat_strand: Literal['+', '-'],
    feat_phase_frame: Union[int, Literal['']],
    feat_attrs: Dict[str, str],
    read_chr: str,
    read_start: int,
    read_end: int,
    read_name: str,
    read_score: Union[float, Literal['']],
    read_strand: Literal['+', '-'],
    overlap_len: int
) → None
```

Initialize a validated record.

All arguments are expected to be already coerced to their target types.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validate_bedtools_intersect.py#L224"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `cross_checks`

```python
cross_checks() → None
```

Validate relationships across fields.

Ensures:

- Feature and read are on the same chromosome and strand.
- start <= end for both feature and read.
- Overlap > 0 (required by `-wo`) and equals the computed overlap.


**Raises:**

 - <b>`FileFormatError`</b>:  If any rule is violated.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validate_bedtools_intersect.py#L131"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>classmethod</kbd> `from_line`

```python
from_line(line: str, sep: str = '\t') → Record
```

Parse and validate a single `bedtools intersect -wo -s` line.


**Arguments:**

- <b>`line`</b>:  Raw line from the output file.
- <b>`sep`</b>:  Field separator (default: `TAB`).


**Returns:**

A validated `Record` instance.


**Raises:**

- <b>`FileFormatError`</b>:  If the column count is wrong or any field is
  invalid.

---
