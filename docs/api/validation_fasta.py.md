<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validation_fasta.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `validation_fasta.py`

Filter FASTA files.

Process both uncompressed and `gzip`-compressed FASTA files by trimming,
filtering, and validating sequence records based on user-defined criteria.

Sequence IDs are trimmed at the first occurrence of any specified characters
in `--trim` to standardize naming conventions. If not character string is
provided, the first white space is used.

To filter the FASTA file by sequence IDs, a text file, with one (trimmed) ID
per line, has to be passed to `--filter`. Whether to keep (`--mode k`) or
discard (`--mode d`) the sequences with those IDs must be specified.

Sequences exceeding a given length threshold (`--remove`) are excluded.

If a path is provided to `--idlist`, the resulting sequence IDs are written
one per line in a separate text file.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validation_fasta.py#L21"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_and_validate_arguments`

```python
parse_and_validate_arguments()
```

Parse and validate command-line arguments.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validation_fasta.py#L124"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `open_fasta`

```python
open_fasta(in_file: Path) → TextIO
```

Open a FASTA or FASTA.GZ for text‐mode reading.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validation_fasta.py#L146"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `write_id_file`

```python
write_id_file(out_file: Path, id_list: List[str]) → None
```

Write the final sequence IDs, one per line.


**Arguments:**

- <b>`out_file`</b>:  Path to the file where to write the IDs.
- <b>`id_list`</b>:  FASTA IDs to be written.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validation_fasta.py#L158"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `compile_trim_pattern`

```python
compile_trim_pattern(trim_str: str) → Pattern[str]
```

Get a compiled regex pattern to trim at a character's first occurrence.


**Arguments:**

- <b>`trim_str`</b>:  Characters used to determine where trimming occurs. If
  empty,  white space is used as the default delimiter.


**Returns:**
A compiled regex pattern that captures (1) everything up to the first match
and (2) the rest of the string.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validation_fasta.py#L174"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `trim_id`

```python
trim_id(seq_rec: SeqRecord, _pattern: Pattern[str]) → SeqRecord
```

Trim a FASTA ID using the first-occurrence of any character in _pattern.

All parameters must be passed by keyword.


**Arguments:**

- <b>`seq_rec`</b>:  A `Bio.SeqRecord.SeqRecord` to be trimmed in place.
- <b>`_pattern`</b>: (internal) a pre-compiled regex from `get_trim_pattern`.


**Returns:**
The same `SeqRecord`, with `.id` and `.description` possibly updated.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/validation_fasta.py#L198"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main`

```python
main(arguments) → None
```

Filter and process a FASTA file.

---
