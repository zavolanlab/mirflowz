<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `mirna_extension.py`

Extend miRNA start and end coordinates and ensure name uniqueness.


This method updates the attributes of precursor and mature miRNA
entries to ensure consistent naming based on their paralog or sequence
variant status.

**For precursors:**

- A suffix indicates distinct genomic loci (paralogs) that express
  identical mature sequences. This is typically extracted from the 'Name'
  or 'ID' attribute.
- Format: `SPECIES-mir-NUMBER[LETTER]-#` (Name) and `ALIAS_#` (ID)
  where:
     - `LETTER` denotes a sequence variant of the mature miRNA (paralogous
       variant with similar but not identical sequences),
     - `#` indicates the paralog number (replica/locus index), included
        when multiple loci express the same or similar miRNAs.

**For mature miRNAs:**

- The replica number is added or replaced as an infix/suffix in the name.
- Formats:
    - `SPECIES-miR-NUMBER[LETTER]-#-ARM`
    - `SPECIES-miR-NUMBER[LETTER]-#`
    - `SPECIES-miR-NUMBER[LETTER]-ARM`

**Cases:**

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
with `_-n` and/or `_+m`, where n and m represent the extensions on the 5' and
3' end, respectively (or 0 otherwise). The modified mature miRNA and precursor
annotations will be written to separate GFF3 files.

**Constraints:**

The script was written for GFF3 files obtained from [miRBase][mirbase].

The following **assumptions** are made:

- The input GFF3 file contains miRNA annotations for a single species.
- The input GFF3 file contains features of type "miRNA_primary_transcript"
  (referred to as "precursors") and "miRNA" (referred to as "mature miRNAs").
- The input GFF3 file contains a 'Name' and 'ID' attribute for each precursor.
- Each precursor contains one or more mature miRNAs.
- Each mature miRNA is a child of exactly one precursor and is completely
  within the boundaries of the precursor.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L299"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_arguments`

```python
parse_arguments()
```

Parse command-line arguments.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L410"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main`

```python
main(args) → None
```

Extend miRNAs start/end coordinates.


---

## <kbd>class</kbd> `AnnotationException`
A custom exception class for `MirnaExtension` class.


---

## <kbd>class</kbd> `MirnaExtension`
Class for updating miRNA annotated coordinates and names.


**Attributes:**

 - <b>`db`</b>:  In-memory database of the input GFF3 records.
 - <b>`db_out`</b>:  In-memory database of the updated GFF3 records.
 - <b>`seq_lengths`</b>:  Dictionary mapping reference sequence IDs to their
   lengths.

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L25"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `__init__`

```python
__init__() → None
```

Initialize class.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L98"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `adjust_names`

```python
adjust_names(precursor: Feature, matures: list[Feature]) → None
```

Adjust miRNA attributes for uniqueness and consistency.

This method updates the attributes of precursor and mature miRNA entries to
ensure consistent naming based on their paralog or sequence variant status.

For precursors:

- A suffix indicates distinct genomic loci (paralogs) that express identical
  mature sequences. This is typically extracted from the `Name` or `ID`
  attribute.
- Format: `SPECIES-mir-NUMBER[LETTER]-#` (Name) and `ALIAS_#` (ID) where:
    - `LETTER` denotes a sequence variant of the mature miRNA (paralogous
      variant with similar but not identical sequences),
    - `#` indicates the paralog number (replica/locus index), included when
      multiple loci express the same or similar  miRNAs.

For mature miRNAs:

- The replica number is added or replaced as an infix/suffix in the name.
- Formats:
    - `SPECIES-miR-NUMBER[LETTER]-#-ARM`
    - `SPECIES-miR-NUMBER[LETTER]-#`
    - `SPECIES-miR-NUMBER[LETTER]-ARM`

Cases:

- If a precursor has multiple genomic instances (paralogs), the first occurrence
  typically lacks a numeric suffix; subsequent  ones are numbered incrementally.
- The `Derives_from` attribute of each mature miRNA is updated to match the
  precursor's `ID`.
- If a precursor has a different sequence variant designation (`LETTER`) than
  its associated matures, the mature miRNA names are updated to match the
  precursor's designation.


The `Alias` attribute remains unchanged.


**Arguments:**

 - <b>`precursor`</b>:  "miRNA primary transcript" feature entry
 - <b>`matures`</b>: list of the corresponding "mature miRNA" feature(s)
   entry(s)

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L187"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `process_precursor`

```python
process_precursor(precursor: Feature, n: int = 6) → list[Feature]
```

Extend miRNAs start and end coordinates and ensure name uniqueness.

This method elongates the start and end coordinates of mature miRNAs by `n`
nucleotides. In the case that this extension makes the start/end coordinates
to exceed the corresponding primary miRNA boundaries, these will be extended
as far as the miRNA coordinates. If provided, the elongation will take into
account the chromosome size.

In addition, the method `adjust_names` is called to ensure uniqueness in the
`Name` attribute for both the precursor and its arms.


**Arguments:**

 - <b>`precursor`</b>:  "miRNA primary transcript" feature entry
 - <b>`n`</b>:  Number of nucleotides to extend miRs start and end coordinates.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L31"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `set_db`

```python
set_db(path: Path) → None
```

Load GFF3 file into `gffutils.FeatureDB`.


**Arguments:**

 - <b>`path`</b>:  Path to a GFF3 file.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L52"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `set_seq_lengths`

```python
set_seq_lengths(path: Optional[Path] = None) → None
```

Set the reference sequence lengths.


**Arguments:**

 - <b>`path`</b>:  Path to a tabulated file containing the names of reference
   sequences and the corresponding total length, in base pairs.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L244"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `update_db`

```python
update_db(n: int = 6) → None
```

Update miRNA annotations in the local database.

Using the method `process_precursor` annotated coordinates are extended by
`n` nucleotides and, if needed, the `Name` attribute is modified to contain the
sequence replica number.


**Arguments:**

 - <b>`n`</b>:  Number of nucleotides to extend miRs start and end coordinates.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_extension.py#L274"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `write_gff`

```python
write_gff(path: Path, feature_type: Optional[str] = None) → None
```

Write features to a GFF3 file.


**Arguments:**

- <b>`path`</b>:  Path to the output file.
- <b>`feature_type`</b>:  Feature type to write. If `None`, all features will
   be written.

---
