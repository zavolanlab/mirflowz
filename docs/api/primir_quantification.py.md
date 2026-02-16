<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/primir_quantification.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `primir_quantification.py`

Tabulate `bedtools intersect -wo -s` output file.

For each intersecting feature from an INTERSECT file, calculate the sum of its
contributions and print the result in tab-delimited format.

**Contribution logic** (controlled by `--collapsed` and `--nh`):

- No flags: contribution = 1
- Only `--nh`: contribution = 1 / NH
- Only `--collapsed`: contribution = #reads / 1
- Both `--collapsed` and `--nh`: contribution = #reads / NH

The number of collapsed reads (`#reads`) and the NH value are inferred from the
read name, which must follow one of these formats:

- plain: READ
- NH only: READ_NH
- collapsed only: READ-#reads
- collapsed + NH: READ-#reads_NH

**EXPECTED INPUT FILE**

The expected INTERSECT file must be the output of:
`bedtools intersect -wo -s -a GFF3/GTF -b BAM`.
The first 10 lines are validated for format. If the INTERSECT file is empty, no
output is produced.

**OUTPUT TABLE FORMAT**

The basic output table has no header and two columns:

1) Feature identifier (taken from the attribute named by `--id`,
   which must match a key in the attributes column of the original
   GFF3/GTF features used in the bedtools call)

2) Feature contribution (float)

**Optional columns** (controlled by flags):

- `--read-ids` appends a semicolon-separated list of read IDs that overlap the
  feature (always added as the last column).
- `--feat-extension` appends two columns with 5′ and 3′ extension sizes. These
  are inferred by splitting the FEATURE-ID value (selected via `--id`)  on
  underscores, assuming the convention: `NAME_5EXT_3EXT`. If that pattern  isn't
  present in the selected attribute, the extension columns will be 0.



Examples
--------

### Example 1 | Contribution when using `--collapsed`

**Use case:**

A single feature with several intersecting reads. The flag `--collapsed` is
used, so contribution equals the # of reads per alignment.

**IN INTERSECT records:**

    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0     19      44413   44434   8-2     255     +       21
    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0     19      44413   44434   24-1    255     +       21

**Alignments:**

    Read ID: 8-2  Number of collapsed reads: 2  Contribution: 2

    Read ID: 24-1  Number of collapsed reads: 1  Contribution: 1

**OUT table:**

    hsa-mir-524_-0_+0      3



### Example 2 | Contribution when using `--nh`

**Use case:**

A single feature with several intersecting reads. The flag `--nh` is used, so
contribution equals 1/NH.

**IN INTERSECT records:**

    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0     19      44413   44434   8_1     255     +       21
    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0     19      44413   44434   24_1    255     +       21

**Alignments:**

    Read ID: 8_1  Number of mapped genomic loci: 1  Contribution: 1
    Read ID: 24_1  Number of mapped genomic loci: 1  Contribution: 1

**OUT table:**

    hsa-mir-524_-0_+0      2



### Example 3 | Contribution when using `--collapsed` and `--nh`

**Use case:**

A single feature with several intersecting reads. The flags `--nh` and
`--contribution` are used, so contribution equals # of reads/NH.

**IN INTERSECT records:**

    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0     19      44413   44434   8-2_1   255     +       21
    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0     19      44413   44434   24-1_1  255     +       21

**Alignments:**

    Read ID: 8-2_1  Number of collapsed reads: 2  Number of mapped genomic loci: 1  Contribution: 2/1 = 2
    Read ID: 23-1_1  Number of collapsed reads: 1  Number of mapped genomic loci: 1  Contribution: 1/1 = 1

**OUT table:**

    hsa-mir-524_-0_+0      3

### Example 4 | Column with intersecting reads; using `--read-ids`

**Use case:**

A single feature with several intersecting reads. Each read contributes with 1.
Read IDs intersecting the feature are appended as the last column.

**IN INTERSECT records:**

    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0     19      44413   44434   8-2_1   255     +       21
    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-0_+0     19      44413   44434   24-1_1  255     +       21

**Alignments:**

    Read ID: 8-2_1  Contribution: 1
    Read ID: 24-1_1  Contribution: 1

**OUT table:**

    hsa-mir-524_-0_+0      2       8-2_1;24-1_1

### Example 5 | Columns with feature shifts; using `--feat-extension`

**Use case:**

A single feature with several intersecting reads. Each read contributes with 1.
Feature start and end coordinates shift are appended as the third and fourth
column respectively.

**IN INTERSECT records:**

    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-1_+3     19      44413   44434   8-2_1   255     +       21
    19      .       miRNA_primary_transcript        44362   44448   .       +       .       ID=MI0003160;Alias=MI0003160;Name=hsa-mir-524_-1_+3     19      44413   44434   24-1_1  255     +       21

**Alignments:**

    Read ID: 8-2_1  Contribution: 1
    Read ID: 24-1_1  Contribution: 1

**Feature:**

    Feature name: hsa-mir-524_-1_+3
    5' shift: -1
    3' shift: +3

**OUT table:**

    hsa-mir-524_-1_+3      2       -1       +3

**Global Variables**
---------------
- **TYPE_CHECKING**

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/primir_quantification.py#L190"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_arguments`

```python
parse_arguments()
```

Command-line arguments parser.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/primir_quantification.py#L263"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_contribution`

```python
get_contribution(
    query_id: str,
    collapsed: bool = False,
    nh: bool = False
) → float
```

Return the contribution of a single alignment.

The read name format depends on the flags:

- collapsed + NH: `READ-#reads_NH`
- collapsed only: `READ-#reads`
- NH only: `READ_NH`
- neither: `READ`


**Arguments:**

 - <b>`query_id`</b>:  Read/query name from the INTERSECT record.
 - <b>`collapsed`</b>:  If `True`, parse #reads from the read name.
 - <b>`nh`</b>:  If `True`, parse NH from the read name.



**Returns:**
 Contribution as (#reads / NH) following the rules above.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/primir_quantification.py#L349"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_initial_data`

```python
get_initial_data(name: str, feat_extension: bool) → list[str]
```

Get the feature name and its extension.


**Arguments:**

- <b>`name`</b>: string with the feature name that can or not include its
  annotation extension in the format `name_5-extension_3-extension`
- <b>`feat_extension`</b>: specify whether the feature annotation extension
  has to be a field  on the final output



**Returns:**
  list with the feature name to be found in the final table and the number of
  extended positions (if asked for)


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/primir_quantification.py#L375"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main`

```python
main(args) → None
```

Tabulate `bedtools intersect -wo -s` output file.


---
