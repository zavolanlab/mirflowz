<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_quantification.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `mirna_quantification.py`

Quantify miRNAs and corresponding isomiRs.

Read the input SAM file, calculate the contribution sum of the intersecting
alignments for each feature and write the result in tab-delimited format. The
following sections will cover what the SAM file must be like, how the counting
is done, and how the output can look like. The final section contain examples
covering the main modes.

**EXPECTED INPUT**

The SAM file must contain only the alignments that intersect with canonical
miRNAs and/or isomiRs. For each alignment, the intersecting feature name(s)
must be stored in a tag. If the tag determined by the CLI argument `--tag` is
empty, the alignment is ignored. In addition, the SAM file must be sorted by
the tag storing the feature names. Optionally, the SAM file can contain
collapsed reads, as, e.g., produced by the `fastx_collapser` tool of the
[FASTX-Toolkit][docs-fastx]. Finally, read names can also contain the alignment
NH vale (see **READ NAME FORMAT** section).

**FEATURE NAME FORMAT**

The name of the intersecting feature(s) has to follow the format:

    feat_name|5p-shift|3p-shift|CIGAR|MD|READ_SEQ

- For isomiRs, the whole name is used
- For canonical miRNAs, the `feat_name` is used

A feature is classified as canonical if the 5p and 3p shifts are both 0 and
the CIGAR and MD strings are the same. If there are different features names
separated by a semi-colon, they are treated as a unique isomiR name. For example,
the feature name `feat_name_1|0|0|23M|23|read_seq` would represent a canonical
miRNA and the feature names `feat_name_2|0|0|22M|0A0A2T17|read_seq` and
`feat_name_3|0|0|22M|22|read_seq;feat_name_4|1|0|23M|17A3T0C0|read_seq` would
represent isomiRs.


**READ NAME FORMAT**

The read name must have one of the following formats:

- `NAME`
- `NAME_NH`
- `NAME-COUNT`
- `NAME-COUNT_NH`

where:

- `NAME` is an arbitrary unique name
- `COUNT` is the number of identical reads that were collapsed
- `NH` is the alignment NH tag.

For example, the query name `my_query_name-13_2` would indicate that the read
represented by this alignment occurred 13 times before collapsing and the
alignment's NH value is 2.

Only if `COUNT` is in the name the CLI flag `--collapsed` can be set. Likewise,
the CLI flag `--nh` can only be set if NH is in the name.


**COUNT TABULATION**

The count for each feature is the sum of all the intersecting alignment's
contribution. The appropriate contribution is the ratio between the number of
reads and the alignment's NH.

The CLI flags `--collapsed` and `--nh` determine if these values are taken
from the alignment's name:

- If `--collapsed` is not set, the number of reads is considered to be 1
- If `--nh` is not set, the NH value is obtained from the NH tag
- If the NH tag is missing, its value is considered to be 1

**TABULATED FEATURES**

Which kind of features will appear in the output file is determined by the CLI
argument `--mir-list`. The three possible options are:

1. `mirna` to only include canonical miRNAs
2. `isomir` to only include isomiRs
3. `mirna` and `isomir` to include both canonical miRNAs and isomiRs


**OUTPUT FORMAT**

The output is tab-delimited file named `mir_counts_LIB`, where `LIB` is the
read's library name set in the CLI option `--lib`.

If the SAM file is empty, an empty file is produced.

By default each row contains the feature name and its partial count.

Three extra columns can be added by using the CLI flags `--count`, `--len`
and `--read-ids`:

- If `--count` is set, the output table will contain the number of best
  alignments for each feature
- If `--len` is set, the output table will contain the read's length
- If both `--count` and `--len` are set, the count will always be followed by
  the read's length
- If `--read-ids` is set, a semicolon-separated list of all alignment IDs
  that overlap with the feature will always be in the last column

Examples
--------

### Example 1 | Feature name stored in the tag `XN`

**Input:**

SAM with intersecting feature names in stored in the tag XN

**Command:**

    mirna_quantification.py SAM --tag XN

**Output:**

    hsa-miR-516b-5p                                     3.0
    hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG 0.6000000000000001


### Example 2 | Quantify mature miRNA

**Input:**

SAM meeting the minimal characteristics

**Command:**

    mirna_quantification.py SAM --mir-list mirna

**Output:**

    hsa-miR-517b-3p     1.3333333333333333
    hsa-miR-518b        1.0

### Example 3 | Quantify isomiR

**Input:**

SAM meeting the minimal characteristics

**Command:**

    mirna_quantification.py SAM --mir-list isomir

**Output:**

    hsa-miR-512-3p|0|1|23M|22C|AAGTGCTGTCATAGCTGAGGTAA  1.6666666666666665
    hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG 0.6000000000000001


### Example 4 | Collapsed reads with `NH` tag in the name

**Input:**

SAM with read names following the format `NAME-COUNT_NH`

**Command:**

    mirna_quantification.py SAM --collpased --nh

**Output:**

    hsa-miR-516b-5p                                     3.0
    hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG 0.6000000000000001


### Example 5 | Collapsed reads

**Input:**

SAM with read names following the format `NAME-COUNT`

**Command:**

    mirna_quantification.py SAM --collapsed

**Output:**

    hsa-miR-516b-5p                                     3.0
    hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG 0.6000000000000001


### Example 6 | Reads with the `NH` tag in their name

**Input:**

SAM with read names following the format `NAME_NH`

**Command:**

    mirna_quantification.py SAM --nh

**Output:**

    hsa-miR-516b-5p                                     3.0
    hsa-miR-517-5p|-1|0|23M|22T|ACCTCTAGATGGAAGCACTGTCG 0.6000000000000001


### Example 7 | Quantify mature miRNAs and add read IDs to the final table

**Input:**

SAM with read names following the format `NAME_NH`

**Command:**

    mirna_quantification.py SAM --read-ids --nh --mir-list mirna

**Output:**

    hsa-miR-498-5p      4.333333333333333   270396_3
    hsa-miR-1323        12.0                673650_2;673650_2


### Example 8 | Quantify isomiRs, add isomiR length and absolute counts

**Input:**

SAM meeting the minimal characteristics

**Command:**

    mirna_quantification.py SAM --count --len --mir-list isomir

**Output:**

    hsa-miR-512-3p|0|1|23M|22C0|AAGTGCTGTCATAGCTGAGGTAA     4.333333333333333  6 23
    hsa-miR-512-3p|0|1|23M|3T18C0|AAGGGCTGTCATAGCTGAGGTAA   12.0               8 19


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_quantification.py#L132"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_arguments`

```python
parse_arguments()
```

Command-line arguments parser.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_quantification.py#L246"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `collapsed_nh_contribution`

```python
collapsed_nh_contribution(aln: AlignedSegment) → float
```

Get the contribution of the alignment to the overall count.

The contribution is computed as the ratio of the number of reads collapsed in
the alignment and the NH value. It is assumed that the alignment query name
contains the number of collapsed reads as well as the NH value in the format
`NAME-COUNT_NH`.


**Arguments:**

- <b>`aln`</b>:  Alignment to which the overall contribution is calculated



**Returns:**
  Contribution of alignment to overall count

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_quantification.py#L279"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `collapsed_contribution`

```python
collapsed_contribution(aln: AlignedSegment) → float
```

Get the contribution of the alignment to the overall count.

The contribution is computed as the ratio of the number of reads collapsed in
the alignment and the value stored in the NH tag. If the tag is missing, the NH
value is 1. It is assumed that the alignment query name contains the number of
collapsed reads in the format `NAME-COUNT`.


**Arguments:**

- <b>`aln`</b>:  Alignment to which the overall contribution is calculated



**Returns:**
  Contribution of alignment to overall count


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_quantification.py#L316"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `nh_contribution`

```python
nh_contribution(aln: AlignedSegment) → float
```

Get the contribution of the alignment to the overall count.

The contribution is computed as the ratio of the number of reads collapsed in
the alignment and the value stored in the NH tag. If the tag is missing, the NH
value is 1. It is assumed that the alignment query name contains the NH value in
the format `NAME_NH`.


**Arguments:**

- <b>`aln`</b>:  Alignment to which the overall contribution is calculated


**Returns:**
  Contribution of alignment to overall count


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_quantification.py#L347"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `contribution`

```python
contribution(aln: AlignedSegment) → float
```

Get the contribution of the alignment to the overall count.

The contribution is computed as the ratio of the number of reads collapsed in
the alignment and the value stored in the NH tag. If the tag is missing, the
overall contribution is 1.


**Argumens:**

- <b>`aln`</b>:  Alignment to which the overall contribution is calculated



**Returns:**
  Contribution of alignment to overall count


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_quantification.py#L368"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_name`

```python
get_name(pre_name: str) → list[str]
```

Get the final name for the species name.

Take a string and processes it to obtain the final name for the species and the
type of miRNA the string belongs to. Only the `feat_name` is returned if the
3p-shift and 5p-shift are 0 and the CIGAR and MD are the same. Otherwise, the
whole input string is returned.


**Arguments:**

- <b>`pre_name`</b>:  string with the format
  `feat_name|5-shift|3-shift|CIGAR|MD|read_seq`


**Returns:**
  list with the species name to be found in the final table and its type


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_quantification.py#L393"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `write_output`

```python
write_output(
    name: str,
    species: list[str],
    mir_list: list[str],
    mirna_out: Path
) → None
```

Write to the output the correct miRNA type.

**Arguments**:

- <b>`name`</b>: miRNA species type
- <b>`species`</b>: List of candidate miRNA species to write in the output file
- <b>`mir_list`</b>: List of miRNA types to have in the output table
- <b>`mirna_out`</b>: Path to the output file

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/mirna_quantification.py#L404"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main`

```python
main(args) → None
```

Quantify miRNAs and corresponding isomiRs.


---
