<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/oligomap_output_to_sam_nh_filtered.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `oligomap_output_to_sam_nh_filtered.py`

Transform oligomap output FASTA file to SAM keeping the best alignments.

Read the input file, pick the best alignment(s) for each read (_i.e_, all of
the alignments have either 0 or 1 error). In addition, if the `--nh-filter`
CLI argument is set, filter the reads with more hits than the number provided.
The following sections will cover what the input file must be like and how the
output can look like.

**EXPECTED INPUT**

The FASTA file must be the output of mapping your library with the tool
[oligomap][oligomap]. Refer to the **Output format** section in its main
`README.md` for more information.

In addition, alignments have to be sorted by the read name. An example of how
an entry would look like in the **EXAMPLE** section.


**OUTPUT FORMAT**

The output consist on the filtered set of alignments from the input file in
SAM format. Only the best alignments per read (_i.e._, either all the alignments
have 0 or 1 error) are written to the standard output. Moreover, if the
`--nh-filter` CLI argument is given, reads with more hits than the provided
value are discarded. If none of the alignments meet the criteria, nothing is
returned. The fields in the SAM entry are:

    field | value
    --------------
    QNAME | Read's name
    FLAG  | Set to 0 if the reference sequence is found in the positive strand,
          | and to 16 otherwise.
    RNAME | Reference sequence name
    POS   | Alignment's first position in the reference sequence
    MAPQ  | Value set to 255 (mapping quality not available)
    CIGAR | Alignment's CIGAR string
    RNEXT | Value set to * (unavailable information)
    PNEXT | Value set to 0 (unavailable information)
    TLEN  | Value set to 0 (unavailable information)
    SEQ   | Mapped sequence
    QUAL  | Value set to * (Phred-scaled base quality not stored)
    NM:i: | Alignment's edit distance (Optional field)
    MD:Z: | Alignment's MD string (Optional field)
    NH:i: | Alignment's NH (Optional field)

Example
-------

**Input format:**

    82-1 (23 nc) 1..23      19      44377..44398
    19
    errors: 1 orientation: +
    CTACAAAGGGAAGCACTTGTCTC
    |||||||||||||||||| ||||
    CTACAAAGGGAAGCACTT-TCTC



    97-1 (22 nc) 1..22      19      5338..5359
    19
    errors: 1 orientation: +
    TCAAAACTGAGGTGCATTTTCT
    |||||||||||| |||||||||
    TCAAAACTGAGGGGCATTTTCT

**Output format:**

    82-1    0       19      44377   255     18M1I4M *       0       0       CTACAAAGGGAAGCACTTGTCTC *       NM:i:1  MD:Z:23         NH:i:1
    97-1    0       19      5338    255     22M     *       0       0       TCAAAACTGAGGTGCATTTTCT  *       NM:i:1  MD:Z:12G9       NH:i:1



Paula Iborra. Zavolan Lab.
Adapted version of Alessandro Crippa script.
Refactored and documented by Iris Mestres-Pascual.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/oligomap_output_to_sam_nh_filtered.py#L97"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_arguments`

```python
parse_arguments()
```

Command-line arguments parser.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/oligomap_output_to_sam_nh_filtered.py#L127"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_cigar_md`

```python
get_cigar_md(
    errors: str,
    sequence: str,
    bars_line: str,
    ref_seq: str
) → tuple[str, str]
```

Get the CIGAR and MD strings.

Given the read and target sequences, the number of errors and the alignment's
bar representation, this function first checks if there is any error in the
mapping. If there is no errors, the CIGAR and MD strings are the read's length
followed by and "M" and the read's length respectively. If there is an error,
the function checks if the read has a deletion (the read sequence contains a
dash), an insertion (the reference sequence contains a dash), or a single-point
mutation. Finally, it checks the error's position in the read sequence and builds
the CIGAR and MD strings.


**Arguments:**

- <b>`errors`</b>: number of mapping errors (either 0 or 1)
- <b>`sequence`</b>:  read sequence
- <b>`bars_line`</b>:  alignment's bar representation
- <b>`ref_seq`</b>:  reference sequence


**Returns:**

- <b>`cigarStr`</b>:  alignment's CIGAR string
- <b>`matchingString`</b>:  alignment's MD string with its corresponding tag


Examples
---------

### Example 1 | Perfect match

    read sequence ->               GAAGGCGCTTCACCTTTGGAGT
    alignment representation ->    ||||||||||||||||||||||
    reference sequence ->          GAAGGCGCTTCACCTTTGGAGT
    errors -> 0


 - `cigarStr`: 21M
 - `matchingString`: MD:Z:21


### Example 2 | One insertion

    read sequence ->               GAAGGCGCTTCACCTTTGGAGT
    alignment representation ->    ||||||||||| ||||||||||
    reference sequence ->          GAAGGCGCTTC-CCTTTGGAGT
    errors -> 1

- `cigarStr`: 11M1I10M
- `matchingString`: MD:Z:21


### Example 3 | One deletion

    read sequence ->               GAAGGCGCTTC-CCTTTGGAGT
    alignment representation ->    ||||||||||| ||||||||||
    reference sequence ->          GAAGGCGCTTCCCCTTTGGAGT
    errors -> 1

- `cigarStr`: 11M1D10M
- `matchingString`: MD:Z:11C10


### Example 4 | One mismatch

    read sequence ->               GAAGGCGCTTCACCTTTGGAGT
    alignment representation ->    ||||||||||| ||||||||||
    reference sequence ->          GAAGGCGCTTCCCCTTTGGAGT
    errors -> 1

- `cigarStr`: 21M
- `matchingString`: MD:Z:11C10


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/oligomap_output_to_sam_nh_filtered.py#L250"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_sam_fields`

```python
get_sam_fields(aln: list[str]) → Fields
```

Create the read's alignment in SAM format.

Given the set of lines oligomap provides as an alignment representation,
this function returns a `Fields` class object with the mandatory fields in a
SAM file and three optional ones.

    field | value
    --------------
    QNAME | Read's name
    FLAG  | Set to 0 if the reference sequence is found in the positive strand,
          | and to 16 otherwise.
    RNAME | Reference sequence name
    POS   | Alignment's first position in the reference sequence
    MAPQ  | Value set to 255 (mapping quality not available)
    CIGAR | Alignment's CIGAR string
    RNEXT | Value set to * (unavailable information)
    PNEXT | Value set to 0 (unavailable information)
    TLEN  | Value set to 0 (unavailable information)
    SEQ   | Mapped sequence
    QUAL  | Value set to * (Phred-scaled base quality not stored)
    NM:i: | Alignment's edit distance (Optional field)
    MD:Z: | Alignment's MD string (Optional field)
    NH:i: | Alignment's NH (Optional field)


**Arguments:**

- <b>`aln`</b>: list made up of the lines that form the output alignment after
  mapping with oligomap.



**Returns:**

- <b>`fields`</b>:  `Fields` class `NamedTuple` with the alignment's data as
  SAM fields.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/oligomap_output_to_sam_nh_filtered.py#L309"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `eval_aln`

```python
eval_aln(
    nhfilter: int,
    d: Dict[str, list],
    min_err_nh: Dict[str, list],
    fields: Fields
) → None
```

Evaluate an alignment to add, discard or write it to the `STDOUT`.

Given a read's alignment, this function first checks if the dictionary storing
the read alignments is empty. If it is, the function checks if the read has an
entry in the dictionary storing the current minimum error and the NH. If there
is not, or the minimum number of errors is higher than the alignment under
evaluation, the entry for that read in both dictionaries is set to the data of
the alignment under evaluation.

If the alignments dictionary is not empty, the function checks if the read
under evaluation is the first in the dictionary. If it is, and the number of
errors in the current alignment is the same as the minimum error for that read,
the number of hits is increased by 1. Then:

- If the number of hits does not exceed the maximum NH or the NH is not provided,
  the alignment under evaluation is appended to the read's entry.
- If a maximum NH is set, and the read exceeds it, it is removed.
- If the minimum number of errors is higher than the number of errors of the
  alignment under evaluation, the entry for that read in both dictionaries is
  set to the data of the current alignment.

If the read is not the first in the dictionary, the function writes all the
alignments for the read in the first position to the standard output, and both
dictionaries are reset to only include the data for the alignments under
evaluation.


**Arguments:**

- <b>`nhfilter`</b>:  maximum number of hits an read can have to be kept
- <b>`d`</b>:  dictionary with read names as keys and a list with `Fields`
  class `NamedTuples` as values
- <b>`min_err_nh`</b>:  dictionary with read names as keys and a list with the
  current minimum error and the number of hits in this order as values
- <b>`fields`</b>:  read's alignment as a `Fields` class `NamedTuple`


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/oligomap_output_to_sam_nh_filtered.py#L413"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main`

```python
main(arguments) → None
```

Convert the alignments in the oligomap output file to SAM format.


---

## <kbd>class</kbd> `Fields`
Class to store an alignment in its different SAM fields.


**Attributes:**

- <b>`read_name`</b>: read's name
- <b>`flags`</b>: bitwise SAM flags
- <b>`ref_seq_name`</b>: reference sequence name
- <b>`position_in_ref`</b>: 1-based leftmost mapping position
- <b>`map_q`</b>: mapping quality
- <b>`cigar_str`</b>: CIGAR string
- <b>`r_next`</b>: reference name of the mate/next read
- <b>`p_next`</b>: position of the mate/next read
- <b>`t_len`</b>: observed template length
- <b>`sequence`</b>: segment sequence
- <b>`qual`</b>: ASCII of Phred-scaled base quality + 33
- <b>`edit_dist`</b>: Alignment's edit distance
- <b>`md_str`</b>: MD string


---
