<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/filter_multimappers.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `filter_multimappers.py`

Filter miRNA reads mapped to multiple locations by indel count.

From all the alignments with the same query name and edit distance, keep the
ones with the higher amount of indels. Supplementary alignments are dismissed
as they are part of a chimeric alignment, which is composed of multiple linear
alignments with minimal overlap.

Additionally, the `NH` and `HI` tags are updated to match the new amount of
alignments and keep their identifier within the new set respectively.

If the CLI flag `--nh` is set, query names are appended the suffix `_#` were `#`
is the updated alignment's NH tag.

The following assumptions are made:

- The input SAM file is sorted by query name.



Examples
---------

### Example 1 | Different number of indels

**IN SAM records:**

    read-1  0       19      77595   255     8M1D14M *       0       0       CTGACATCAGTGATTCTCCTGC  *       MD:Z:3G1T2^A14  NH:i:2  NM:i:3  XA:Z:Q  XI:i:1
    read-1  0       19      330456  255     4M1D1M1I3M1D13M *       0       0       CTGACATCAGTGATTCTCCTGC  *       MD:Z:4^G4^A13   NH:i:2  NM:i:3  XA:Z:Q  XI:i:0

**Alignments:**

    CTGACATC-AGTGATTCTCCTGC
    ||| | || |||||||||||||| (1 indel, 2 mismatches, discarded)
    CTGGCTTCAAGTGATTCTCCTGC

    CTGA-CATCA-GTGATTCTCCTGC
    |||| | ||| ||||||||||||| (3 indels, 0 mismatches, retained)
    CTGAGC-TCAAGTGATTCTCCTGC

**Command:**

    filter_multimappers.py SAM > out_SAM

**OUT SAM record:**

    read-1  0       19      330456  255     4M1D1M1I3M1D13M *       0       0       CTGACATCAGTGATTCTCCTGC  *       MD:Z:4^G4^A13   NH:i:1  HI:i:1  NM:i:3  XA:Z:Q  XI:i:0


### Example 2 | Equal number of indels

**IN SAM records**:

    read-2  0       19      142777  255     5M1I15M *       0       0       GCTTCAAGCCTCCCACCTAGC   *       MD:Z:14A0G4     NH:i:3  NM:i:3  XA:Z:Q  XI:i:0
    read-2  0       19      270081  255     6M1I14M *       0       0       GCTTCAAGCCTCCCACCTAGC   *       MD:Z:14G0G4     NH:i:3  NM:i:3  XA:Z:Q  XI:i:2
    read-2  0       19      545543  255     6M1I14M *       0       0       GCTTCAAGCCTCCCACCTAGC   *       MD:Z:14A0G4     NH:i:3  NM:i:3  XA:Z:Q  XI:i:1

**Alignments:**

    GCTTCAAGCCTCCCACCTAGC
    ||||| |||||||||  |||| (1 Indel, 2 mismatches, retained)
    GCTTC-AGCCTCCCAAGTAGC

    GCTTCAAGCCTCCCACCTAGC
    |||||| ||||||||  |||| (1 Indel, 2 mismatches, retained)
    GCTTCA-GCCTCCCAGGTAGC

    GCTTCAAGCCTCCCACCTAGC
    |||||| ||||||||  |||| (1 Indel, 2 mismatches, retained)
    GCTTCA-GCCTCCCAAGTAGC

**Command:**

    filter_multimappers.py SAM > out_SAM

**OUT SAM records:**

    read-2_3        0       19      142777  255     5M1D15M *       0       0       GCTTCAAGCCTCCCACCTAGC   *       MD:Z:14A0G4     NH:i:3  HI:i:1  NM:i:3  XA:Z:Q  XI:i:0
    read-2_3        0       19      270081  255     6M1I14M *       0       0       GCTTCAAGCCTCCCACCTAGC   *       MD:Z:14G0G4     NH:i:3  HI:i:2  NM:i:3  XA:Z:Q  XI:i:2
    read-2_3        0       19      545543  255     6M1I14M *       0       0       GCTTCAAGCCTCCCACCTAGC   *       MD:Z:14A0G4     NH:i:3  HI:i:3  NM:i:3  XA:Z:Q  XI:i:1


### Example 3 | Add NH tag as read's name suffix

**IN SAM record:**

    read-3  0       19      5338    1       15M1I7M *       0       0       TCAAAACTGAGGGGCTATTTTCT *       HI:i:1  NH:i:1  NM:i:1  MD:Z:22 RG:Z:A1 YZ:Z:0

**Alignments:**

    TCAAAACTGAGGGGCTATTTTCT
    ||||||||||||||| ||||||| (1 Indel, 0 mismatches, retained)
    TCAAAACTGAGGGGC-ATTTTCT

**Command:**

    filter_multimappers.py SAM --nh > out_SAM

**OUT SAM record:**

    read-3_1        0       19      5338    1       15M1I7M *       0       0       TCAAAACTGAGGGGCTATTTTCT *       HI:i:1  NH:i:1  NM:i:1  MD:Z:22 RG:Z:A1 YZ:Z:0


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/filter_multimappers.py#L95"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_arguments`

```python
parse_arguments()
```

Command-line arguments parser.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/filter_multimappers.py#L129"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `count_indels`

```python
count_indels(aln: AlignedSegment) → int
```

Count the number of indels in an alignment based on its CIGAR string.

This function counts the number of indels in the alignment based on the
alignment CIGAR string returned by `.cigartuples`. Insertions are encoded as 1s
whereas the deletions are encoded as 2. Refer to the
[`.cigartuples` documentation][cigartuples] for more information.



**Arguments:**

 - <b>`aln`</b>:  alignment to count insertions and deletions for.



**Returns:**

 - <b>`int`</b>:  sum of insertions and deletions in that alignment.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/filter_multimappers.py#L150"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `find_best_alignments`

```python
find_best_alignments(
    alns: List[AlignedSegment],
    nh: bool = False
) → List[AlignedSegment]
```

Find alignments with more indels.

Using the function `count_indels` each alignment is assigned with its number of
indels. Then, only those alignments with the higher amount of indels are
returned. In addition, the `NH` and `HI` tags are updated to match the new
amount of alignments and keep its identifier within the set respectively.

If `nh` is set to `True`, all query names are appended a suffix `_#` were `#`
is the new alignment's NH tag.


**Arguments:**

 - <b>`alns`</b>:  alignments with the same query name
 - <b>`nh`</b>:  whether to suffix read names with the alignment's NH tag



**Returns:**

 - <b>`best_alignments`</b>:  alignments with the more indels


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/filter_multimappers.py#L195"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `write_output`

```python
write_output(alns: List[AlignedSegment]) → None
```

Write the output to the standard output (stdout).



**Arguments:**

 - <b>`alns`</b>:  alignments with the same query name


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/filter_multimappers.py#L205"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main`

```python
main(args) → None
```

Filter multimappers by indels count.



---
