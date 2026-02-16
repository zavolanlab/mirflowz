<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/annotate_sam_with_intersecting_features.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `annotate_sam_with_intersecting_features.py`

Annotate SAM alignments with their intersecting feature(s).

Add a custom tag ("YW") to each alignment in the SAM file if an intersecting
feature is found in the INTERSECT file. The INTERSECT file must be the result
of the call `bedtools intersect -wo -s -a GFF3/GTF -b BAM`. If either the
INTERSECT or the SAM file is empty, only the SAM file header is returned.

Each matching feature is used to build a tag with the following format:

 FEATURE_ID|5p-shift|3p-shift|CIGAR|MD|READ_SEQ

Where:

- `FEATURE_ID`: Extracted from the specified attribute in the INTERSECT file
  (default: "name")
- `5p-shift`: Difference between the (possibly adjusted) feature start and
  the alignment start
- `3p-shift`: Difference between the alignment end and the (possibly  adjusted)
  feature end
- `CIGAR` and `MD`: The alignment's CIGAR string and MD tag respectively
- `READ_SEQ`: The read sequence from the alignment as stored in the SAM file

Optional adjustment:

`--extension`: Adjust the feature's start and end coordinates by the given value
(start is increased and end decreased). In addition, requires both the 5' and
3' shift values to be within +/- this value to include the feature tag in the
alignment.

If an alignment has multiple intersecting features, the tag values are
concatenated using a semicolon as the separator. If there are no intersecting
features or none of the intersecting features pass the shift filter, the
alignment is skipped.

**Examples**
--------

### Example 1 | Feature intersects alignment; coordinates adjustment and shift allowed

**Use case:**

Prior to checking if the feature is intersecting the alignment, its
coordinates are adjusted by the value specified in `--extension`. In
addition, the same value is used to specify the +/- shift allowed between
the feature and the read alignment start and end coordinates.

**Command:**

    annotate_sam_with_intersecting_features.py -i INTERSECT -s SAM --extension 5

**IN INTERSECT record:**

    19      .       miRNA   5332    5365    .       +       .       ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786     19      5337    5358    read_1  255     +       21

**IN SAM record:**

    read_1  0       19      5338    255     21M     *       0       0       TCAAAACTGAGGGGCATTTTC   *       MD:Z:21 NH:i:1  NM:i:0

**Intersection before coordinates adjustment visualization:**


        ---|===============================|--- (feature)
        ---------|===================|--------- (read)

**Intersection after coordinates adjustment visualization:**


        --------|=====================|-------- (feature)
        ---------|===================|--------- (read)

**OUT SAM record:**

    read_1  0       19      5338    255     21M     *       0       0       TCAAAACTGAGGGGCATTTTC   *       MD:Z:21 NH:i:1  NM:i:0  YW:Z:hsa-miR-1323|1|-1|21M|21|TCAAAACTGAGGGGCATTTTC

**Description:**

The feature start and end coordinates after the adjustment are 5337 and 5360
respectively. The read alignment starts at position 5338. As the read has length
21, its end position is 5359. The feature intersects the read alignment with an
overhang within the specified shift range (+/- 5) so it is added as a new tag in
the output SAM record.



### Example 2 | Feature intersects alignment; no coordinates adjustment or shift allowed

**Use case:**

The feature and read alignment coordinates must perfectly match.

**Command:**

    annotate_sam_with_intersecting_features.py -i INTERSECT -s SAM

**IN INTERSECT record:**

    19      .       miRNA   5338    5359    .       +       .       ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786     19      5337    5358    read_2  255     +       21

**IN SAM record:**

    read_2  0       19      5338    255     21M     *       0       0       TCAAAACTGAGGGGCATTTTC   *       MD:Z:21 NH:i:1  NM:i:0

**Intersection visualization:**


        ---------|===================|--------- (feature)
        ---------|===================|--------- (read)

**OUT SAM record:**

    read_2  0       19      5338    255     21M     *       0       0       TCAAAACTGAGGGGCATTTTC   *       MD:Z:21 NH:i:1  NM:i:0  YW:Z:hsa-miR-1323|0|0|21M|21|TCAAAACTGAGGGGCATTTTC

**Description:**

The feature start and end coordinates are 5338 and 5359 respectively. The read
alignment starts at position 5338. As the read has length 21, its end position
is 5359. The feature perfectly intersects the read alignment so it is added as
a new tag in the output SAM record.


### Example 3 | Non-intersecting feature; shift filter not passed

**Use case:**

Prior to checking if the feature is intersecting the alignment, its coordinates
are adjusted by the value specified in `--extension`. In  addition, the same
value is used to specify the +/- shift allowed  between the feature and the
read alignment start and end coordinates.

**Command:**

    annotate_sam_with_intersecting_features.py -i INTERSECT -s SAM --extension 1

**IN INTERSECT record:**

    19      .       miRNA   5332    5365    .       +       .       ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786     19      5337    5358    read_3  255     +       21

**IN SAM record:**

    read_3  0       19      5338    255     21M     *       0       0       TCAAAACTGAGGGGCATTTTC   *       MD:Z:21 NH:i:1  NM:i:0

**Intersection before coordinates adjustment visualization:**


        ---|===============================|--- (feature)
        ---------|===================|--------- (read)

**Intersection after coordinates adjustment visualization:**


        ----|=============================|---- (feature)
        ---------|===================|--------- (read)

**Description:**

The feature start and end coordinates after the adjustment are 5333 and 5364
respectively. The read alignment starts at position 5338. As the read has
length 21, its end position is 5359. There is a 5-nucleotide overhang on both
ends. Thus, the feature is  not considered to intersect the read alignment and
the alignment is  not written in the output file.



### Example 4 | Feature intersects alignment; using feature's "Alias"

**Use case:**

The feature and read alignment coordinates must perfectly match.

**Command:**

    annotate_sam_with_intersecting_features.py -i INTERSECT -s SAM --id alias

**IN INTERSECT record:**

    19      .       miRNA   5338    5359    .       +       .       ID=MIMAT0005795;Alias=MIMAT0005795;Name=hsa-miR-1323;Derives_from=MI0003786     19      5337    5358    read_4  255     +       21

**IN SAM record:**

    read_4  0       19      5338    255     21M     *       0       0       TCAAAACTGAGGGGCATTTTC   *       MD:Z:21 NH:i:1  NM:i:0

**Intersection visualization:**


        ---------|===================|--------- (feature)
        ---------|===================|--------- (read)

**OUT SAM record:**

    read_4  0       19      5338    255     21M     *       0       0       TCAAAACTGAGGGGCATTTTC   *       MD:Z:21 NH:i:1  NM:i:0  YW:Z:MIMAT0005795|0|0|21M|21|TCAAAACTGAGGGGCATTTTC

**Description:**

The feature start and end coordinates are 5338 and 5359 respectively. The read
alignment starts at position 5338. As the read has length 21, its end position
is 5359. The feature perfectly intersects the read alignment so it is added as
a new tag in the output SAM record. In this case, instead of using the feature
`Name` (default), the `Alias` is used.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/annotate_sam_with_intersecting_features.py#L203"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_arguments`

```python
parse_arguments()
```

Command-line arguments parser.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/annotate_sam_with_intersecting_features.py#L261"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_intersect_output`

```python
parse_intersect_output(
    intersect_file: Path,
    feat_id: str = 'name',
    extension: int = 0
) → Optional[Dict[str, List[Tuple[str, int, int]]]]
```

Parse `bedtools intersect -wo -s` output file.

Given an INTERSECT file generated by intersecting a GFF3/GTF file (`-a`) with a
BAM file (`-b`) using bedtools intersect, create a dictionary where the
alignment names are the keys. The values are lists containing the feature name,
start and end positions.

The `feat_id` argument specifies the feature name to use, and the extension
argument adjusts the feature coordinates by adding the given value and
subtracts it from the end position.

If the INTERSECT file is empty, `None` is returned.


**Arguments:**

- <b>`intersect_file`</b>:  Path to the INTERSECT file.
- <b>`feat_id`</b>:  ID used to identify the feature. Defaults to "name".
- <b>`extension`</b>:  Number of nucleotides the start and end coordinates have
  to be adjusted. Defaults to 0.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/annotate_sam_with_intersecting_features.py#L298"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_tags`

```python
get_tags(
    intersecting_feat: list,
    alignment: AlignedSegment,
    extend: int = 0
) → set
```

Construct a custom tag for an alignment based on intersecting features.

For each intersecting feature (given as a tuple of identifier, start, and end),
compute the 5'-shift and 3'-shift values relative to the alignment.

Tag format:

    FEATURE-ID|5p-shift|3p-shift|CIGAR|MD|READ_SEQ

Only features where both shift values are within +/- `extend` are included.
If multiple features are valid, the resulting tags are collected in a set.

Coordinate conventions
----------------------

- `feat_start` / `feat_end` (from INTERSECT annotations) are 1-based, inclusive
- `pysam` exposes `alignment.reference_start` as 0 based, inclusive and
  `alignment.reference_end` as 0-based, exclusive

To compare these consistently, the feature coordinates are implicitly converted
to 0-based half-open:

    new_feat_start = feat_start - 1
    new_feat_end = feat_end

Therefore, the implemented shifts are:

    shift_5p = alignment.reference_start - (feat_start - 1)  = alignment.reference_start - feat_start + 1
    shift_3p = alignment.reference_end - feat_end


**Arguments:**

 - <b>`intersecting_feat`</b>:  list of tuples (feature identifier, start and end).
 - <b>`alignment`</b>:  a pysam.AlignedSegment object
 - <b>`extend`</b>:  maximum allowed absolute shift for both ends (default: 0)


**Returns:**
 A set of tag strings constructed from valid intersecting features.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/annotate_sam_with_intersecting_features.py#L361"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main`

```python
main(args) → None
```

Annotate alignments in a SAM file with intersecting feature tags.

---
