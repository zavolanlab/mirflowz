<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/nh_filter.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `nh_filter.py`

Filter alignments in a SAM file by NH tag.


For each alignment, check its NH tag, and if the value is higher than the one
specified in `--max_nh`, the aligned read is removed.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/nh_filter.py#L10"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_arguments`

```python
parse_arguments()
```

Parse command-line arguments.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/nh_filter.py#L51"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main`

```python
main(arguments) → None
```

Filter alignments by its NH tag value.

---
