<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/blocksort.sh#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>script</kbd> `blocksort.sh`

Sort oligomap alignments based on their numerical names.

Usage
-----

```bash
bash blocksort.sh [input_file number_of_threads output_file | -h]
```

Arguments
---------

- <b>`$1`</b> | <b>`input_file`</b> (string): path to the input file with the
  oligomap alignments
- <b>`$2`</b> | <b>`number_of_threads`</b> (int): number of threads to run the
  sorting with
- <b>`$3`</b> | <b>`output_file`</b> (string): path to the output file where to
  write the sorted alignments

Options
-------

- <b>`-h`</b> | <b>`--help`</b>: Display help.

Exit codes
-----------

- <b>`0`</b>: If successful
- <b>`1`</b>: If input is wrong. Either the amount, order, or type

---
