<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/get_lines_w_pattern.sh#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>script</kbd> `get_lines_w_pattern.sh`

Retrieve all lines with the matching pattern `--pattern` within the requested
column (`--column`).

Usage
-----

```bash
bash get_lines_w_pattern.sh [[[-f file ] [-c column] [-p pattern] [-o outname]] | [-h]]
```

Options
-------

- <b>`-h`</b> | <b>`--help`</b>: Display help.
- <b>`-f`</b> | <b>`--file`</b> (string): Path to the input file
- <b>`-o`</b> | <b>`--out`</b> (string): Path to the output file
- <b>`-c`</b> | <b>`--column`</b> (int): Column index to where to look for the
  pattern
- <b>`-p`</b> | <b>`--pattern`</b> (string): Character pattern to look for

Exit codes
-----------

- <b>`0`</b>: If successful
- <b>`1`</b>: If input is wrong. Either the amount or type

---
