<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/merge_tables.R#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>script</kbd> `merge_tables.R`

Merge miRNAs quantification tables.

Usage
-----

```bash
Rscript merge_tables.R [--help] [--verbose] [OPTIONS] --input_dir <PATH>
```

Arguments
---------

- <b>`--input_dir=DIRECTORY` (required)</b>: Absolute path from where input
  files shall be read.

Options
-------

- <b>`--output_file=FILE`</b>: Path to the output file

  (default: working-directory/counts.tab).
- <b>`--prefix`</b>: Prefix for reading input files (default: `NULL`).
- <b>`-h`</b> | <b>`--help`</b>: Show this information and die.
- <b>`-u`</b> | <b>`--usage`</b>: Show this information and die.
- <b>`-v`</b> | <b>`--verbose`</b>: Print log messages to `STDOUT`.

Dependencies
------------

- <b>R version</b>: `>= 3.6.0`
- <b>R packages</b>:
    - <b>`optparse`</b>: `>= 1.6.2`
    - <b>`dplyr`</b>: `>= 1.1.4`

---


<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/merge_tables.R#L134"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `get_table`


```r
get_table <- function(tbl_pth, prefix)
```

Read and process input table

`get_table()` uses `tryCatch()` to read the file in `tbl_pth`. If the table
is empty and an error is raised, the returned data frame consist of one row
with a `NA` in both fields.

**Arguments:**

- <b>`tbl_pth`</b>: Path to the input table.
- <b>`prefix`</b>: String to be removed from the input file name. It must be
  present in all the tables to be merged.

**Returns:**

A data frame containing the miRNA species to be counted in first column, named
`ID`, and their counts in that file in the second one. The name of the second
column in the data frame is obtained by removing the `prefix` from the input
file name. If no `prefix` is given, the whole file name is used.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/merge_tables.R#L179"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `merge_tables`

```r
merge_tables <- function(cwd, prefix)
```

Merge tables with the same prefix

`merge_tables()` takes all the files in `cwd` that start with `prefix` and
merges them keeping all the miRNA species present in each of the tables.

The function `get_table()` is used to make sure that even if an empty input
file is given, the merge can still be done by creating a data frame with a
single row made of `MA`s. Therefore, prior to the returning of the merged table,
if there is a row with a `NA` in the `ID` filed, it is removed.

The function `dplyr::full_join()` is used for the merge. This implies that if
a miRNA species in `ID` is missing in any of the tables being joined, its value
is set to `NA` in that column.

**Arguments:**

- <b>`cwd`</b>: Path to the directory containing the input tables.
- <b>`prefix`</b>: String used in all the tables to be selected for the merge.
  If not provided, all the files in `cwd` are used.

**Returns:**

A single data frame, `mat`, with all the miRNA species present in the input
tables in the first column, `ID`, and their counts. Each input file has it own
column.

If all the input tables are empty, the output only consist of the table's
header, and if no files starting with `prefix` are found, nothing is returned.

---
