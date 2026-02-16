<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/map_chromosomes.pl#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>script</kbd> `map_chromosomes.pl`

Map/rename chromosome identifiers in a delimited text file using a tab-delimited
mapping table.

Reads a delimited text file, replaces the chromosome identifier in a
user-selected column according to a mapping table, and writes the transformed
lines to an output file.

The mapping table is read from a tab-delimited file. Lines starting with `#`
are ignored. Each mapping line is split on `TAB`; the first field is
treated as the source chromosome name and the second field as the target name.

If a mapping line contains only one field, the target is treated as `<remove>`,
which causes matching lines to be dropped.

Lines whose chromosome value is not present in the mapping table are removed
(with a message to `STDOUT`). Header / comment lines (first field starts with
`#`) are preserved in the output.

Output behavior
---------------

- If the chromosome field is defined and exists in the mapping table: it is
  replaced.
- If the mapped value is `<remove>`: the line is dropped and a message is
  printed to `STDOUT`.
- If the first field starts with `#`: the line is written unchanged to the
  output.
- Otherwise: the line is dropped and a message is printed to `STDOUT`
  indicating an invalid/unmapped chromosome.

Usage
-----

```bash
perl map_chromosomes.pl [INPUT] [COL] [DELIMITER] [MAP] [OUTPUT]
```

Arguments (positional)
----------------------

The script expects exactly 5 arguments in a specific order:

1. <b>`[INPUT]`</b>: Path to the input text file to be processed
2. <b>`[COL]`</b>: 1-based column index indicating which column contains the
   chromosome identifier to map. Whitespace is stripped. Internally, the script
   converts this to 0-based indexing. Must be numeric (validated via
  `looks_like_number`)
3. <b>`[DELIMITER]`</b>: Delimiter identifier for splitting and re-joining
   columns. Must be one of: `TAB`, `COMMA`, `DASH`, `UNDERSCORE`, `PIPE`,
  `DOT`, `SPACE`. These correspond to the split/join characters used by the
  script
4. <b>`[MAP]`</b>: Tab-delimited mapping file. Each non-comment line is split
  on `TAB`: `<from>\t<to>`
5. <b>`[OUTPUT]`</b>: Output file path. The script writes transformed lines to
  this file

Requirements
------------

- <b>Perl version</b>: `>= 5.40.2`
- <b>Modules:</b>
    - <b>`Scalar::Util`</b>: `>= 5.42.0`

---
