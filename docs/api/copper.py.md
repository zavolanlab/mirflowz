<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `copper.py`

COlor-coded PileuP gEneratoR.

Generate color-coded pileups in HTML and their corresponding CSS style file
either for a whole directory or a single file.

If a directory is provided as input, an HTML file is generated for each
file with just one CSS file.

For a proper HTML creation, the counts column title must be specified in
the CLI argument `--counts_id` (see the "Constraints" section for a more
detailed explanation on the input format)

The final color-coded pileup has a two-line header:

    - The first line starts with a label used to categorize samples based
        on experimental context specified in `--group_id`, followed by the
        sample name.
    - The second line contains the mature miRNA the pileup is made for,
        followed by the +/- nucleotide shift allowed at either end of the
        miRNA sequence.

All the nucleotides outside the canonical sequence use the `--gap`
background. Therefore, only the positions aligning with the canonical
sequence are colored.

If the input pileup contains the canonical sequence, the output pileup
will display it in bold.


Constraints:
- The tabulated file has to follow the format specified [here][pileup-format].

The following assumptions are made:

- The input file name must be `LIB.MIRNA.#-shift.tab` where `LIB` is the
  library name, `MIRNA` is the mature miRNA the pileup is made for, and
  `#` is the maximum allowed shift on both ends of the canonical
  sequence.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L426"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `parse_arguments`

```python
parse_arguments()
```

Parse command-line arguments.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L599"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `create_dir_pileups`

```python
create_dir_pileups(in_dir: Path, out_dir: Path, obj: ColorPileup) → None
```

Create color-coded pileups for each file in a directory.

This function assumes that ASCII-style alignment pileup files have the
extension ".tab".


**Arguments:**

 - <b>`in_dir`</b>:  Path to the input directory where pileup files are stored.
 - <b>`out_dir`</b>:  Path to the directory where color-coded pileups are
   written.
 - <b>`obj`</b>:  Initialized `ColorPileup` class with the active color scheme.


---


<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L623"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>function</kbd> `main`

```python
main(args)
```

Generate color-coded HTML alignment pileups.


---

## <kbd>class</kbd> `ColorPileup`

Class for code-coloring pileups.


**Attributes:**

 - <b>`input`</b>:  Path to the ASCII-style alignment pileup file.
 - <b>`html`</b>:  Path to the HTML output color-coded pileup file.
 - <b>`css`</b>:  Path to the CSS output file.
 - <b>`group_id`</b>:  String to categorize samples based on experimental
   context.
 - <b>`counts`</b>:  String used as the counts column title.
 - <b>`keep_info`</b>:  Whether to keep the sequence representation name and
   genomic coordinates on the final display.
 - <b>`color_code`</b>:  Dictionary mapping pileup characters to `ColorValues`.


<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L104"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `__init__`

```python
__init__() → None
```

Initialize class.


---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L169"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `create_css`

```python
create_css() → None
```

Write the CSS stylesheet used by all generated pileup HTML files.

This method creates the stylesheet referenced by the HTML pileup output. The
stylesheet is composed of three sections:

1. Global layout rules for page headers and container elements.
2. Per-line formatting rules for each pileup row and its associated text.
3. Per-character styling rules for nucleotide, gap, and generic symbols.

The character-specific CSS classes are generated dynamically from
`self.color_code`. Each entry defines the background and text color for a
pileup character, allowing the HTML output to render color-coded sequence
positions consistently.

The final stylesheet is written to `self.css`.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L345"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `create_html`

```python
create_html() → None
```

Create the color-coded pileup HTML file.

This method reads a single ASCII-style pileup file from `self.input` and
converts it into an HTML representation written to `self.html`. The HTML
document includes a header with the sample group and pileup title, and a body
where each pileup line is rendered as a sequence of color-coded character boxes
followed by its associated values or annotations.

The pileup title and allowed sequence shift are inferred from the input file
name. The stylesheet referenced in the HTML output is taken from `self.css`.
Sequence characters are converted to HTML using `self.get_char_seq()`.

Input lines are interpreted according to their number of tab-delimited fields:
    - 1 field: sequence-only lines or the counts label line
    - 2 fields: sequence and value lines
    - 3 fields: sequence, value, and additional annotation lines

The final HTML document is written to the path stored in `self.html`.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L136"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `get_char_seq`

```python
get_char_seq(seq: str, shift: int) → str
```

Create HTML string for the aligned sequence.

Bases outside the canonical representation are treated as insertions and
colored as such.


**Arguments:**

 - <b>`seq`</b>: Pileup character string.
 - <b>`shift`</b>: Number of nucleotides flanking the canonical sequence.



**Returns:**

 HTML string with color-coded characters.

---

<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L116"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>function</kbd> `get_color_scheme`

```python
get_color_scheme(scheme: ColorScheme) → None
```

Create the pileup character-to-color mapping.

Populate `self.color_code` by mapping pileup characters to their background and
text hexadecimal color codes using the provided color scheme.


**Arguments:**

 - <b>`scheme`</b>:  Background color names to use for nucleotide, gap, and
   generic pileup symbols.


---

## <kbd>class</kbd> `ColorScheme`

Background color names used for pileup symbols.

**Attributes:**

 - <b>`adenine`</b>: String for the adenine background color.
 - <b>`cytosine`</b>: String for the cytosine background color.
 - <b>`guanine`</b>: String for the guanine background color.
 - <b>`thymine`</b>: String for the thymine background color.
 - <b>`gap`</b>: String for the deletions background color.
 - <b>`generic`</b>: String for the generic symbols background color.

---

## <kbd>class</kbd> `ColorValues`

Background and text colors for a pileup symbol.

 - <b>`bg`</b>: Hexadecimal string for background color.
 - <b>`text`</b>: Hexadecimal string for text color.

---

## <kbd>class</kbd> `ParsedPileupName`

Metadata derived from the pileup input file name.

- <b>`parts`</b>: List of the substrings subtracted from the pileup file name.
- <b>`title`</b>: miRNA name for which the pileup has been created.
- <b>`shift`</b>: Number of allowed bases at either end of the reference
  sequence.

---

## <kbd>class</kbd> `RowRender`

Rendered HTML fragment and current sequence length.

- <b>`html`</b>: Row HTML string representation.
- <b>`seq_len`</b>: Lenght of the sequence representation

---

