<a href="https://github.com/zavolanlab/mirflowz/blob/dev/workflow/scripts/copper.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `copper.py`

COlor-coded PileuP gEneratoR.


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

---

## <kbd>class</kbd> `ColorValues`

Background and text colors for a pileup symbol.

---

## <kbd>class</kbd> `ParsedPileupName`

Metadata derived from the pileup input file name.


---

## <kbd>class</kbd> `RowRender`

Rendered HTML fragment and current sequence length.

---

