# Output Files

In this page you can find an overview of the output files.

## Expected output files

Upon successful execution of _MIRFLOWZ_, the tool automatically removes all
intermediate files generated during the process.

The final outputs comprise:

1. A SAM file containing alignments intersecting a pri-miR locus. These
   alignments intersect with (extended) start and positions specified in the
   provided pri-miR annotations. Please note that they may not contribute to
   the final counting and may not appear in the final table.

2. A SAM file containing alignments intersecting a mature miRNA locus. Similar
   to the previous file, these alignments intersect with (extended) start and
   end positions specified in the provided miRNA annotations. Please note that
   they may not contribute to the final counting and may not appear in the
   final table.

3. A BAM file containing the set of alignments contributing to the final
   counting and its corresponding index file (`.bam.bai`).

4. Table(s) containing the counting data from all libraries for (iso)miRs
   and/or pri-miRs. Each row corresponds to a miRNA species, and each column
   represents a sample library. Each read is counted towards all the annotated
   miRNA species it aligns to, with `1/n`, where `n` is the number of genomic
   and/or transcriptomic loci that read aligns to.

5. **OPTIONAL**. ASCII-style pileups of read alignments and publication-style
   color-coded pileups produced for individual libraries, combinations of
   libraries and/or all libraries of a given run. The exact number and nature
   of the outputs depends on the workflow inputs/parameters. See the
   [pileups section](../workflow/modules/pileups.md) for a detailed
   description.

??? tip "I want to keep the intermediate files"

    To retain all intermediate files, include `--no-hooks` in the workflow
    call.

    ```sh
    snakemake \
    --snakefile="path/to/Snakefile" \
    --cores 4  \
    --configfile="path/to/config.yaml" \
    --software-deployment-method conda \
    --printshellcmds \
    --rerun-incomplete \
    --no-hooks \
    --verbose
    ```

    After a successful execution of the workflow, the intermediate files are
    going to be found in the `results/intermediates` directory.

The following file structure is expected:

```
test/results/
├── pileups
│   ├── color_coded
│   │   ├── check_file.txt
│   │   ├── all
│   │   │   ├── pileup_style.css
│   │   │   ├── all-samples_hsa-mir-A_6-shift.html
│   │   │   └── all-samples_hsa-mir-B_6-shift.html
│   │   ├── group_A
│   │   │   ├── pileup_style.css
│   │   │   ├── group_A_hsa-mir-A_6-shift.html
│   │   │   └── group_A_hsa-mir-B_6-shift.html
│   │   └── test_lib
│   │       ├── pileup_style.css
│   │       ├── test_lib_hsa-mir-A_6-shift.html
│   │       └── test_lib_hsa-mir-B_6-shift.html
│   ├── mod
│   │   ├── all
│   │   │   ├── check_file.txt
│   │   │   ├── all-samples.hsa-mir-A.6-shift.tab
│   │   │   └── all-samples.hsa-mir-B.6-shift.tab
│   │   ├── group_A
│   │   │   ├── check_file_group_A.txt
│   │   │   ├── group-A.hsa-mir-A.6-shift.tab
│   │   │   └── group-A.hsa-mir-B.6-shift.tab
│   │   └── test_lib
│   │       ├── check_file.txt
│   │       ├── lib-name.hsa-mir-A.6-shift.tab
│   │       └── lib-name.hsa-mir-B.min.6-shift.tab
│   └── raw
│       ├── all
│       │   ├── check_file.txt
│       │   ├── all_samples.hsa-mir-A.min.1.pileup.tab
│       │   └── all_samples.hsa-mir-B.min.1.pileup.tab
│       ├── group_A
│       │   ├── check_file_group_A.txt
│       │   ├── group_A.hsa-mir-A.min.1.pileup.tab
│       │   └── group_A.hsa-mir-B.min.1.pileup.tab
│       └── test_lib
│           ├── check_file.txt
│           ├── lib_name.hsa-mir-A.min.1.pileup.tab
│           └── lib_name.hsa-mir-B.min.1.pileup.tab
├── TABLES
│   ├── all_mirna_counts.tab
│   └── all_pri-mir_counts.tab
└── lib_name
    ├── alignments_intersecting_mirna.sam
    ├── alignments_intersecting_mirna_uncollapsed_sorted.bam
    ├── alignments_intersecting_mirna_uncollapsed_sorted.bam.bai
    └── alignments_intersecting_primir.sam
```

