###############################################################################
### Functions
###############################################################################


def get_sample(column_id: str, sample_id: int = None) -> str:
    """Get relevant per sample information."""
    if sample_id:
        return str(
            samples_table[column_id][samples_table.index == sample_id].iloc[0]
        )
    else:
        return str(samples_table[column_id].iloc[0])


def convert_lib_format(lib_format: str) -> str:
    """Convert library file format."""
    formats = {
        "fa": "fa",
        "fasta": "fa",
        "FASTA": "fa",
        "fq": "fastq",
        "fastq": "fastq",
        "FASTQ": "fastq",
    }
    return formats[lib_format]
