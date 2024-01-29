from pathlib import Path

###############################################################################
### Functions
###############################################################################


def get_sample(column_id: str, sample_id: int = None) -> str:
    """Get relevant per sample information."""
    if sample_id:
        return str(samples_table[column_id][samples_table.index == sample_id].iloc[0])
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


def create_empty_bed_file(config_file: dict, dir: Path) -> Path:
    """Create empty bed file and store its path in the configuration file."""
    empty_bed = f"{dir}/empty_file.bed"
    config_file["bed_file"] = empty_bed

    return empty_bed
