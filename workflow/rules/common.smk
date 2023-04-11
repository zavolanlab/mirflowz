###############################################################################
### Functions
###############################################################################


def get_sample(column_id: str, sample_id: int = None) -> str:
    """Get relevant per sample information."""
    if sample_id:
        return str(
            samples_table[column_id][samples_table.index == sample_id][0]
        )
    else:
        return str(samples_table[column_id][0])
