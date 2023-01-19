import numpy as np
import pandas as pd
from typing import List

def filtering_files(dataframes: List[pd.DataFrame],
                    start_seq_length: int,
                    end_seq_length: int,
                    start_frag_len: int,
                    end_frag_len: int,
                    start_mz: float,
                    end_mz: float,
                    start_int: float,
                    end_int: float,
                    ion_filter: List[bool]) -> List[pd.DataFrame]:
    """

    """
    ion_filter_translation = ["a", "b", "c", "x", "y", "z"]
    ions_considered = []

    for i, val in enumerate(ion_filter):
        if val:
            ions_considered.append(ion_filter_translation[i])

    fragments_df = dataframes[0]
    spectra_df = dataframes[1]

    # filter by ion type
    fragments_df = fragments_df[fragments_df["frag_code"].str.contains("|".join(ions_considered))]

    fragments_df_filtered = fragments_df.query(f"frag_length > {start_frag_len} & frag_length < {end_frag_len} & frag_mz > {start_mz} & frag_mz < {end_mz} & frag_intensity > {start_int} & frag_intensity < {end_int}")
    spectra_df_filtered = spectra_df.query(f"peptide_length > {start_seq_length} & peptide_length < {end_seq_length}")

    return [fragments_df_filtered, spectra_df_filtered]
