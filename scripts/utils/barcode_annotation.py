#!/usr/bin/env python3
"""
DataFrame-based barcode annotation for interactive analysis.

Provides high-level pandas integration for annotating cells by experimental
conditions encoded in combinatorial barcodes. Designed for use in Jupyter
notebooks and interactive analysis.

For pipeline/CLI use, see: scripts/CapWGS/assign_conditions.py
"""

import pandas as pd
from barcode_core import (
    BC_SLICES, POSITION_COLUMNS, build_barcode_lookup,
    extract_barcode_position, assign_conditions
)


def annotate_barcodes(df, barcode_table, column_name, conditions, barcode_set='A'):
    """
    Annotate a dataframe whose index contains combinatorial cell barcodes.

    Cell barcodes have the format (45 bases):
        {bcD}AGGA{bcC}ACTC{bcB}AAGG{bcA}A

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame whose index values are 45-base cell barcode strings.
    barcode_table : pd.DataFrame
        8x12 DataFrame (rows A-H, columns 1-12) of barcode sequences.
    column_name : str
        Name of the new annotation column to add.
    conditions : dict
        Maps condition labels to wells. Wells can be:
          - row letters: 'A', 'B', etc.
          - column numbers: 1, 2, etc.
          - specific wells: 'A1', 'B3', etc.
    barcode_set : str
        Which barcode position to match against: 'A', 'B', 'C', or 'D'.
        Default 'A' since conditions are typically encoded there.

    Returns
    -------
    pd.DataFrame
        Copy of df with the new column added.

    Examples
    --------
    >>> # Annotate by rows in position A
    >>> df = annotate_barcodes(
    ...     df, bcs,
    ...     column_name='amplification_method',
    ...     conditions={
    ...         'PTA': ['A', 'B', 'C', 'D'],
    ...         'mMDA': ['E', 'F', 'G', 'H'],
    ...     },
    ...     barcode_set='A'
    ... )

    >>> # Annotate by specific wells in position D
    >>> df = annotate_barcodes(
    ...     df, bcs,
    ...     column_name='replicate',
    ...     conditions={
    ...         'Rep1': ['A10', 'B10'],
    ...         'Rep2': ['A11', 'B11'],
    ...     },
    ...     barcode_set='D'
    ... )
    """
    # Build lookup using core utilities
    lookup = build_barcode_lookup(conditions, barcode_table, barcode_set, warn=False)

    # Assign conditions to all cells
    cell_barcodes = list(df.index)
    assignments = assign_conditions(cell_barcodes, lookup, barcode_set)

    # Add as new column
    result = df.copy()
    result[column_name] = [assignments[bc] for bc in result.index]

    return result


def add_qc(df, path):
    """Add QC metrics from compiled CSV to DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with cell barcodes as index
    path : str
        Path to compiled_qc_metrics.csv

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with QC columns added
    """
    qc_df = pd.read_csv(path, index_col='sample')
    qc_df = qc_df.loc[df.index]
    qc_df.rename_axis(index="sample", inplace=True)

    return pd.concat([df, qc_df], axis=1)


def ginis(df, path):
    """Load gini coefficients from individual text files.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with cell barcodes as index
    path : str
        Directory containing {barcode}_gini.txt files

    Returns
    -------
    list
        Gini coefficients in same order as df.index
    """
    ginis = []
    for name in df.index:
        with open(f"{path}/{name}_gini.txt", 'r') as fin:
            ginis.append(float(fin.readline().strip()))
    return ginis
