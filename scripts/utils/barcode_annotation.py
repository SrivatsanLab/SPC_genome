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
    """
    # Define barcode positions within the 45-base string
    # {bcD}AGGA{bcC}ACTC{bcB}AAGG{bcA}A
    # 0-7  8-11 12-19 20-23 24-31 32-35 36-43 44
    bc_slices = {
        'D': slice(0, 8),
        'C': slice(12, 20),
        'B': slice(24, 32),
        'A': slice(36, 44),
    }
    
    bc_slice = bc_slices[barcode_set]
    
    # Map barcode_set letter to column range in the table
    set_to_cols = {
        'A': [1, 2, 3],
        'B': [4, 5, 6],
        'C': [7, 8, 9],
        'D': [10, 11, 12],
    }
    valid_cols = set_to_cols[barcode_set]
    
    # Build lookup: barcode_seq -> condition
    barcode_to_condition = {}
    
    for condition, wells in conditions.items():
        seqs = set()
        for w in wells:
            if isinstance(w, int):
                if w not in valid_cols:
                    raise ValueError(
                        f"Column {w} is not in barcode set {barcode_set} "
                        f"(valid columns: {valid_cols})"
                    )
                seqs.update(barcode_table[w].values)
            elif isinstance(w, str) and len(w) == 1:
                for c in valid_cols:
                    seqs.add(barcode_table.loc[w, c])
            elif isinstance(w, str):
                row = w[0]
                col = int(w[1:])
                if col not in valid_cols:
                    raise ValueError(
                        f"Well {w} column {col} is not in barcode set {barcode_set} "
                        f"(valid columns: {valid_cols})"
                    )
                seqs.add(barcode_table.loc[row, col])
            else:
                raise ValueError(f"Unrecognized well spec: {w}")
        
        for seq in seqs:
            if seq in barcode_to_condition:
                raise ValueError(
                    f"Barcode {seq} assigned to both "
                    f"'{barcode_to_condition[seq]}' and '{condition}'"
                )
            barcode_to_condition[seq] = condition
    
    # Extract the relevant barcode from each cell and look up condition
    result = df.copy()
    result[column_name] = [
        barcode_to_condition.get(cell[bc_slice])
        for cell in result.index
    ]
    return result

def add_qc(df, path):
    qc_df = pd.read_csv(path, index_col='sample')
    qc_df = qc_df.loc[df.index]
    qc_df.rename_axis(index="sample", inplace=True)
    
    return pd.concat([df, qc_df], axis=1)

def ginis(df, path):
    ginis = []
    for name in df.index:
        with open(f"{path}/{name}_gini.txt", 'r') as fin:
            ginis.append(float(fin.readline().strip()))
    return ginis