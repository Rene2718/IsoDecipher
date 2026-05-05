import pandas as pd
import pytest
import os

def test_panel_coordinate_ordering():
    """
    Verify that polyA_groups are ordered correctly:
    - Strand (+): Coordinates must increase.
    - Strand (-): Coordinates must decrease.
    """
    panel_file = "results/panel_features_v2.csv"
    assert os.path.exists(panel_file), "Panel file not found!"
    
    df = pd.read_csv(panel_file)
    
    for gene, group in df.groupby('gene'):
        if len(group) > 1:
            # Sort by group index to check physical coordinate order
            coords = group.sort_values('polyA_group')['rep_coord'].tolist()
            strand = group['strand'].iloc[0]
            
            if strand == '+':
                # Check if strictly increasing
                assert all(x < y for x, y in zip(coords, coords[1:])), \
                    f"Logic Error in {gene} (+): Coords {coords} are not increasing."
            else:
                # Check if strictly decreasing
                assert all(x > y for x, y in zip(coords, coords[1:])), \
                    f"Logic Error in {gene} (-): Coords {coords} are not decreasing."

def test_immunoglobulin_labels():
    """
    Verify that IGHM labels are correctly assigned for clinical relevance.
    """
    df = pd.read_csv("results/panel_features_v2.csv")
    ighm = df[df['gene'] == 'IGHM']
    
    if not ighm.empty:
        # Group 0 should be the proximal 'Secreted' form
        label = ighm[ighm['polyA_group'] == 0]['user_label'].iloc[0]
        assert label == "Secreted", f"IGHM Group 0 should be Secreted, found {label}"