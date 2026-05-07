import pandas as pd
import pytest
import os

PANEL_FILE = "results/panel_features_v2.csv"


@pytest.fixture(scope="module")
def panel():
    assert os.path.exists(PANEL_FILE), f"Panel file not found: {PANEL_FILE}"
    return pd.read_csv(PANEL_FILE)


def test_panel_coordinate_ordering(panel):
    """
    Verify that polyA_groups are ordered correctly by rep_coord:
    - Strand (+): rep_coord must increase across groups.
    - Strand (-): rep_coord must decrease across groups.
    """
    for gene, group in panel.groupby('gene'):
        if len(group) > 1:
            coords = group.sort_values('polyA_group')['rep_coord'].tolist()
            strand = group['strand'].iloc[0]

            if strand == '+':
                assert all(x < y for x, y in zip(coords, coords[1:])), \
                    f"Ordering error in {gene} (+): rep_coords {coords} are not strictly increasing."
            else:
                assert all(x > y for x, y in zip(coords, coords[1:])), \
                    f"Ordering error in {gene} (-): rep_coords {coords} are not strictly decreasing."


def test_coord_spread_within_tolerance(panel):
    """
    Verify that coord_spread (max - min within each group) never exceeds
    the clustering tolerance used. Default tolerance = 10bp.

    If you ran build_panel_features.py with a different --tolerance value,
    update TOLERANCE below to match.
    """
    TOLERANCE = 10

    if 'coord_spread' not in panel.columns:
        pytest.skip("coord_spread column not found — re-run build_panel_features.py to generate it.")

    violations = panel[panel['coord_spread'] > TOLERANCE]
    assert violations.empty, (
        f"{len(violations)} groups exceed tolerance of {TOLERANCE}bp:\n"
        f"{violations[['gene', 'polyA_group', 'rep_coord', 'coord_min', 'coord_max', 'coord_spread']].to_string()}"
    )


def test_no_overlapping_groups(panel):
    """
    Verify that coord ranges (coord_min to coord_max) of different groups
    within the same gene do not overlap.
    """
    if 'coord_min' not in panel.columns or 'coord_max' not in panel.columns:
        pytest.skip("coord_min/coord_max columns not found — re-run build_panel_features.py.")

    for gene, group in panel.groupby('gene'):
        if len(group) < 2:
            continue
        sorted_groups = group.sort_values('polyA_group')
        ranges = list(zip(sorted_groups['coord_min'], sorted_groups['coord_max']))

        for i in range(len(ranges) - 1):
            lo_max = max(ranges[i])
            hi_min = min(ranges[i + 1])
            assert lo_max < hi_min, (
                f"Overlapping groups in {gene}: "
                f"group {i} ends at {lo_max}, group {i+1} starts at {hi_min}"
            )


def test_immunoglobulin_labels(panel):
    """
    Verify that IGHM Group 0 is labelled 'Secreted' (proximal site).
    """
    ighm = panel[panel['gene'] == 'IGHM']

    if ighm.empty:
        pytest.skip("IGHM not found in panel — check gene list.")

    label = ighm[ighm['polyA_group'] == 0]['user_label'].iloc[0]
    assert label == "Secreted", \
        f"IGHM Group 0 should be 'Secreted', found '{label}'"


def test_ighm_has_two_groups(panel):
    """
    IGHM should always produce exactly 2 groups (Secreted + Membrane).
    If tolerance is too wide, the two sites may incorrectly merge into 1.
    """
    ighm = panel[panel['gene'] == 'IGHM']

    if ighm.empty:
        pytest.skip("IGHM not found in panel — check gene list.")

    n_groups = len(ighm)
    assert n_groups == 2, \
        f"IGHM should have 2 polyA groups, found {n_groups}. " \
        f"Check --tolerance setting (too wide merges Secreted+Membrane)."


def test_no_missing_required_columns(panel):
    """
    Verify all expected columns are present in the panel CSV.
    """
    required = {
        'gene', 'polyA_group', 'rep_coord', 'strand', 'chrom',
        'avg_spliced_utr', 'avg_genomic_utr', 'num_transcirpts',
        'transcript_ids', 'transcript_names', 'user_label',
        'coord_min', 'coord_max', 'coord_spread'
    }
    missing = required - set(panel.columns)
    assert not missing, f"Missing columns in panel: {missing}"


def test_no_null_rep_coords(panel):
    """
    rep_coord must never be null — it is the anchor for read assignment.
    """
    nulls = panel[panel['rep_coord'].isnull()]
    assert nulls.empty, \
        f"{len(nulls)} rows have null rep_coord:\n{nulls[['gene', 'polyA_group']].to_string()}"
