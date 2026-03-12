import pandas as pd

# Load the master panel
panel = pd.read_csv("results/panel_features.csv")

# Generate the unique feature names used by IsoDecipher
var_info = panel.copy()
var_info['feature'] = var_info.apply(lambda r: f"{r.gene}_G{r.polyA_group}_{r.user_label}", axis=1)

# Save as the master metadata file
var_info.set_index('feature').to_csv("results/feature_metadata.csv")
print("✅ Master feature metadata saved to results/feature_metadata.csv")