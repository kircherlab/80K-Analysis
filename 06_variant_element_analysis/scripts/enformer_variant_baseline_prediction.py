import pandas as pd
import numpy as np
import yaml
# for the correlation plotting
import seaborn as sns
import matplotlib.pyplot as plt
# for the lasso regression
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler
# -----

# load helpful functions
import sys
sys.path.append('../../00_helpful_functions')
import helpful_functions as hf

config_path = "../../global80K_config.yaml"
# load config file
with open(config_path, "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

run_name = 'all_pred'
# load variant table and enformer prediction table
variant_table = pd.read_csv(config['files']['final_design']['variant_table'], sep="\t")
variant_map = variant_table.drop(columns=['Region']).copy()
variant_map.head()



# add variant effect from bcalm:
bcalm_bbmap_35_df = pd.read_csv(config['files']['creating']['bc_MPRAlm'], sep="\t")
bcalm_bbmap_35_df.head() # 38,364
variant_map_effect = variant_map.merge(bcalm_bbmap_35_df[['variant_id', 'logFC']], left_on='ID', right_on='variant_id', how='inner')
variant_map_effect.drop(columns=['variant_id'], inplace=True)
# print(variant_map_effect.shape[0]) # 38205


# currently focus on 680/5k columns
# enformer_sequence_predictions_path = '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/projects/enformer_predictions/element_enformer_prediction/sequences_for_enformer_predictions_only_sequences_with_padding_680.tsv.gz'
# enformer_sequence_predictions_path = '/data/cephfs-2/unmirrored/groups/kircher/MPRA/IGVF_Y1_design/projects/enformer_predictions/element_enformer_prediction/sequences_for_enformer_predictions_only_sequences_with_padding_all_columns.tsv.gz'
enformer_sequence_predictions_path = '/home/kisa/coding/80K_MPRA/enformer_data/baseline_activity/sequences_for_enformer_predictions_only_sequences_with_padding_all_columns.tsv.gz'
all_enformer_sequence_predictions_df = pd.read_csv(enformer_sequence_predictions_path, sep="\t") # 80215
all_enformer_sequence_predictions_df.columns.values[0] = 'sequence_ID'

# Merge for REF
ref_merged = variant_map_effect.merge(all_enformer_sequence_predictions_df, left_on='REF', right_on='sequence_ID', how='inner')
ref_merged.drop(columns=['sequence_ID'], inplace=True)
ref_merged = ref_merged.rename(columns=lambda x: f"REF_{x}" if x not in variant_map_effect.columns else x)
ref_merged
# # # Merge for ALT
alt_merged = ref_merged.merge(all_enformer_sequence_predictions_df, left_on='ALT', right_on='sequence_ID', how='inner')
alt_merged.drop(columns=['sequence_ID'], inplace=True)
alt_merged = alt_merged.rename(columns=lambda x: f"ALT_{x}" if x not in ref_merged.columns else x)
# alt_merged = alt_merged.drop(columns=['ALT_ID'])


# Compute absolute differences
ref_cols = [col for col in alt_merged.columns if col.startswith('REF_')]
alt_cols = [col for col in alt_merged.columns if col.startswith('ALT_')]
diffs = pd.DataFrame()

for ref, alt in zip(ref_cols, alt_cols):
    diff_col = f"diff_{ref.split('_', 1)[1]}"
    diffs[diff_col] = alt_merged[alt] - alt_merged[ref]

# Add differences back to the table
alt_merged = pd.concat([alt_merged, diffs], axis=1)

# Focus on diff columns
alt_merged = alt_merged.loc[:, ~alt_merged.columns.str.startswith(('REF_', 'ALT_'))]
# alt_merged # 38180


# Features (diff_* columns) and target
X = alt_merged.filter(like="diff_")
y = alt_merged['logFC']

# Standardize features
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# Perform Lasso regression with cross-validation
lasso = LassoCV(cv=5, random_state=42, max_iter=10000).fit(X_scaled, y)

# Extract feature importance
feature_importance = pd.Series(lasso.coef_, index=X.columns)
important_features = feature_importance[feature_importance != 0]

# Display feature importance
print("Feature Importance:")
print(important_features)

# Plot feature importance
plt.figure(figsize=(10, 6))
important_features.sort_values(ascending=False).plot(kind='bar')
plt.title("Feature Importance from Lasso Regression")
plt.ylabel("Coefficient Value")
plt.xlabel("Features")

# Evaluate model performance
r2_score = lasso.score(X_scaled, y)
print(f"R^2 Score: {r2_score:.3f}")

# Paths for saving the files
important_features_path = f"/home/kisa/coding/80K_MPRA/modeling/enformer_variant_baseline_sequence_predictions/{run_name}_important_features.csv"  # Adjust path as needed
r2_score_path = f"/home/kisa/coding/80K_MPRA/modeling/enformer_variant_baseline_sequence_predictions/{run_name}_r2_score.txt"
plot_path = f"/home/kisa/coding/80K_MPRA/modeling/enformer_variant_baseline_sequence_predictions/{run_name}_lasso_feature_importance.png"

# Save feature importance to CSV
important_features.to_csv(important_features_path, header=["Coefficient"])

# Save R^2 score to a text file
with open(r2_score_path, "w") as file:
    file.write(f"R^2 Score: {r2_score:.3f}")

# Create and save the feature importance plot
plt.figure(figsize=(10, 6))
important_features.sort_values(ascending=False).plot(
    kind='bar',
    color=(important_features > 0).map({True: 'green', False: 'red'})
)
plt.title("Feature Importance from Lasso Regression (Positive and Negative)")
plt.ylabel("Coefficient Value")
plt.xlabel("Features")
plt.axhline(0, color='black', linestyle='--', linewidth=0.8)  # Horizontal line at 0
plt.tight_layout()
plt.savefig(plot_path)
plt.close()

print(f"Feature importance saved to: {important_features_path}")
print(f"R^2 score saved to: {r2_score_path}")
print(f"Feature importance plot saved to: {plot_path}")