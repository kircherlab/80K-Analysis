## Load the regions bed files and perform the enrichment analysis vs remap2022

library(dplyr)
# install.packages("tidyr")
library("tidyr")
library(ggplot2)
library(regioneReloaded)
# library(BSgenome.Hsapiens.UCSC.hg38)
input_rds_path = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/sampling_1000_ngn2_with_background_full_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"
undiffWTC11_input_rds_path = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/sampling_1000_undiffWTC11_with_background_full_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"
ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix <- readRDS(input_rds_path)
undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix <- readRDS(undiffWTC11_input_rds_path)


output_directory = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/"
ngn2_and_undiffWTC11_input_rds_path = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/ngn2_undiffWTC11_sampling_1000_with_background_full_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"
ngn2_and_undiffWTC11_increasing_decreasing_CrosswiseMatrix <- readRDS(ngn2_and_undiffWTC11_input_rds_path)



# 1. Investigate the object and understand the structure
# get the names
ngn2_and_undiffWTC11_increasing_decreasing_CrosswiseMatrix@multiOverlaps %>%
    names()

ngn2_and_undiffWTC11_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_upregulating_dCREs %>%
    names()

# get the type
ngn2_and_undiffWTC11_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_upregulating_dCREs %>%
    class()

# no significant results at all: (no background used)
ngn2_and_undiffWTC11_input_rds_path = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/ngn2_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"

# --- Load Data ---
output_directory <- "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/"
ngn2_and_undiffWTC11_input_rds_path = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/ngn2_dCREs_with_and_without_epigenetic_marks_sampling_1000_with_background_full_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"

ngn2_and_undiffWTC11_input_rds_path = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/ngn2_undiffWTC11_sampling_1000_with_background_full_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"
ngn2_and_undiffWTC11_input_rds_path = "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/sampling_1000_ngn2_with_background_full_remap_TF_chip_increasing_decreasing_CrosswiseMatrix.rds"
ngn2_and_undiffWTC11_CrosswiseMatrix <- readRDS(ngn2_and_undiffWTC11_input_rds_path)

# Significance level threshold
sig_level <- 0.1


# --- Step 1-3: Iterate through multiOverlaps and Identify Significant Entries ---

# --- Extract all possible TFBS names (from the original matrix column names) ---
# This ensures our final table includes ALL TFBS, not just those that were significant.
matrix_data <- regioneReloaded::getMatrix(ngn2_and_undiffWTC11_CrosswiseMatrix)
all_tfbs_names <- tibble(Name = colnames(matrix_data)) # Create a base table of all TFBS names

# --- Step 1-3: Iterate through multiOverlaps and Identify Significant Entries ---

# Get the names of the region sets in multiOverlaps (Alist names)
alist_names <- names(ngn2_and_undiffWTC11_CrosswiseMatrix@multiOverlaps)

# Initialize an empty data frame to collect all significant entries for pivot table.
all_significant_hits <- tibble() # Using tibble() creates a 0-row, 0-column tibble


# Initialize a list to store significant results for each Alist name
significant_results_by_group_list <- list()


for (alist_name in alist_names) {
    # Get the data frame for the current Alist region set
    current_df <- ngn2_and_undiffWTC11_CrosswiseMatrix@multiOverlaps[[alist_name]]

    # Filter for significant entries based on adj.p_value AND THEN add the Alist_Region_Set column
    significant_df_with_context <- current_df %>%
        filter(adj.p_value < sig_level) %>%
        mutate(
            Alist_Region_Set = alist_name,
            Enrichment_Status = case_when(
                `norm_zscore` > 0 ~ "Enriched", # Use `z-score` or `normZ` as per your data
                `norm_zscore` < 0 ~ "Depleted",
                TRUE ~ "Neutral" # Should not happen with significant p-values, but good for completeness
            )
        ) # This safely adds the column, even if 0 rows

    # Store the significant TFBS (names from Blist) for current Alist region set in our new list
    if (nrow(significant_df_with_context) > 0) {
        significant_results_by_group_list[[alist_name]] <- significant_df_with_context$name
    } else {
        significant_results_by_group_list[[alist_name]] <- character(0) # Store empty character vector if no hits
    }

    # Bind to the collection data frame for the overall pivot table and combined TFBS list
    all_significant_hits <- bind_rows(all_significant_hits, significant_df_with_context)
}

# --- Step 4: Generate Pivot Table of Significant Results ---

# Get the counts from the significant hits
counts_df <- all_significant_hits %>%
    group_by(Alist_Region_Set) %>%
    summarise(
        NumberOfSignificantTFBS = n(),
        NumberOfEnrichedTFBS = sum(Enrichment_Status == "Enriched"),
        NumberOfDepletedTFBS = sum(Enrichment_Status == "Depleted")
    ) %>%
    ungroup()

# Create a data frame of all original Alist names
all_alist_names_df <- tibble(Alist_Region_Set = alist_names)

# Join with the counts, filling NA with 0 for region sets with no significant hits
pivot_table_significant_counts <- all_alist_names_df %>%
    left_join(counts_df, by = "Alist_Region_Set") %>%
    replace_na(list(
        NumberOfSignificantTFBS = 0,
        NumberOfEnrichedTFBS = 0,
        NumberOfDepletedTFBS = 0
        )) %>%
    # Arrange for better readability
    arrange(desc(NumberOfSignificantTFBS), Alist_Region_Set)


# Print the pivot table
print("Pivot Table: Number of Significant TFBS per Alist Region Set")
print(pivot_table_significant_counts)

# --- Store and Print the list of all unique significantly enriched TFBS (overall) ---
list_of_all_unique_significant_TFBS <- unique(all_significant_hits$name)

# --- Store and Print the list of significant TFBS for EACH group ---
print("\nList of significant TFBS for EACH Alist Region Set (grouped):")
if (length(names(significant_results_by_group_list)) > 0) {
    print(significant_results_by_group_list)
} else {
    print("No significant TFBS found for any group at the specified significance level.")
}


############################ CODE TO get the depleted enriched information for the motifs:

# 1. Prepare data for pivoting: Select relevant columns and add 'value' based on enrichment status
prepared_for_pivot_enriched <- all_significant_hits %>%
    filter(Enrichment_Status == "Enriched") %>% # Filter for enriched
    select(Name = name, Alist_Region_Set) %>%
    mutate(is_enriched = 1) # Assign 1 to indicate enriched presence

prepared_for_pivot_depleted <- all_significant_hits %>%
    filter(Enrichment_Status == "Depleted") %>% # Filter for depleted
    select(Name = name, Alist_Region_Set) %>%
    mutate(is_depleted = 1) # Assign 1 to indicate depleted presence

# 2. Pivot wider for Enriched: Transform rows into columns
combined_enriched_table <- prepared_for_pivot_enriched %>%
    pivot_wider(
        id_cols = Name, # TFBS names will be the rows
        names_from = Alist_Region_Set, # Alist region sets will become new columns
        values_from = is_enriched, # Fill with the 'is_enriched' value (1)
        values_fill = 0, # Fill non-significant cells with 0
        names_prefix = "enriched_in_" # Add a prefix to the new column names
    )

# 3. Pivot wider for Depleted: Transform rows into columns
combined_depleted_table <- prepared_for_pivot_depleted %>%
    pivot_wider(
        id_cols = Name, # TFBS names will be the rows
        names_from = Alist_Region_Set, # Alist region sets will become new columns
        values_from = is_depleted, # Fill with the 'is_depleted' value (1)
        values_fill = 0, # Fill non-significant cells with 0
        names_prefix = "depleted_in_" # Add a prefix to the new column names
    )

# 4. Join the enriched and depleted tables back together with all original TFBS names
final_combined_enrichment_depletion_table <- all_tfbs_names %>%
    left_join(combined_enriched_table, by = "Name") %>%
    left_join(combined_depleted_table, by = "Name") %>%
    # Replace NAs with 0 for columns that didn't have a match in either enriched or depleted sets
    replace_na(list(
        !!!setNames(rep(0, length(alist_names)), paste0("enriched_in_", alist_names)),
        !!!setNames(rep(0, length(alist_names)), paste0("depleted_in_", alist_names))
    )) %>%
    arrange(Name) # Arrange by TFBS name for readability

# Print the final combined table
print("\nCombined Table: TFBS Enrichment/Depletion Status (1 = present, 0 = absent)")
print(final_combined_enrichment_depletion_table)



write.csv(final_combined_enrichment_depletion_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_and_undiffWTC11_enriched_depleted_remap2022_combined_tfbs_presence.csv", row.names = FALSE)
# write.csv(final_combined_enrichment_depletion_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_with_and_without_epigenetic_marks_enriched_depleted_combined_tfbs_presence.csv", row.names = FALSE)
# write.csv(final_combined_enrichment_depletion_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_high_low_enriched_depleted_remap2022_combined_tfbs_presence.csv", row.names = FALSE)


# write.csv(final_combined_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_with_and_without_epigenetic_marks_combined_tfbs_presence.csv", row.names = FALSE)
# write.csv(final_combined_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_high_low_remap2022_combined_tfbs_presence.csv", row.names = FALSE)
# write.csv(final_combined_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_and_undiffWTC11_remap2022_combined_tfbs_presence.csv", row.names = FALSE)


######################## Code to get only the Absent / Present info for the TFBS ChIP in the sequence sets
# --- Generate the Combined Table of TFBS Presence/Absence ---

# 1. Prepare data for pivoting: Select relevant columns and add a 'value' indicating presence
prepared_for_pivot <- all_significant_hits %>%
    select(Name = name, Alist_Region_Set) %>%
    mutate(is_significant = 1) # Assign 1 to indicate presence

# 2. Pivot wider: Transform rows into columns
combined_presence_table <- prepared_for_pivot %>%
    pivot_wider(
        id_cols = Name, # TFBS names will be the rows
        names_from = Alist_Region_Set, # Alist region sets will become new columns
        values_from = is_significant, # Fill with the 'is_significant' value (1)
        values_fill = 0, # Fill non-significant cells with 0
        names_prefix = "is_in_" # Add a prefix to the new column names
    )

# 3. Ensure all original TFBS names are present:
# Perform a left_join with the master list of all TFBS names.
# This will add TFBS that were never significant in any group, filling their new columns with 0.
final_combined_table <- all_tfbs_names %>%
    left_join(combined_presence_table, by = "Name") %>%
    # If a TFBS was never significant, the columns added by left_join will be NA.
    # Replace those NAs with 0.
    replace_na(list(
        # Dynamically create the list of columns to replace NA for
        !!!setNames(rep(0, length(alist_names)), paste0("is_in_", alist_names))
    )) %>%
    arrange(Name) # Arrange by TFBS name for readability

# Print the final combined table
print("\nCombined Table: TFBS Presence in Significant Sets (1 = present, 0 = absent)")
print(final_combined_table)


# write.csv(final_combined_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_with_and_without_epigenetic_marks_combined_tfbs_presence.csv", row.names = FALSE)
# write.csv(final_combined_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_high_low_remap2022_combined_tfbs_presence.csv", row.names = FALSE)
# write.csv(final_combined_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_and_undiffWTC11_remap2022_combined_tfbs_presence.csv", row.names = FALSE)

# TFBS_of_Neuron_upregulating_dCREs = ["ASCL1","ASXL3","ATF1","ATF2","ATF3","BAP1","BHLHE40","BMI1","BRD3","CHD2","CHD7","CREB1","CTCF","CTNNB1","DPF1","DRAP1","E2F6","EBF3","EP300","EZH2","FLI1","FOSL2","FOXK2","FOXP2","GFI1B","GRHL2","GTF2F1","HAND2","HDGF","ISL1","JUND","KAT2B","KDM3A","KLF10","KLF13","KLF9","LIN54","MAFB","MAX","MAZ","MITF","MNT","MRTFB","MTA2","MXI1","MYC","MYCN","NELFA","NEUROD1","NEUROG2","NFE2","NFYA","NFYB","NFYC","NKX3-1","NR2C1","OLIG2","ONECUT1","ONECUT2","PAX7","PDX1","PHF21A","PHF8","PITX3","PKNOX1","PLAG1","POU2F1","POU2F3","POU4F2","POU5F1" ,"RAD51","RARA","RBBP4","RBBP5","RBM22","RBM39","RCOR1","REST","RFX1","RFX5","SETDB1","SIRT6","SMC3","SOX13","SOX2","SP2","SP4","SP5","SP7","T","TAF1","TBL1XR1","TBP","TBX2","TP53","TRIM28","UBTF","USF1","USF2","YY1","YY2","ZBTB42","ZFP42","ZNF263","ZNF341" ,"ZNF350","ZNF423","ZNF449","ZNF467","ZNF528","ZNF554"]

# TFBS_of_Neuron_downregulating_dCREs = ["ASCL1","ATF1","BMI1","BRD9","CHD2","CHD7","DDX20","DNMT3B","DPF1","E2F4","E4F1","EOMES","EP300","FOSL2","FOXA2","FOXK2","FOXP2","GATA6","HAND2","HDGF","HNF1B"  ,"HOXB8","ID3","JUND","KLF9","MAX","MAZ","MITF","MNT","NANOG","NEUROD1","NFYA","NFYB","NFYC","ONECUT2","PLAG1","POU5F1","RAD21","RAD51","REST","RFX1","RXRA","SMAD2","SMARCA2","SMARCA4","SMARCC1","SOX2","T","TBL1XR1","TBX2","TCF12","TERF2","TP53","TRIM28","USF1","USF2","ZNF207","ZNF263","ZNF281","ZNF462","ZNF467","ZNF592","ZZZ3"]

# TFBS_of_undiffWTC11_upregulating_dCREs = ["ATF2","BCL11A","BCOR","BRD9","CHD2","CHD7","CREB1","CTBP2","CTNNB1","DDX20","E2F6","EGR1","EOMES","EP300","ETS1","FOS","FOXP2","GABPA","GATA6","GFI1B","GLIS3","GTF3C2","HDAC8","HDGF","JUN","JUND","KDM4B","KLF4","LEF1","MAX","MORC2","NANOG","NIPBL","PCGF1","POU3F1","POU5F1","PRDM14","RAD51","RXRA","RYBP","SALL3","SIN3A","SMAD2","SMAD2-3","SMAD3","SMARCA4","SMARCB1","SMARCC1","SOX2","SOX3","SP1","SP3","SP4","T","TAF1","TAL1","TBP","TCF12","TCF3","TEAD4","TFAP2A","TP53","TRIM28","VEZF1","XRCC5","YAP1","YY1AP1","ZNF114","ZNF148","ZNF207","ZNF398","ZNF462","ZNF639"]


# Check the overlap between the significant TFBS results
# We'll focus on "Neuron_upregulating_dCREs" and "undiffWTC11_upregulating_dCREs"
neuron_up_tfbs <- significant_results_by_group_list[["Neuron_upregulating_dCREs"]]
undiffWTC11_up_tfbs <- significant_results_by_group_list[["undiffWTC11_upregulating_dCREs"]]
neuron_down_tfbs <- significant_results_by_group_list[["Neuron_downregulating_dCREs"]]

# Prepare data for UpSetR plot
list_for_upset <- list(
    Neuron_up = neuron_up_tfbs,
    undiffWTC11_up = undiffWTC11_up_tfbs
)

# Remove empty sets if any, to avoid issues with UpSetR
list_for_upset <- list_for_upset[sapply(list_for_upset, length) > 0]

if (length(list_for_upset) > 0) {
    print("Generating UpSet plot for 'upregulating' TFBS sets:")
    tryCatch(
        {
            # Plotting using UpSetR
            # If sets are too small, UpSetR might show fewer intersections.
            # min.size and max.size can be adjusted, or just let default.
            upset(fromList(list_for_upset),
                order.by = "freq",
                mainbar.y.label = "Intersection Size",
                sets.y.label = "Set Size"
            )
        },
        error = function(e) {
            message("Could not generate UpSet plot, possibly due to insufficient data for intersections or other error: ", e$message)
            message("Counts for plotting:")
            print(sapply(list_for_upset, length))
        }
    )
} else {
    print("No data available to generate UpSet plot for 'upregulating' TFBS sets (empty sets).")
}

# # Calculate intersections and unique sets for TFBS selection
common_up_tfbs <- intersect(neuron_up_tfbs, undiffWTC11_up_tfbs)
unique_neuron_up_tfbs <- setdiff(neuron_up_tfbs, undiffWTC11_up_tfbs)
unique_undiffWTC11_up_tfbs <- setdiff(undiffWTC11_up_tfbs, neuron_up_tfbs)

# Print counts for verification
cat("\n--- TFBS Set Sizes for Selection ---:\n")
cat("Neuron_upregulating_dCREs (total significant):", length(neuron_up_tfbs), "\n")
cat("undiffWTC11_upregulating_dCREs (total significant):", length(undiffWTC11_up_tfbs), "\n")
cat("TFBS common to both 'upregulating' sets:", length(common_up_tfbs), "\n")
cat("TFBS unique to 'Neuron_upregulating_dCREs':", length(unique_neuron_up_tfbs), "\n")
cat("TFBS unique to 'undiffWTC11_upregulating_dCREs':", length(unique_undiffWTC11_up_tfbs), "\n")
cat("Neuron_downregulating_dCREs (total significant):", length(neuron_down_tfbs), "\n")


# --- Step 1.b: Select 5 TFBS from each category based on MOST ENRICHED (norm_zscore) ---
num_to_select <- 5
selected_tfs_for_plot <- c()

# Helper function to get top TFBS by norm_zscore for a given set and TFBS list
# This function will filter 'all_significant_hits' by the relevant Alist_Region_Set
# and the specific TFBS names for the category (e.g., unique TFBS).
get_top_enriched_tfbs <- function(tfbs_names_in_category, alist_region_set_name, n_select, all_sig_hits_df, filter_multi_cell_type_significant = FALSE) {
    if (length(tfbs_names_in_category) == 0) {
        return(character(0))
    }

    temp_all_sig_hits_df <- all_sig_hits_df # Work on a copy to avoid modifying global df

    # Step 1: Identify TFBS significant in multiple cell types, IF the filter option is TRUE
    if (filter_multi_cell_type_significant) {
        # Count how many distinct Alist_Region_Set each TFBS is significant in
        tfbs_multi_significant_counts <- all_sig_hits_df %>% # Use original all_sig_hits_df for full context
            group_by(name) %>%
            summarise(num_sets_significant = n_distinct(Alist_Region_Set)) %>%
            ungroup()

        # Identify TFBS that are significant in more than one set
        multi_set_tfbs <- tfbs_multi_significant_counts %>%
            filter(num_sets_significant > 1) %>%
            pull(name)

        # Filter out these multi-set TFBS from the working dataframe for this specific selection
        temp_all_sig_hits_df <- temp_all_sig_hits_df %>%
            filter(!(name %in% multi_set_tfbs))

        # Also, adjust the list of TFBS names for this specific category
        # to only include those that are *not* significant in multiple sets.
        tfbs_names_in_category <- tfbs_names_in_category[!(tfbs_names_in_category %in% multi_set_tfbs)]

        # If, after filtering, there are no TFBS left in the category, return empty
        if (length(tfbs_names_in_category) == 0) {
            message(paste0("  Note: All TFBS in category '", alist_region_set_name, "' were filtered out as they were significant in multiple cell types."))
            return(character(0))
        }
    }

    # Step 2: Apply original filtering, sorting, and selection
    filtered_and_sorted <- temp_all_sig_hits_df %>%
        filter(
            Alist_Region_Set == alist_region_set_name,
            name %in% tfbs_names_in_category # This list has been updated if filter_multi_cell_type_significant was TRUE
        ) %>%
        arrange(desc(norm_zscore)) %>% # Sort by norm_zscore in descending order (most enriched)
        head(n_select) %>%
        pull(name) # Extract just the 'name' column as a vector

    return(filtered_and_sorted)
}


# --- Step 1.b: Select 5 TFBS from each category based on MOST ENRICHED (norm_zscore) ---
message("\n--- Selecting Most Enriched TFBS for Plotting ---")

# 1. TFBS unique for "Neuron_upregulating_dCREs" (apply multi-cell type filter)
selected_tfs_for_plot$unique_Neuron_up <- get_top_enriched_tfbs(unique_neuron_up_tfbs, "Neuron_upregulating_dCREs", num_to_select, all_significant_hits, filter_multi_cell_type_significant = TRUE)
message(paste0("Selected ", length(selected_tfs_for_plot$unique_Neuron_up), " most enriched TFBS unique to 'Neuron_upregulating_dCREs' (excluding multi-cell type significant):"))
print(selected_tfs_for_plot$unique_Neuron_up)


# 2. TFBS unique for "undiffWTC11_upregulating_dCREs" (apply multi-cell type filter)
selected_tfs_for_plot$unique_undiffWTC11_up <- get_top_enriched_tfbs(unique_undiffWTC11_up_tfbs, "undiffWTC11_upregulating_dCREs", num_to_select, all_significant_hits, filter_multi_cell_type_significant = TRUE)
message(paste0("Selected ", length(selected_tfs_for_plot$unique_undiffWTC11_up), " most enriched TFBS unique to 'undiffWTC11_upregulating_dCREs' (excluding multi-cell type significant):"))
print(selected_tfs_for_plot$unique_undiffWTC11_up)


# 3. TFBS common to both "upregulating" sets (DO NOT apply multi-cell type filter)
if (length(common_up_tfbs) > 0) {
    selected_tfs_for_plot$common_up <- all_significant_hits %>%
        filter(
            Alist_Region_Set %in% c("Neuron_upregulating_dCREs", "undiffWTC11_upregulating_dCREs"),
            name %in% common_up_tfbs
        ) %>%
        group_by(name) %>%
        summarise(max_norm_zscore = max(norm_zscore, na.rm = TRUE)) %>%
        ungroup() %>%
        arrange(desc(max_norm_zscore)) %>%
        head(num_to_select) %>%
        pull(name)
    message(paste0("Selected ", length(selected_tfs_for_plot$common_up), " most enriched TFBS common to both upregulating sets (ranked by max norm_zscore):"))
    print(selected_tfs_for_plot$common_up)
} else {
    message("No common TFBS found for upregulating sets to select from.")
    selected_tfs_for_plot$common_up <- character(0) # Ensure it's an empty vector
}


# 4. TFBS significant for "Neuron_downregulating_dCREs" (apply multi-cell type filter for specificity)
selected_tfs_for_plot$Neuron_down <- get_top_enriched_tfbs(neuron_down_tfbs, "Neuron_downregulating_dCREs", num_to_select, all_significant_hits, filter_multi_cell_type_significant = TRUE)
message(paste0("Selected ", length(selected_tfs_for_plot$Neuron_down), " most enriched TFBS from 'Neuron_downregulating_dCREs' (excluding multi-cell type significant):"))
print(selected_tfs_for_plot$Neuron_down)


# Ensure uniqueness of the final list of TFBS for plotting
selected_tfs_for_plot <- unique(selected_tfs_for_plot)

if (length(selected_tfs_for_plot) == 0) {
    stop("No TFBS were selected for plotting after applying all filters. Please check your criteria.")
}
message(paste0("\nTotal unique TFBS selected for plotting: ", length(selected_tfs_for_plot)))
print(selected_tfs_for_plot)


# --- Step 2: Subset the matrix and set it into an existing copy of a CrosswiseMatrix ---
message("\nPreparing CrosswiseMatrix for plotting selected significant TFBS...")

selected_tfs_for_plot_filtered <- selected_tfs_for_plot[selected_tfs_for_plot %in% colnames(matrix_data)]

if (length(selected_tfs_for_plot_filtered) == 0) {
    stop("None of the specifically selected TFBS were found in the main matrix data. Cannot generate plot.")
}

subset_matrix_data <- matrix_data[, selected_tfs_for_plot_filtered]

message(paste0("Subsetted matrix dimensions: ", paste(dim(subset_matrix_data), collapse = "x")))
message(paste0("Plotting ", ncol(subset_matrix_data), " TFBS across ", nrow(subset_matrix_data), " dCRE sets."))

ngn2_and_undiffWTC11_CrosswiseMatrix_subset <- ngn2_and_undiffWTC11_CrosswiseMatrix
ngn2_and_undiffWTC11_CrosswiseMatrix_subset@matrix$GMat <- subset_matrix_data

if (!is.null(ngn2_and_undiffWTC11_CrosswiseMatrix_subset@matrix$FitCol)) {
    ngn2_and_undiffWTC11_CrosswiseMatrix_subset@matrix$FitCol <- NULL
    message("Removed old column clustering information (FitCol) to allow re-clustering on subset.")
}

# --- Step 3: Plot the resulting genoMatriXeR using plotCrosswiseMatrix() ---
message("\nGenerating plotCrosswiseMatrix for the selected significant TFBS subset...")

subset_plot <- plotCrosswiseMatrix(ngn2_and_undiffWTC11_CrosswiseMatrix_subset)

print(subset_plot)

ggsave("/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/ngn2_and_undiffWTC11_subset_TFBS_crosswise_matrix_plot_max_enrichment.pdf", plot = subset_plot, width = 10, height = 8)
# ggsave("/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/ngn2_and_undiffWTC11_subset_TFBS_crosswise_matrix_plot_not_sorted.pdf", plot = subset_plot, width = 10, height = 8)

# TODO: For the significant results: check the effect size of the enrichment / depletion

# --- NEW: Prepare data for pivoting, split by enrichment/depletion ---

# --- Modified Step 4: Generate Pivot Table of Significant Results with Enrichment/Depletion Split ---

message("\nGenerating Pivot Table of Significant Results (Enrichment/Depletion Split)...")

# Get the counts from the significant hits, splitting by direction
counts_df_with_direction <- all_significant_hits %>%
  group_by(Alist_Region_Set) %>%
  summarise(
    # Total counts (regardless of direction)
    Total_Significant_TFBS = n(),

    # Counts for enriched TFBS (norm_zscore > 0)
    Num_Enriched_TFBS = sum(norm_zscore > 0, na.rm = TRUE), # na.rm=TRUE is for safety

    # Counts for depleted TFBS (norm_zscore < 0)
    Num_Depleted_TFBS = sum(norm_zscore < 0, na.rm = TRUE), # na.rm=TRUE is for safety
  ) %>%
  ungroup() %>%
  # Calculate proportions within each set (handle division by zero for Total_Significant_TFBS)
  mutate(
    Prop_Enriched_in_Set = ifelse(Total_Significant_TFBS > 0, Num_Enriched_TFBS / Total_Significant_TFBS, 0),
    Prop_Depleted_in_Set = ifelse(Total_Significant_TFBS > 0, Num_Depleted_TFBS / Total_Significant_TFBS, 0)
  )


# Create a data frame of all original Alist names (to ensure all sets are in the table, even if no hits)
all_alist_names_df <- tibble(Alist_Region_Set = alist_names)

# Join with the new counts, filling NA with 0 for region sets with no significant hits at all
pivot_table_significant_counts_direction <- all_alist_names_df %>%
  left_join(counts_df_with_direction, by = "Alist_Region_Set") %>%
  # Dynamically replace NAs for all calculated count and proportion columns with 0
  replace_na(list(
    Total_Significant_TFBS = 0,
    Num_Enriched_TFBS = 0,
    Num_Depleted_TFBS = 0,
    Prop_Enriched_in_Set = 0,
    Prop_Depleted_in_Set = 0
  )) %>%
  # Arrange for better readability (e.g., by total significant TFBS count)
  arrange(desc(Total_Significant_TFBS), Alist_Region_Set)


# Print the new pivot table
print("\nPivot Table: Number of Significant TFBS per Alist Region Set (Enrichment/Depletion Split)")
print(pivot_table_significant_counts_direction)

write.csv(pivot_table_significant_counts_direction, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_and_undiffWTC11_remap2022_combined_tfbs_pivot_presence_enrichment_depletion.csv", row.names = FALSE)

# Filter for enriched TFBS (significant AND norm_zscore > 0)
enriched_data_for_pivot <- all_significant_hits %>%
    filter(norm_zscore > 0) %>% # Filter for positive z-score
    select(Name = name, Alist_Region_Set) %>%
    mutate(is_enriched = TRUE) # Indicate presence for enrichment

# Filter for depleted TFBS (significant AND norm_zscore < 0)
depleted_data_for_pivot <- all_significant_hits %>%
    filter(norm_zscore < 0) %>% # Filter for negative z-score
    select(Name = name, Alist_Region_Set) %>%
    mutate(is_depleted = TRUE) # Indicate presence for depletion

# --- Pivot wider for enriched data ---
pivot_enriched_table <- enriched_data_for_pivot %>%
    pivot_wider(
        id_cols = Name,
        names_from = Alist_Region_Set,
        values_from = is_enriched,
        values_fill = FALSE, # Fill non-enriched cells with FALSE
        names_prefix = "enriched_in_"
    )

# --- Pivot wider for depleted data ---
pivot_depleted_table <- depleted_data_for_pivot %>%
    pivot_wider(
        id_cols = Name,
        names_from = Alist_Region_Set,
        values_from = is_depleted,
        values_fill = FALSE, # Fill non-depleted cells with FALSE
        names_prefix = "depleted_in_"
    )

# --- Combine the enriched and depleted pivot tables ---
# Use full_join to ensure TFBS that are only enriched or only depleted are kept
combined_direction_table <- pivot_enriched_table %>%
    full_join(pivot_depleted_table, by = "Name")

# --- Ensure all original TFBS names are present and fill NAs ---
# Perform a left_join with the master list of all TFBS names.
# This will add TFBS that were never significant in ANY direction.
final_combined_table <- all_tfbs_names %>%
    left_join(combined_direction_table, by = "Name")

# Dynamically generate the list of all expected 'enriched_in_' and 'depleted_in_' columns
all_expected_direction_cols <- c(
    paste0("enriched_in_", alist_names),
    paste0("depleted_in_", alist_names)
)

# Replace any remaining NAs (for TFBS never significant) with FALSE
# We filter to only those columns that actually exist in the table after the join
existing_direction_cols <- all_expected_direction_cols[all_expected_direction_cols %in% colnames(final_combined_table)]

# Create a named list for replace_na, where names are columns and values are FALSE
na_replace_list <- as.list(setNames(rep(FALSE, length(existing_direction_cols)), existing_direction_cols))

final_combined_table <- final_combined_table %>%
    replace_na(na_replace_list) %>%
    arrange(Name) # Arrange by TFBS name for readability

# --- Calculate Sums and Proportions ---
# Identify the enriched and depleted columns for calculation
enriched_specific_cols <- grep("^enriched_in_", colnames(final_combined_table), value = TRUE)
depleted_specific_cols <- grep("^depleted_in_", colnames(final_combined_table), value = TRUE)

final_combined_table <- final_combined_table %>%
    rowwise() %>% # Perform calculations row by row for each TFBS
    mutate(
        Total_Enriched_Occurrences = sum(c_across(all_of(enriched_specific_cols))),
        Total_Depleted_Occurrences = sum(c_across(all_of(depleted_specific_cols))),
        Total_Significant_Occurrences = Total_Enriched_Occurrences + Total_Depleted_Occurrences,
        Prop_Enriched = ifelse(Total_Significant_Occurrences > 0, Total_Enriched_Occurrences / Total_Significant_Occurrences, 0),
        Prop_Depleted = ifelse(Total_Significant_Occurrences > 0, Total_Depleted_Occurrences / Total_Significant_Occurrences, 0)
    ) %>%
    ungroup() %>% # Important: Ungroup after rowwise operations
    # Reorder columns for better readability: Name, then sums/proportions, then specific enriched/depleted columns
    select(
        Name,
        Total_Significant_Occurrences,
        Total_Enriched_Occurrences,
        Total_Depleted_Occurrences,
        Prop_Enriched,
        Prop_Depleted,
        starts_with("enriched_in_"),
        starts_with("depleted_in_")
    )


# Print the final combined table
print("\nCombined Table: TFBS Enrichment/Depletion and Summaries")
print(final_combined_table)
write.csv(final_combined_table, "/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/regioneReloaded_sampling1000_ngn2_and_undiffWTC11_remap2022_combined_tfbs_presence_enrichment_depletion.csv", row.names = FALSE)



get_top_enriched_tfbs_abs_zscore <- function(tfbs_names_in_category, alist_region_set_name, n_select, all_sig_hits_df, filter_multi_cell_type_significant = FALSE) {
    if (length(tfbs_names_in_category) == 0) {
        return(character(0))
    }

    temp_all_sig_hits_df <- all_sig_hits_df

    if (filter_multi_cell_type_significant) {
        tfbs_multi_significant_counts <- all_sig_hits_df %>%
            group_by(name) %>%
            summarise(num_sets_significant = n_distinct(Alist_Region_Set)) %>%
            ungroup()

        multi_set_tfbs <- tfbs_multi_significant_counts %>%
            filter(num_sets_significant > 1) %>%
            pull(name)

        temp_all_sig_hits_df <- temp_all_sig_hits_df %>%
            filter(!(name %in% multi_set_tfbs))

        tfbs_names_in_category <- tfbs_names_in_category[!(tfbs_names_in_category %in% multi_set_tfbs)]

        if (length(tfbs_names_in_category) == 0) {
            message(paste0("  Note: All TFBS in category '", alist_region_set_name, "' were filtered out as they were significant in multiple cell types."))
            return(character(0))
        }
    }

    # Step 2: Apply original filtering, sorting by ABSOLUTE norm_zscore, and selection
    filtered_and_sorted <- temp_all_sig_hits_df %>%
        filter(
            Alist_Region_Set == alist_region_set_name,
            name %in% tfbs_names_in_category
        ) %>%
        arrange(desc(abs(norm_zscore))) %>% # <--- CHANGED: Sort by absolute norm_zscore
        head(n_select) %>%
        pull(name)

    return(filtered_and_sorted)
}


# --- Step 1.b: Select 5 TFBS from each category based on MOST ENRICHED (absolute norm_zscore) ---
num_to_select <- 5
selected_tfs_for_plot_by_group <- list()

message("\n--- Selecting Most Enriched TFBS for Plotting (ranked by absolute norm_zscore) ---")

# 1. TFBS unique for "Neuron_upregulating_dCREs"
selected_tfs_for_plot_by_group$unique_Neuron_up <- get_top_enriched_tfbs_abs_zscore(unique_neuron_up_tfbs, "Neuron_upregulating_dCREs", num_to_select, all_significant_hits, filter_multi_cell_type_significant = TRUE)
message(paste0("Selected ", length(selected_tfs_for_plot_by_group$unique_Neuron_up), " most enriched TFBS unique to 'Neuron_upregulating_dCREs' (excluding multi-cell type significant):"))
print(selected_tfs_for_plot_by_group$unique_Neuron_up)


# 2. TFBS unique for "undiffWTC11_upregulating_dCREs"
selected_tfs_for_plot_by_group$unique_undiffWTC11_up <- get_top_enriched_tfbs_abs_zscore(unique_undiffWTC11_up_tfbs, "undiffWTC11_upregulating_dCREs", num_to_select, all_significant_hits, filter_multi_cell_type_significant = TRUE)
message(paste0("Selected ", length(selected_tfs_for_plot_by_group$unique_undiffWTC11_up), " most enriched TFBS unique to 'undiffWTC11_upregulating_dCREs' (excluding multi-cell type significant):"))
print(selected_tfs_for_plot_by_group$unique_undiffWTC11_up)


# 3. TFBS common to both "upregulating" sets
# Ranking by the maximum ABSOLUTE norm_zscore across the two 'upregulating' sets
if (length(common_up_tfbs) > 0) {
    selected_tfs_for_plot_by_group$common_up <- all_significant_hits %>%
        filter(
            Alist_Region_Set %in% c("Neuron_upregulating_dCREs", "undiffWTC11_upregulating_dCREs"),
            name %in% common_up_tfbs
        ) %>%
        group_by(name) %>%
        summarise(max_abs_norm_zscore = max(abs(norm_zscore), na.rm = TRUE)) %>% # <--- CHANGED: Max of absolute norm_zscore
        ungroup() %>%
        arrange(desc(max_abs_norm_zscore)) %>% # <--- CHANGED: Arrange by max absolute score
        head(num_to_select) %>%
        pull(name)
    message(paste0("Selected ", length(selected_tfs_for_plot_by_group$common_up), " most enriched TFBS common to both upregulating sets (ranked by max absolute norm_zscore):"))
    print(selected_tfs_for_plot_by_group$common_up)
} else {
    message("No common TFBS found for upregulating sets to select from.")
    selected_tfs_for_plot_by_group$common_up <- character(0)
}


# 4. TFBS significant for "Neuron_downregulating_dCREs"
selected_tfs_for_plot_by_group$Neuron_down <- get_top_enriched_tfbs_abs_zscore(neuron_down_tfbs, "Neuron_downregulating_dCREs", num_to_select, all_significant_hits, filter_multi_cell_type_significant = TRUE)
message(paste0("Selected ", length(selected_tfs_for_plot_by_group$Neuron_down), " most enriched TFBS from 'Neuron_downregulating_dCREs' (excluding multi-cell type significant):"))
print(selected_tfs_for_plot_by_group$Neuron_down)


# Combine all selected TFBS from all groups into one unique list for the plot
selected_tfs_for_plot <- unique(unlist(selected_tfs_for_plot_by_group))

if (length(selected_tfs_for_plot) == 0) {
    stop("No TFBS were selected for plotting after applying all filters. Please check your criteria.")
}
message(paste0("\nTotal unique TFBS selected for plotting: ", length(selected_tfs_for_plot)))
print(selected_tfs_for_plot)


# --- Step 2: Subset the matrix and set it into an existing copy of a CrosswiseMatrix ---
message("\nPreparing CrosswiseMatrix for plotting selected significant TFBS...")

selected_tfs_for_plot_filtered <- selected_tfs_for_plot[selected_tfs_for_plot %in% colnames(matrix_data)]

if (length(selected_tfs_for_plot_filtered) == 0) {
    stop("None of the specifically selected TFBS were found in the main matrix data. Cannot generate plot.")
}

subset_matrix_data <- matrix_data[, selected_tfs_for_plot_filtered]

message(paste0("Subsetted matrix dimensions: ", paste(dim(subset_matrix_data), collapse = "x")))
message(paste0("Plotting ", ncol(subset_matrix_data), " TFBS across ", nrow(subset_matrix_data), " dCRE sets."))

ngn2_and_undiffWTC11_CrosswiseMatrix_subset <- ngn2_and_undiffWTC11_CrosswiseMatrix
ngn2_and_undiffWTC11_CrosswiseMatrix_subset@matrix$GMat <- subset_matrix_data

if (!is.null(ngn2_and_undiffWTC11_CrosswiseMatrix_subset@matrix$FitCol)) {
    ngn2_and_undiffWTC11_CrosswiseMatrix_subset@matrix$FitCol <- NULL
    message("Removed old column clustering information (FitCol) to allow re-clustering on subset.")
}

# --- Step 3: Plot the resulting genoMatriXeR using plotCrosswiseMatrix() ---
message("\nGenerating plotCrosswiseMatrix for the selected significant TFBS subset...")

subset_plot <- plotCrosswiseMatrix(ngn2_and_undiffWTC11_CrosswiseMatrix_subset)

print(subset_plot)
ggsave("/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/ngn2_and_undiffWTC11_subset_TFBS_crosswise_matrix_plot_abs_enrichment.pdf", plot = subset_plot, width = 10, height = 8)


##### OLD code:
# plot a subset of the data in the matrix form (working)

# Get the matrix of Z-scores (or transformed values)
matrix_data <- regioneReloaded::getMatrix(ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix)
matrix_data <- regioneReloaded::getMatrix(undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix)

dim(matrix_data)

colnames(matrix_data) # These should be your TF names
rownames(matrix_data) # These should be your Alist names

# Get the enrichment values for sequence set 1
enrichment_set1 <- matrix_data["Neuron_upregulating_dCREs", ]
head(enrichment_set1)
# Get the enrichment values for sequence set 2
enrichment_set2 <- matrix_data["Neuron_downregulating_dCREs", ]
head(enrichment_set2)

# Order TFs by enrichment for set 1 (descending)
sorted_tfs_set1 <- names(sort(enrichment_set1, decreasing = TRUE))

# Order TFs by enrichment for set 2 (descending)
sorted_tfs_set2 <- names(sort(enrichment_set2, decreasing = TRUE))

# Select the top 15 enriched TFs for each set
top15_set1_tfs <- head(sorted_tfs_set1, 15)
top15_set2_tfs <- head(sorted_tfs_set2, 15)


# Combine and unique the TF names to avoid duplicates if a TF is in both top 15 lists
selected_tfs <- unique(c(top15_set1_tfs, top15_set2_tfs))

# You might want to ensure that these selected_tfs are indeed in your matrix_data colnames
selected_tfs <- selected_tfs[selected_tfs %in% colnames(matrix_data)]

# Subset the matrix data by selected TFs (columns)
subset_matrix_data <- matrix_data[, selected_tfs]

# Check dimensions of the new matrix
dim(subset_matrix_data)

# Make a copy of the original genoMatriXeR object
ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset <- ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix
undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset <- undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix

# Replace the 'matrix' slot with the subsetted matrix data
# The @matrix slot itself is a list, and the actual matrix is usually under 'GMat' or 'LZM'
# Let's inspect the structure of the @matrix slot first
# str(ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix@matrix)

# Assuming the matrix is stored in the 'GMat' element of the @matrix slot
ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$GMat <- subset_matrix_data
undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$GMat <- subset_matrix_data

# If you also had clustering information stored in @matrix (e.g., FitCol for column clustering),
# you'll need to either remove it (to force plotCrosswiseMatrix to re-cluster if you set clusterize_cols=TRUE)
# or update it if you manually re-cluster. For simplicity, let's remove it and let plotCrosswiseMatrix handle it.

# Remove existing column clustering information if it exists, to allow new clustering on subset
# This depends on how the clustering information is stored. It's often in FitCol/FitRow.
if (!is.null(ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$FitCol)) {
    ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$FitCol <- NULL
}
if (!is.null(undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$FitCol)) {
    undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset@matrix$FitCol <- NULL
}

subset_plot <- plotCrosswiseMatrix(ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix_subset)
ggsave("/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/ngn2_subset_TFBS_crosswise_matrix_plot.pdf", plot = subset_plot, width = 10, height = 8)

# subset_plot <- plotCrosswiseMatrix(undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix_subset)
# ggsave("/home/kisa/coding/80K_MPRA/remap_chipseq_pia/80k_results/undiffWTC11_subset_TFBS_crosswise_matrix_plot.pdf", plot = subset_plot, width = 10, height = 8)


# Extract the two datasets
up_data <- ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_upregulating_dCREs %>%
    mutate(regulation = "increasing transcription") # Add a column to label the dataset
down_data <- ngn2_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_downregulating_dCREs %>%
    mutate(regulation = "decreasing transcription") # Add a column to label the dataset
# up_data <- undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_upregulating_dCREs %>%
#     mutate(regulation = "increasing transcription") # Add a column to label the dataset
# down_data <- undiffWTC11_remap_chip_increasing_decreasing_CrosswiseMatrix@multiOverlaps$Neuron_downregulating_dCREs %>%
#     mutate(regulation = "decreasing transcription") # Add a column to label the dataset

head(up_data)

# Select top x other enriched experiments (not in specific_experiments) for each dataset
n_top = 15
sig_level = 0.1
top5_up <- up_data %>%
    filter(adj.p_value < sig_level & is.finite(norm_zscore)) %>%
    arrange(desc(norm_zscore)) %>%
    slice_head(n = n_top) %>%
    pull(name)

top5_down <- down_data %>%
    filter(adj.p_value < sig_level & is.finite(norm_zscore)) %>%
    arrange(desc(norm_zscore)) %>%
    slice_head(n = n_top) %>%
    pull(name)

head(top5_up)
head(top5_down)

# Filter for specific experiments
filtered_up <- up_data %>%
    filter(name %in% c(top5_up))
filtered_down <- down_data %>%
    filter(name %in% c(top5_down))

# Combine all data
combined_data <- bind_rows(filtered_up, filtered_down)
head(combined_data)
x_limits <- range(c(filtered_up$norm_zscore, filtered_down$norm_zscore), na.rm = TRUE)


# Plot the data
pdf_name_enhancer <- "ngn2_dotplot_enhancer_remap2022_TF_enrichment_sampling.pdf"
pdf_name_enhancer <- "undiffWTC11_dotplot_enhancer_remap2022_TF_enrichment_sampling.pdf"
pdf(file.path(output_directory, pdf_name_enhancer), width = 8, height = 6) # Adjust width and height
ggplot(filtered_up, aes(x = norm_zscore, y = reorder(name, norm_zscore), fill = regulation)) +
    geom_point(aes(size = -log10(adj.p_value)), shape = 21, alpha = 0.6, stroke = NA) + # Ensure fill is mapped and color is NA
    labs(
        x = "Enrichment (normalized z-score)",
        y = "Experiment",
        title = "Enrichment of CREs"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "top",
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 10, margin = margin(r = 5)), # Reduce y-axis text size for density
        axis.text.x = element_text(size = 12) # Adjust x-axis text size
    ) +
    scale_fill_manual(
        values = c("Upregulating" = "#33a02c"), # , "Downregulating" = "#e31a1c"),  # Green for Enhancer, Red for Silencing
        labels = c("Upregulating" = "Enhancer"), # , "Downregulating" = "Silencer"),  # Rename labels
        name = "Regulation"
    ) +
    scale_size_continuous(
        range = c(2, 5), # Adjust dot size range for better visibility
        name = "Significance (-log10 p-value)"
    ) +
    guides(
        fill = guide_legend(override.aes = list(size = 6, shape = 21)) # Increase legend dot size
    ) +
    scale_x_continuous(limits = x_limits) +
    coord_cartesian(clip = "off") # Ensure no points are cut off
dev.off()

# Plot the data
pdf_name_silencer <- "ngn2_dotplot_silencer_remap2022_TF_enrichment_sampling.pdf"
pdf_name_silencer <- "undiffWTC11_dotplot_silencer_remap2022_TF_enrichment_sampling.pdf"
pdf(file.path(output_directory, pdf_name_silencer), width = 8, height = 6) # Adjust width and height
ggplot(filtered_down, aes(x = norm_zscore, y = reorder(name, norm_zscore), fill = regulation)) +
    geom_point(aes(size = -log10(adj.p_value)), shape = 21, alpha = 0.6, stroke = NA) + # Ensure fill is mapped and color is NA
    labs(
        x = "Enrichment (normalized z-score)",
        y = "Experiment",
        title = "Enrichment of CREs"
    ) +
    theme_minimal(base_size = 12) +
    theme(
        legend.position = "top",
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.major.x = element_blank(),
        axis.text.y = element_text(size = 10, margin = margin(r = 5)), # Reduce y-axis text size for density
        axis.text.x = element_text(size = 12) # Adjust x-axis text size
    ) +
    scale_fill_manual(
        values = c("Downregulating" = "#e31a1c"), # Green for Enhancer, Red for Silencing
        labels = c("Downregulating" = "Silencer"), # Rename labels
        name = "Regulation"
    ) +
    scale_size_continuous(
        range = c(2, 5), # Adjust dot size range for better visibility
        name = "Significance (-log10 p-value)"
    ) +
    guides(
        fill = guide_legend(override.aes = list(size = 6, shape = 21)) # Increase legend dot size
    ) +
    scale_x_continuous(limits = x_limits) +
    coord_cartesian(clip = "off") # Ensure no points are cut off
dev.off()
