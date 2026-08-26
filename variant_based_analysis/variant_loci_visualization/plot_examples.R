source("scripts/plot_variants.R")

print("Basis for Figure 4b of the manuscript:")
plot_variant(rsid = 'rs1128287',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
print("Other examples mentioned in the manuscript:")
plot_variant(rsid = 'rs1257445811', brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs11635753',   brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs11789013',   brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)

# potentially interesting based on GWAS associations
print("potentially interesting based on GWAS associations")
plot_variant(rsid = 'rs7970847',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs17572795',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs12203592',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs10903341',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs7236461',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)

# potentially interesting based on PheWAS associations
print("potentially interesting based on PheWAS associations")
plot_variant(rsid = 'rs45575136',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)

print("Other examples of interest")
plot_variant(rsid = 'rs776572617',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs11688390',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)