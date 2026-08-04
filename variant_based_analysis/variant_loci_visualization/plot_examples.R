source("src/plot_variants_v3.R")

# print: 2606: in text
print("2606: in text")
plot_variant(rsid = 'rs1257445811', brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs11635753',   brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs11789013',   brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs1128287',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)

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

# Variant mentioned as interesting from Nick but is rare
print("Variant mentioned as interesting from Nick")
plot_variant(rsid = 'rs776572617',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)
plot_variant(rsid = 'rs11688390',    brain_atac = c('astrocyte','cortical','hippocampus','motor'), show_ngn2_atac = FALSE, show_gtex = FALSE, show_ems = FALSE, arc_lwd = 2.5)


# rs1128287: NC_000005.10:177301456:G:T => associated to NSD1

# GPRIN1: 176610156       176610157
# UIMC1:  177022633       177022634
# PRELID1: 177303799       177303800
# MXD3: 177312757       177312758
# SLC34A1: 177379235       177379236
# PRR7: 177446445       177446446
# TSPAN17: 176647387       176647388
# FGFR4: 177086905       177086906
# NSD1: 177133025       177133026
# RAB24: 177303744       177303745
# PRELID1: 177303799       177303800
# MXD3: 177312757       177312758
# RGS14: 177357924       177357925

# tissues:
# basalganglia_EUR
# cerebellum_EUR
# cortex_AFR
# cortex_EAS
# cortex_EUR
# hippocampus_EUR
# spinalcord_EUR


# 176365487       176365488       0.072033        0.034790887336864       ARL10   cortex_EUR
# 176388751       176388752       -0.078603       0.0212417193329594      HIGD2A  cortex_EUR
# 176447628       176447629       0.080651        0.0180875579597154      FAF2    cortex_EUR
# 177086905       177086906       0.17671 1.9946211478428163e-07  FGFR4   cortex_EUR      Tru
# 177303799       177303800       0.178606        1.4706648810486074e-07  PRELID1 cortex_EUR
# 177312757       177312758       -0.553511       2.813879885546997e-65   MXD3    cortex_EUR
# 177357924       177357925       -0.207742       9.080888899038995e-10   RGS14   cortex_EUR
# 177409576       177409577       -0.081275       0.0172118631949172      F12     cortex_EUR
# 177511274       177511275       -0.077065       0.0239143436238473      DOK3    cortex_EUR

#  177303799       177303800       0.306849        0.0070995852914762      PRELID1 bas
#  177312757       177312758       -0.562023       2.653181539500062e-07   MXD3    bas
#  177446445       177446446       0.273394        0.0168460974486398      PRR7    bas
#  177086905       177086906       0.300775        0.0001031980406812      FGFR4   cer
#  177303744       177303745       0.378673        8.097303155510625e-07   RAB24   cer
#  177303799       177303800       0.265155        0.0006461294657527      PRELID1 cer
#  177312757       177312758       -0.727461       2.4832070068988e-24     MXD3    cer
#  177312757       177312758       -0.394983       0.0020472088051572      MXD3    cor
#  176365487       176365488       0.254256        0.0083994242072561      ARL10   cor
#  177022696       177022697       0.238291        0.0136980092694015      ZNF346  cor
#  177303799       177303800       0.364848        0.00011779684071        PRELID1 cor
#  177357924       177357925       0.230446        0.0172339175415791      RGS14   cor
#  177379235       177379236       0.230445        0.0172346551053893      SLC34A1 cor
#  177446445       177446446       0.322773        0.0007234890848182      PRR7    cor
#  178113532       178113533       0.266928        0.0055769326596707      N4BP3   cor
#  178204533       178204534       0.232344        0.0163127817310758      HNRNPAB cor
#  176447628       176447629       0.080651        0.0180875579597154      FAF2    cor
#  177086905       177086906       0.17671 1.9946211478428163e-07  FGFR4   cortex_EUR
#  177303799       177303800       0.178606        1.4706648810486074e-07  PRELID1 cor
#  177312757       177312758       -0.553511       2.813879885546997e-65   MXD3    cor
#  177357924       177357925       -0.207742       9.080888899038995e-10   RGS14   cor
#  177409576       177409577       -0.081275       0.0172118631949172      F12     cor
#  177133025       177133026       0.305098        0.0168412741066221      NSD1    hip
#  177312757       177312758       -0.45892        0.0002378890900098      MXD3    hip
#  178130996       178130997       -0.364294       0.0040429950796816      RMND5B  hip
#  178204533       178204534       0.342134        0.0071000473078433      HNRNPAB hip