#!/bin/bash

# geral settings:
n_threads = 16

#train the regression model
#-y 3 selects regression (default of -y 0 does classification)
#-t controls the choice of kernel. -t 2 (instead of the default of 4)
# avoids upweighting the central positions.
#The parameter -p controls the “epsilon” used in the support vector
# regression objective (i.e. the tolerated margin of error).
# HEADS UP: don’t confuse this with the -e parameter, which determines
# the epsilon used for the convergence criterion (they are both
# called ‘epsilon’, sigh).
#Depending on the problem, the -c parameter (controlling the cost of
# misprediction) will likely have to be adjusted, i.e. the defaults (10) may
# not work out-of-the-box.

# path to gkmsvm script:
gkmsvm_script = /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/notebooks/rlsgkm/bin/gkmtrain

# train data:
training_fasta = /home/kisa/coding/80K_MPRA/element_variant_analysis_output/gkm_SVM/downsample_without_activity_sequences_n_19k/gkm_svm_train.fa
training_label = /home/kisa/coding/80K_MPRA/element_variant_analysis_output/gkm_SVM/downsample_without_activity_sequences_n_19k/gkm_svm_training_label.tsv
regression_output_name = /home/kisa/coding/80K_MPRA/element_variant_analysis_output/gkm_SVM/downsample_without_activity_sequences_n_19k/regression_downsampled_without_activity
# test data for the predictions:
test_fasta = /home/kisa/coding/80K_MPRA/80K-Analysis/06_variant_element_analysis/notebooks/gkm_svm_test.fa
regression_model_path = "${regression_output_name}.model.txt"
test_predictions_output =
# training:
$gkmsvm_script -T $n_threads -y 3 -t 2 -c 10 -p 0.3 $training_fasta $training_label $regression_output_name

# make predictions on test:
!rlsgkm/bin/gkmpredict -T 16 $test_fasta regression_removed_no_activity.model.txt regression_test_preds_no_activity.txt
# store predictions: