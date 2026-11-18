
########## mean ##########
grep Acc PPLM-PPI.inter_intra_attn_embed.mean.round_1.cv_*.recording | sort -rn -k14 | head -n 5 | awk -F '.recording' '{print $1}' > 000
grep Acc PPLM-PPI.inter_intra_attn_embed.mean.round_1.cv_*.recording | sort -rn -k14 | head -n 5 | awk -F 'Epoch: ' '{print $2}' | awk -F '/' '{print "model_"$1".pkl"}' > 001
paste -d '/' 000 001 > PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list_best_human_all
python /mnt/rna01/junl/PPLM/PPLM/ppi_training/training.PPI.single_gpu.testing_both.20250912.multi-model.py PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list_best_human_all --mode PPI_inter_intra_attn_embed_mean > PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list_best_human_all.results

grep Acc PPLM-PPI.inter_intra_attn_embed.mean.round_2.cv_*.recording | sort -rn -k14
grep === PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list_best_human_all.results  | awk '{print $8" "$9" "$10" "$11" "$15" "$16}'


########## max ##########
grep Acc PPLM-PPI.inter_intra_attn_embed.max.round_1.cv_*.recording | sort -rn -k14 | head -n 5 | awk -F '.recording' '{print $1}' > 000
grep Acc PPLM-PPI.inter_intra_attn_embed.max.round_1.cv_*.recording | sort -rn -k14 | head -n 5 | awk -F 'Epoch: ' '{print $2}' | awk -F '/' '{print "model_"$1".pkl"}' > 001
paste -d '/' 000 001 > PPLM-PPI.inter_intra_attn_embed.max.round_1.top5_list_best_human_all
python /mnt/rna01/junl/PPLM/PPLM/ppi_training/training.PPI.single_gpu.testing_both.20250912.multi-model.py PPLM-PPI.inter_intra_attn_embed.max.round_1.top5_list_best_human_all --mode PPI_inter_intra_attn_embed_max > PPLM-PPI.inter_intra_attn_embed.max.round_1.top5_list_best_human_all.results

grep === PPLM-PPI.inter_intra_attn_embed.max.round_1.top5_list_best_human_all.results  | awk '{print $8" "$9" "$10" "$11" "$15" "$16}'


########## min ##########
grep Acc PPLM-PPI.inter_intra_attn_embed.min.round_1.cv_*.recording | sort -rn -k14 | head -n 5 | awk -F '.recording' '{print $1}' > 000
grep Acc PPLM-PPI.inter_intra_attn_embed.min.round_1.cv_*.recording | sort -rn -k14 | head -n 5 | awk -F 'Epoch: ' '{print $2}' | awk -F '/' '{print "model_"$1".pkl"}' > 001
paste -d '/' 000 001 > PPLM-PPI.inter_intra_attn_embed.min.round_1.top5_list_best_human_all
python /mnt/rna01/junl/PPLM/PPLM/ppi_training/training.PPI.single_gpu.testing_both.20250912.multi-model.py PPLM-PPI.inter_intra_attn_embed.min.round_1.top5_list_best_human_all --mode PPI_inter_intra_attn_embed_min > PPLM-PPI.inter_intra_attn_embed.min.round_1.top5_list_best_human_all.results

grep === PPLM-PPI.inter_intra_attn_embed.min.round_1.top5_list_best_human_all.results  | awk '{print $8" "$9" "$10" "$11" "$15" "$16}'


########################################## combine multiple pooling strateties ##########################################
########## mean_max ##########
python /mnt/rna01/junl/PPLM/PPLM/ppi_training/training.PPI.single_gpu.testing_both.20250912.multi-model.multi_pooling.py --mean_path PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list_best_human_all --max_path PPLM-PPI.inter_intra_attn_embed.max.round_1.top5_list_best_human_all > PPLM-PPI.inter_intra_attn_embed.combine_mean_max.round_1.top5_list_best_human_all.results

grep === PPLM-PPI.inter_intra_attn_embed.combine_mean_max.round_1.top5_list_best_human_all.results  | awk '{print $8" "$9" "$10" "$11" "$15" "$16}'


########## mean_max_min ##########
python /mnt/rna01/junl/PPLM/PPLM/ppi_training/training.PPI.single_gpu.testing_both.20250912.multi-model.multi_pooling.py --mean_path PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list_best_human_all --max_path PPLM-PPI.inter_intra_attn_embed.max.round_1.top5_list_best_human_all --min_path PPLM-PPI.inter_intra_attn_embed.min.round_1.top5_list_best_human_all > PPLM-PPI.inter_intra_attn_embed.combine_mean_max_min.round_1.top5_list_best_human_all.results

grep === PPLM-PPI.inter_intra_attn_embed.combine_mean_max_min.round_1.top5_list_best_human_all.results  | awk '{print $8" "$9" "$10" "$11" "$15" "$16}'


########## mean_max ##########
python /mnt/rna01/junl/PPLM/PPLM/ppi_training/training.PPI.single_gpu.testing_both.20250912.multi-model.multi_pooling.py --mean_path PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list_best_human_all --max_path PPLM-PPI.inter_intra_attn_embed.max.round_1.top5_list_best_human_all > PPLM-PPI.inter_intra_attn_embed.combine_mean_max.round_1.top5_list_best_human_all.results

grep === PPLM-PPI.inter_intra_attn_embed.combine_mean_max.round_1.top5_list_best_human_all.results  | awk '{print $8" "$9" "$10" "$11" "$15" "$16}'


########## mean_min ##########
python /mnt/rna01/junl/PPLM/PPLM/ppi_training/training.PPI.single_gpu.testing_both.20250912.multi-model.multi_pooling.py --mean_path PPLM-PPI.inter_intra_attn_embed.mean.round_1.top5_list_best_human_all --min_path PPLM-PPI.inter_intra_attn_embed.min.round_1.top5_list_best_human_all > PPLM-PPI.inter_intra_attn_embed.combine_mean_min.round_1.top5_list_best_human_all.results

grep === PPLM-PPI.inter_intra_attn_embed.combine_mean_min.round_1.top5_list_best_human_all.results  | awk '{print $8" "$9" "$10" "$11" "$15" "$16}'
