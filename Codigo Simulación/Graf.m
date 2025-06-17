
K = 20;
p = 100;

Setup_results_file=['CFmMIMO_MR_MMSE_UL_K_' num2str(K) '_p_' num2str(p) '.mat'];
% save(Setup_results_file,'SE_MMSE_original','SINR_MMSE_original','SE_MR_original','SINR_MR_original','SE_MMSE_DCC','SINR_MMSE_DCC','SE_MR_DCC','SINR_MR_DCC');
load(Setup_results_file);

