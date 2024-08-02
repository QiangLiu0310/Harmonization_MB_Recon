clear
close all
clc;
%%
addpath utils
addpath imagine
addpath(genpath('./utils/mtimesx_20110223/'))
addpath(genpath('./utils/SENSE_LSQR_Toolbox/'))
addpath(genpath('./utils/Recon_Jon/mrir_toolbox/'))
addpath(genpath('./utils/read_meas_dat__20140924112147/'))

%%
load 'toCL.mat'

%%
k = fftc(k_ap,4);
k = permute(k,[1 2 4 3]);

% detect ky lines of each shot
[~,ky_idx,kz_idx]=detect_kykzLines(k);


[img_recon] = recon_3DSense(k, ky_idx,kz_idx, sens_svd);

