% Split-GRAPPA
% Qiang Liu based on 
% Berkin (XA30), Scott (Split-grappa)
% For Yunqi: coil-by-coil k-space data 
% For Linbo: coil-by-coil image data (complex)

clear;close all;clc;
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/functions_recon_nomapVBVD'))
addpath(genpath('/data/pnl/home/ql087/data_processing/read_meas_dat__20140924112147'))
addpath(genpath('/data/pnl/home/ql087/data_processing/FID-A-master'))
addpath(genpath('/data/pnlx/home/ql087/code/arrShow-develop'))

file_path = '/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_09_14_bwh_prisma_sub5/2024_09_14_bwh_prisma_sub5_raw/';
file_name='meas_MID00696_FID40232_Diffusion_SMS2_R2_AP.dat';
save_path='/data/pnlx/home/ql087/data_processing/2025_for_Yogesh/yunqi/';

tic
    dat = mapVBVD([file_path, file_name]);
toc

tic
    meas.prot = read_meas_prot_bb_v0([file_path, file_name]);
toc

prot = meas.prot;

if isempty(prot.aflRegridADCDuration)
    prot.aflRegridADCDuration = dat{2}.hdr.Meas.aflRegridADCDuration(1);
end

[ k_image ]=recon_sms_std_xa30_3(dat,prot);

save('ap_k_46to62.mat','k_image','-v7.3')

% check the k-space
as(abs(k_image).^0.2)
% check the image
tmp=sqrt(sum(abs(fif(k_image)).^2,3));
as(tmp)