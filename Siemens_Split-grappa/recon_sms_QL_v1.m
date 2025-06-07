% Split-GRAPPA
% Qiang Liu based on 
% Berkin (XA30), Scott (Split-grappa)
% For Linbo: coil-by-coil image data (complex)

clear;close all;clc;
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/functions_recon_nomapVBVD'))
addpath(genpath('/data/pnl/home/ql087/data_processing/read_meas_dat__20140924112147'))
addpath(genpath('/data/pnl/home/ql087/data_processing/FID-A-master'))
addpath(genpath('/data/pnlx/home/ql087/code/arrShow-develop'))

file_path = '/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_09_14_bwh_prisma_sub5/2024_09_14_bwh_prisma_sub5_raw/';
file_name='meas_MID00696_FID40232_Diffusion_SMS2_R2_AP.dat';
save_path='/data/pnlx/home/ql087/data_processing/2025_for_Yogesh/linbo/';

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

[ image ]=recon_sms_std_xa30_1(dat,prot);

% cd (save_path)
% save('cplx_image_multicoil_46to62.mat','image','-v7.3') % this line might give an error, but it could be saved sperately...
% 
% % check the image
% tmp=abs(image);
% as(tmp)