% Split-GRAPPA
% Qiang Liu based on 
% Berkin (XA30), Scott (Split-grappa)
% For Linbo: coil-by-coil image data (complex)

clear;close all;clc;
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/functions_recon_nomapVBVD'))
addpath(genpath('/data/pnl/home/ql087/data_processing/read_meas_dat__20140924112147'))
addpath(genpath('/data/pnl/home/ql087/data_processing/FID-A-master'))
addpath(genpath('/data/pnlx/home/ql087/code/arrShow-develop'))

file_path = '/data/pnlx/home/ql087/data_bwh/2025-06-14-3/';
file_name='2025-06-14-223317.dat';
save_path='/scratch/home/ql087/data_processing/for_yogesh/noise_06_14_2025_bwh_invivo/';

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

[ image ]=recon_sms_std_xa30_7(dat,prot);

