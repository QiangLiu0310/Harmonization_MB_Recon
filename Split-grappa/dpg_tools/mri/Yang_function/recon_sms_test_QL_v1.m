% test for XA30
% read rawdata and convert to cfl for bart
% Qiang Liu based on Berkin's code

clear;close all;clc;
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/functions_recon_nomapVBVD'))
addpath(genpath('/data/pnl/home/ql087/data_processing/read_meas_dat__20140924112147'))
addpath(genpath('/data/pnl/home/ql087/data_processing/FID-A-master'))

file_path = '/rfanfs/pnl-zorro/home/ql087/2024_04_09_bwh_prisma_sub3/';
file_name='meas_MID00648_FID37591_Diffusion_SMS2_R2_PA.dat';

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

[raw,coil]=recon_sms_std_xa30_tmp(dat,prot);


% write to cfl for bart

% kspace_data: N-dimensional matrix with slices in the 13th dimension
% coil_sens_maps: N-dimensional matrix with slices in the 13th dimension

kspace_data = permute(raw,[1 3 2 5 6 7 8 9 10 11 12 13 4]);
coil_sens_maps=permute(coil,[1 2 4 5 6 7 8 9 10 11 12 13 3]);

% Write k-space data to CFL file
writecfl('kspace_data', kspace_data);

% Write coil sensitivity maps to CFL file
writecfl('coil_sens_maps', coil_sens_maps);


















