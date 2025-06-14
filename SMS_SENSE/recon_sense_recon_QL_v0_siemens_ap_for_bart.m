% test for XA30
% Qiang Liu based on Congyu and Berkin's code
% SENSE recon, ref data from EPI
% AP data
% March 14 2024

clear;close all;clc;
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/functions_recon_nomapVBVD'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/read_meas_dat__20140924112147'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/FID-A-master'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Joint_Loraks_Toolbox'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Espirit_matlab_only_toolbox'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Harmonization_MB_Recon/SMS_SENSE'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/sms_bart/bart-master/matlab'))

%% extract the ref and img k-space data with mapVBVD, after EPI correction

file_path = '/data/pnlx/home/ql087/data_bwh/For_Hun/Siemens_Prisma_Phantom_Feb_09_2025/2025_02_09_bwh_phantom/';

file_name='meas_MID00059_FID139203_Diffusion_SMS2_R2_AP.dat';


dat = mapVBVD([file_path, file_name]);
meas.prot = read_meas_prot_bb_v0([file_path, file_name]);
prot = meas.prot;

if isempty(prot.aflRegridADCDuration)
    prot.aflRegridADCDuration = dat{2}.hdr.Meas.aflRegridADCDuration(1);
end

[ ref_k, img_k ]=recon_sms_std_xa30_2(dat,prot);

%% extract the sense map from the ref data itself using ESPIRIT

refscan_first_tmp=ref_k;
refscan_first_tmp=padarray(refscan_first_tmp,[0, (146-24)/2], 'both');
img_patref=ifft2c3(refscan_first_tmp);

num_acs = 24;% 36
kernel_size = [6,6];
eigen_thresh = 0.8;			% for mask size 0.8

receive = zeross(size(img_patref));
delete(gcp('nocreate'))
tic
parfor slc_select = 1:s(img_patref,3)
    disp(num2str(slc_select))

    [maps, weights] = ecalib_soft( fft2c( sq(img_patref(:,:,slc_select,:)) ), num_acs, kernel_size, eigen_thresh );

    receive(:,:,slc_select,:) = permute(dot_mult(maps, weights >= eigen_thresh ), [1,2,4,3]);
end
toc
delete(gcp('nocreate'))
sens_gre=permute(receive,[1 2 4 3]);

%% SMS-SENSE recon
PhaseShiftBase=0;
img_k=permute(img_k,[1 3 2 4 5]);
AccY = prot.lAccelFactPE;
AccZ = 2;
NRep=size(img_k,5);
NSlc=size(img_k,4)*AccZ;
kspace_cor=zeros(size(img_k,1), size(img_k,2)*AccY, size(img_k,3),size(img_k,4), size(img_k,5));
kspace_cor(:,2:AccY:end,:,:,:)=img_k;
kspace_cor_tmp=zeros(size(img_k,1), size(img_k,2)*AccY+36, size(img_k,3),size(img_k,4), size(img_k,5));
kdata_dwi=zeros(146*2,146,size(img_k,3), size(img_k,4), size(img_k,5));

for iReps =1:NRep
    pf = prot.ucPhasePartialFourier;
    PE_raw= size(kspace_cor,2);
    PE = max(ceil(PE_raw/AccY)*AccY,prot.lPhaseEncodingLines);
    kspace_cor_tmp(:,:,:,:,iReps)  = mrir_zeropad(kspace_cor(:,:,:,:,iReps),[0 PE-PE_raw 0 0 0 0 0 0 0 0 0],'pre');
    [kspace_cor_tmp(:,:,:,:,iReps), ky_idx]=detect_kyLines_QL_v1(kspace_cor_tmp(:,:,:,:,iReps) );
    delete(gcp('nocreate'))
    tic
    show_mercy = 2;
    [kdata,sens] = recon_SMS_data_XCLv3_parfor_QL_bart(squeeze(kspace_cor_tmp(:,:,:,:,iReps)), ky_idx, sens_gre,AccY, AccZ,PhaseShiftBase,show_mercy);
    toc
    kdata_dwi(:,:,:,:,iReps)=kdata;
    disp(iReps);
end

cd /rfanfs/pnl-zorro/home/ql087/sms_bart/rawdata/

sens=permute(sens,[1 2 14 3 5:13 4]);
writecfl('sens',sens);

kdata_dwi=single(kdata_dwi);
kdata_dwi=permute(kdata_dwi,[1 2 14 3 5:13 4]); 

[FE, PE, ~, COIL, DWI, ~, ~, ~, ~, ~, ~, ~, ~, SLICE] = size(kdata_dwi);
for dwi_idx = 1:DWI
    kdata_slice = kdata_dwi(:, :, :, :, dwi_idx, :, :, :, :, :, :, :, :, :);
    writecfl(sprintf('kdata_%d', dwi_idx), kdata_slice);
end














