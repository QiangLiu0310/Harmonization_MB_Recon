% test for XA30
% Qiang Liu based on Congyu and Berkin's code
% SENSE recon, ref data from Gre
% March 14 2024

clear;close all;clc;
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/functions_recon_nomapVBVD'))
addpath(genpath('/data/pnl/home/ql087/data_processing/read_meas_dat__20140924112147'))
addpath(genpath('/data/pnl/home/ql087/data_processing/FID-A-master'))
addpath(genpath('/data/pnl/home/ql087/Joint_Loraks_Toolbox'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Espirit_matlab_only_toolbox'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/SMS_SENSE'))

%% extract the ref and img k-space data with mapVBVD, after EPI correction
file_path = '/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_03_15_bwh_prisma_phantom/raw/';
file_name='meas_MID00863_FID23317_Diffusion_SMS2_R2_PA.dat';
dat = mapVBVD([file_path, file_name]);
meas.prot = read_meas_prot_bb_v0([file_path, file_name]);
prot = meas.prot;

if isempty(prot.aflRegridADCDuration)
    prot.aflRegridADCDuration = dat{2}.hdr.Meas.aflRegridADCDuration(1);
end

[ ref_k, img_k ]=recon_sms_std_xa30_2(dat,prot); % ref size 146 24 90 32

%% Gre reference scan
clear ref_k
file_name='meas_MID00849_FID23303_gre.dat';
dat1 = mapVBVD([file_path, file_name]);
ref_k = squeeze( dat1{end}.image()); % oversampled 2
ref_k=permute(ref_k,[1 3 2 4 ]);
tmp=[2:2:90,1:2:90];
[~,or]=sort(tmp);
ref_k=ref_k(:,:,:,or); 
clear or tmp
% remove OS
ref_k = ifft2call(ref_k);
ref_k = fft2call(ref_k(1+end/4:3*end/4,:,:,:,:));
[n(1), n(2), ~, ~] = size(ref_k);
N=[146, 146];
% zero pad patref k-space
Img_gre = ifft2call(padarray(ref_k, (N-n)/2));
ref_k = padarray(ref_k, (N-n)/2);
ref_k=rot90(ref_k,2);
%% extract the sense map from the ref data itself using ESPIRIT

refscan_first_tmp=ref_k;
img_patref=ifft2c3(refscan_first_tmp);
img_patref=permute(img_patref,[ 1 2 4 3]);

num_acs = 24;% 36
kernel_size = [6,6];
eigen_thresh = 0.8;			% for mask size

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
sens_gre=flipdim(sens_gre,2);
%% SMS-SENSE recon
PhaseShiftBase=2*pi/3; % or pi
img_k=permute(img_k,[1 3 2 4]);
AccY = prot.lAccelFactPE;
AccZ = 2;
kspace_cor=zeros(size(img_k,1), size(img_k,2)*AccY, size(img_k,3),size(img_k,4));
kspace_cor(:,2:AccY:end,:,:)=img_k;

for iReps =1:1
    % zero-pad for the partial Fourier part
    pf = prot.ucPhasePartialFourier;
    PE_raw= size(kspace_cor,2);
    PE = max(ceil(PE_raw/AccY)*AccY,prot.lPhaseEncodingLines);
    kspace_cor = mrir_zeropad(kspace_cor,[0 PE-PE_raw 0 0 0 0 0 0 0 0 0],'pre');

    % detect ky lines of each shot
    [kspace_cor,ky_idx]=detect_kyLines_QL_v1(kspace_cor);

    delete(gcp('nocreate'))
    tic
    show_mercy = 2;
    [img_recon] = recon_SMS_data_XCLv3_parfor_QL_v1(kspace_cor, ky_idx, sens_gre,AccY, AccZ,PhaseShiftBase,show_mercy);
    toc

    % apodization to improve SNR and reduce gibbs ringing
    apodization_para = 0.2;
    if apodization_para > 0
        kcomb = mrir_fDFT_freqencode(mrir_fDFT_phasencode(img_recon));
        kapodize = mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(kcomb, 1, apodization_para),  2, apodization_para) ;
        img_recon = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kapodize));
    end

    img_recon = single(img_recon);

end

















