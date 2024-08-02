%% ________________________________________________________________________
%% SMS-SENSE reconstruction for EPI data
%% Congyu Liao, PhD,  cyliao@stanford.edu
%% $03/20/2021
%% ______________________________________________________________________
clear
close all
clc;

addpath utils
addpath imagine
addpath(genpath('./utils/mtimesx_20110223/'))
addpath(genpath('./utils/SENSE_LSQR_Toolbox/'))
addpath(genpath('./utils/Recon_Jon/mrir_toolbox/'))
addpath(genpath('./utils/read_meas_dat__20140924112147/'))

%%
% save path
save_path = './recon/';
mkdir(save_path);

%data path
file_path = './data/';
file_name= 'meas_MID00432_FID18794_fleet_sms3pat3_1p5iso_shift0p75_64dir_avg1.dat';

f_DiffData = [file_path, file_name];


[meas.prot, meas.evp] = read_meas_prot(f_DiffData); % load protocol from "diffusion data"

AccZ = 3;                          % sms factor
AccY = meas.prot.lAccelFactPE;     % PAT factor

PhaseShiftBase = pi;% 0- no caipi shfit; pi-> FOV/2 caipi shift

voxel_size = [meas.prot.dReadoutFOV,meas.prot.dPhaseFOV, meas.prot.dThickness] ./ ...
    [meas.prot.lBaseResolution, meas.prot.lPhaseEncodingLines, meas.prot.sSliceArray_lSize]
esp = meas.prot.iEffectiveEpiEchoSpacing * 1e-6   % echo spacing
start_line = meas.evp.NFirstLin;
N = [meas.prot.lBaseResolution , meas.prot.lPhaseEncodingLines]  % image size

%% do noise prewhitening + coil compression using PAT ref data

if (~exist (strcat(save_path,'sens_gre.mat'),'file'))
    
    iReps =1;
    [meas] = load_SMS_data_CLv2(f_DiffData,iReps,AccZ,meas.evp);
    
    %pat reference
    PAT_ref = meas.patrefscan(:,:,:,:,:,:,:,:,:,:);
    % navigator of PAT ref
    PAT_nav = meas.patrefscan_phascor(:,:,:,:,:,:,:,:,:,:);
    
    
    Ncc= 16; % compress the number of coil to 16.
    PAT_img  = ghost_correct_pat_ref_v1_FLEET(meas.prot, PAT_ref,PAT_nav, AccY);
    [PAT_img, ~,whmcc]= ref_prewht(PAT_img,meas.noiseadjscan, N,Ncc);  % load ref data, coil compression and noise prewhitening
    
    sens_gre = sq(CoilSense_ESPIRIT3d( permute(PAT_img, [1,2,4,3]))); % calculate coil sensitivity using ESPIRIT
    
    sens_gre= permute(sens_gre, [1,2,4,3]);
    sens_gre= single(squeeze(sens_gre));
    save (strcat(save_path,'sens_gre.mat'), 'sens_gre','whmcc', '-V7.3');
else
    load (strcat(save_path,'sens_gre.mat'));
end

%% SMS-SENSE reconstruction
%
num_Rep = meas.evp.RawRep;

for iReps =1:num_Rep;
    [meas] = load_SMS_data_CLv2(f_DiffData,iReps,AccZ,meas.evp);
    clear EPI_data EPI_nav kspace_cor
    % undersampeld EPI Data
    EPI_data = meas.data;
    % navigator of EPI data
    EPI_nav = meas.data_phascor1d;
    
    % ghost correction for epi data
    [~,kspace_cor] = ghost_correct_v2_STD(meas.prot, EPI_data, EPI_nav, AccY,start_line);
    % noise prewhitening +  coil compression
    [kspace_cor] = coilcomp_prewht(kspace_cor,whmcc);
    
    % caipirinha_deblur
    [kspace_cor]= caipi_deblur_CLv1(kspace_cor, meas.prot, meas.evp, AccZ,PhaseShiftBase);
    
    % kyshift and zero-pad for the partial Fourier part
    [kspace_cor]= PF_preprocess(kspace_cor,meas.prot,0);
    
    % detect ky lines of each shot
    [kspace_cor,ky_idx]=detect_kyLines(kspace_cor);
    
    delete(gcp('nocreate'))
    tic
    % non-parallel computation
    % [img_recon] = recon_SMS_data_XCLv3(kspace_cor, ky_idx, sens_gre,AccY, AccZ,PhaseShiftBase);
    % with parallel computation
    show_mercy = 4;
    [img_recon] = recon_SMS_data_XCLv3_parfor(kspace_cor, ky_idx, sens_gre,AccY, AccZ,PhaseShiftBase,show_mercy);
    toc
    
    % apodization to improve SNR and reduce gibbs ringing
    apodization_para = 0.2;
    if apodization_para > 0
        kcomb = mrir_fDFT_freqencode(mrir_fDFT_phasencode(img_recon));
        kapodize = mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(kcomb, 1, apodization_para),  2, apodization_para) ;
        img_recon = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kapodize));
    end
    img_recon = single(img_recon);
    
    
    save (strcat(save_path,'img_recon_Rep',num2str(iReps),'.mat'), 'img_recon');
    
end




