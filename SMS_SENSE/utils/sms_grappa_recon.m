
clear all
close all
clc;

addpath utils
addpath imagine
addpath(genpath('./utils/Recon_Jon/mrir_toolbox/'))
addpath(genpath('./utils/read_meas_dat__20140924112147/'))

%%
% save path
save_path = './recon/test/';
mkdir(save_path);

%data path
file_path = '/local_mount/space/kunkka/1/users/congyu/dwi_recon/';
file_name= 'meas_MID00432_FID18794_fleet_sms3pat3_1p5iso_shift0p75_64dir_avg1.dat';
f_DiffData = [file_path, file_name];


[meas.prot, meas.evp] = read_meas_prot(f_DiffData); % load protocol from "diffusion data"

AccZ = 3;                          % sms factor
AccY = meas.prot.lAccelFactPE;     % PAT factor

PhaseShiftBase = pi;% 0- no capi shfit; pi-> FOV/2 shift

voxel_size = [meas.prot.dReadoutFOV,meas.prot.dPhaseFOV, meas.prot.dThickness] ./ ...
    [meas.prot.lBaseResolution, meas.prot.lPhaseEncodingLines, meas.prot.sSliceArray_lSize]
esp = meas.prot.iEffectiveEpiEchoSpacing * 1e-6   % echo spacing
start_line = meas.evp.NFirstLin;
N = [meas.prot.lBaseResolution *2, meas.prot.lPhaseEncodingLines]  % readout oversampling

%%

%
iReps =1;
[meas] = load_SMS_data(f_DiffData,iReps,AccZ,meas.evp);

% undersampeld EPI Data
EPI_data = meas.data(:,:,:,:,:,:,:,:,:,:);
% navigator of EPI data
EPI_nav = meas.data_phascor1d(:,:,:,:,:,:,:,:,:,:);

%pat reference
PAT_ref = meas.patrefscan(:,:,:,:,:,:,:,:,:,:);
% navigator of PAT ref
PAT_nav = meas.patrefscan_phascor(:,:,:,:,:,:,:,:,:,:);
% SMS reference
SMS_ref = meas.smsrefscan(:,:,:,:,:,:,:,:,:,:);
% navigator of SMS ref
SMS_nav = meas.smsrefscan_phascor(:,:,:,:,:,:,:,:,:,:);

%% do noise prewhitening + coil compression using PAT ref data
if (~exist (strcat(save_path,'sens_gre.mat'),'file'))
    Ncc= 16; % compress the number of coil to 16.
    PAT_img  = ghost_correct_pat_ref_v1_FLEET(meas.prot, PAT_ref,PAT_nav, AccY);
    [PAT_img, PAT_k,whmcc]= ref_prewht(PAT_img,meas.noiseadjscan, N,Ncc);  % load ref data, coil compression and noise prewhitening
    
    sens_gre = sq(CoilSense_ESPIRIT3d( permute(PAT_img, [1,2,4,3]))); % calculate coil sensitivity using ESPIRIT
    
    sens_gre= permute(sens_gre, [1,2,4,3]);
    sens_gre= single(squeeze(sens_gre));
    save (strcat(save_path,'sens_gre.mat'), 'sens_gre','whmcc', '-V7.3');
else
    load (strcat(save_path,'sens_gre.mat'));
end
% noise prewhitening +  coil compression
[EPI_data] = coilcomp_prewht(EPI_data,whmcc);
[EPI_nav] = coilcomp_prewht(EPI_nav,whmcc);
[PAT_ref] = coilcomp_prewht(PAT_ref,whmcc);
[PAT_nav] = coilcomp_prewht(PAT_nav,whmcc);
[SMS_ref] = coilcomp_prewht(SMS_ref,whmcc);
[SMS_nav] = coilcomp_prewht(SMS_nav,whmcc);

%% step1: SMS-GRAPPA reconstruction 

kSize_GRAPPA = [5,5];  % grappa kernel
lambda_tik= eps;  % lambda for grappa regularization
num_acs = [200,140];          % ACS size


prot = meas.prot ;
prot.kSize_GRAPPA = kSize_GRAPPA ;
prot.lambda_tik = lambda_tik;
prot.num_acs = num_acs;
prot.PhaseShiftBase = PhaseShiftBase;
prot.iReps = iReps ;

recon_SMS_GRAPPA(EPI_data, EPI_nav,PAT_k,SMS_ref,SMS_nav,sens_gre,AccZ,meas.evp,prot);


% %% Step 2: background phase estimation
% 
% disp('Step 2:  Estimation backbround phase from slice VCC-GRAPPA reconstruction data')
% [Irecon_SMS_bkgEst,Krecon_SMS_bkgEst,bkgPhase,PhaseDiff_LPF]=estimate_background_phase_v2 (Irecon_EPI_Gfree_fs_comb,Irecon_SMS_ref_fs_comb,Irecon_SMS_ref_fs,SensMap,FilterSize);
% 
% %2.2 Estimation backbround phase of PAT Ref------------------------------------------
% % for feet acquisition, use estimated phase from sms-ref
% [Irecon_PAT_bkgEst,Krecon_PAT_bkgEst]=estimate_background_phase_fleet (Irecon_EPI_Gfree_fs_comb,Irecon_PAT_acs_comb,Irecon_PAT_acs,SensMap,PhaseDiff_LPF,Irecon_SMS_bkgEst);
% 





