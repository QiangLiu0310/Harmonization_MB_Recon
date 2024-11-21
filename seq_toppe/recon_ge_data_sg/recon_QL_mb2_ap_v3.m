
% GE SA file recon
% full version
% PA data
% split-GRAPPA
% Jul/11/2024
% use 3-line scan as the ref scan
% Qiang Liu

clear;close all;clc;

addpath(genpath('/data/pnl/home/ql087/orchestra-sdk-2.1-1.matlab'))
addpath(genpath('/data/pnl/home/ql087/arrShow-develop'));
addpath(genpath('/data/pnl/home/ql087/Bruker_2022'));
addpath(genpath('/data/pnl/home/ql087/functions_recon'));
addpath(genpath('/data/pnl/home/ql087/Pulseq_Mprage_Recon_Toolbox'));
addpath(genpath('/data/pnl/home/ql087/freesurfer'));
addpath(genpath('/data/pnl/home/ql087/qiang_gSlider_data/lq/Harmonization_MB_Recon/SMS_SENSE'))


data_path='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_08_21_ge_sub1/Exam18233/Series3/'; % image

currentDir = pwd;

if ~strcmp(currentDir, data_path)
    cd(data_path);
end
disp('Currently in directory ref scan');


%% read GE MB2 data

tmp=strcat(data_path,'ScanArchive_LONGWOOD30MR2_20240821_214341454.h5');
pfile = fullfile(tmp);
archive = GERecon('Archive.Load', pfile);

% Scan parameters
SA_head =  archive.DownloadData.rdb_hdr_rec;
mydiff_head.sliceorder1 = archive.DownloadData.rdb_hdr_series.start_ras;
mydiff_head.sliceorder2 = archive.DownloadData.rdb_hdr_series.end_ras;
mydiff_head.MB = SA_head.rdb_hdr_mb_factor;
mydiff_head.total_slice = SA_head.rdb_hdr_nslices;
mydiff_head.pass = SA_head.rdb_hdr_npasses;
mydiff_head.diffdir= SA_head.rdb_hdr_numdifdirs;
mydiff_head.logic_slice = mydiff_head.total_slice / mydiff_head.pass;

for i = 1:mydiff_head.logic_slice
    info = GERecon('Archive.Info', archive, 1, i);
    corners(i)= info.Corners;
    a(i,2) =i;
    a(i,1)= corners(i).UpperLeft(3);
    b= sortrows(a,1);
end

mydiff_head.Orientation = info.Orientation;
% now need sort the cornerpoint in order
for i = 1:mydiff_head.logic_slice
    mydiff_head.corners(i)= corners(b(i,2));
end

xRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
yRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres - 1;
stop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.stop_rcv;
start = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.start_rcv;
nChannels = stop - start + 1;
phs = archive.Passes; % number of phases

pass = 1;
zRes = archive.SlicesPerPass(pass);


% load('/data/pnlx/home/ql087/data_processing/2024_recon_reproducibility/sub5/ge/scan1/split_grappa/ap_scan1_sg.mat')
% img=rot90(img);
% img=flipdim(img,2);

load('/data/pnlx/home/ql087/data_processing/2024_recon_reproducibility/sub1/ge/bart_L1_wavelet/scan1_acs_debug/ap_scan1_bart.mat')
img=flipdim(abs(img_bart),1);

% apodization to improve SNR and reduce gibbs ringing
apodization_para = 0.2;
kcomb = mrir_fDFT_freqencode(mrir_fDFT_phasencode(img));
kapodize = mrir_filter_raw_apodize_1d( mrir_filter_raw_apodize_1d(kcomb, 1, apodization_para),  2, apodization_para) ;
im = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kapodize));
img=single(abs(img));



ap=load_nifti('/data/pnlx/home/ql087/data_processing/2024_recon_reproducibility/sub1/ge/bart_L1_wavelet/scan1_acs_debug/scan1_ap.nii.gz');
ref=ap.vol;

gradwarpImage=zeros(146,146,90,62);

for phase=1:62
    for slice=1:90
        magnitudeImage=img(:,:,slice,phase);
        kissoffViews = archive.DownloadData.rdb_hdr_rec.rdb_hdr_kissoff_views;
        magnitudeImage(:,1:kissoffViews) = 0;
        magnitudeImage(:,(end-kissoffViews+1):end) = 0;
        corners= mydiff_head.corners(slice);
        fieldAdjusted = GERecon('Epi.RealtimeFieldAdjustment', magnitudeImage, slice, phase+1);
        gradwarpImage(:,:,slice,phase) = GERecon('Gradwarp', fieldAdjusted, corners, 'XRMW');
    end
end

gradwarpImage=flipdim(gradwarpImage,2);
gradwarpImage=flipdim(gradwarpImage,1);

ap.vol=gradwarpImage*0.8e-3; % different between scan1 and scan2?
save_nifti(ap,'scan1_ap_new.nii.gz')










