% GE ScanArchive file
% Qiang Liu based on Congyu and Berkin's code
% SENSE recon, ref data from EPI
% AP data
% July 17th 2024

clear;close all;clc;
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/functions_recon_nomapVBVD'))
addpath(genpath('/data/pnl/home/ql087/Joint_Loraks_Toolbox'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Espirit_matlab_only_toolbox'))
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/Harmonization_MB_Recon/SMS_SENSE'))
addpath(genpath('/data/pnl/home/ql087/orchestra-sdk-2.1-1.matlab'))
addpath(genpath('/data/pnl/home/ql087/arrShow-develop'));
addpath(genpath('/data/pnl/home/ql087/Bruker_2022'));
addpath(genpath('/data/pnl/home/ql087/functions_recon'));
addpath(genpath('/data/pnl/home/ql087/Pulseq_Mprage_Recon_Toolbox'));
addpath(genpath('/rfanfs/pnl-zorro/home/ql087/sms_bart/bart-master'));



%% extract the ref and img k-space data
data_path='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series8/'; % image
%% read GE MB2 data

cd(data_path);
disp('Currently in directory rawdata');
tmp=strcat(data_path,'ScanArchive_LONGWOOD30MR2_20240409_213322272.h5');
pfile = fullfile(tmp);
archive = GERecon('Archive.Load', pfile);
xRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
yRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres - 1;
stop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.stop_rcv;
start = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.start_rcv;
nChannels = stop - start + 1;
phs = archive.Passes; % number of phases

pass = 1;
zRes = archive.SlicesPerPass(pass);

kspace = complex(zeros(xRes,  nChannels, yRes, zRes));

for i = 1:archive.ControlCount

    control = GERecon('Archive.Next', archive);

    if isfield(control, 'Data')
        kspace(:,:,1:1:end,i) = squeeze(control.Data);
    else
        fprintf('%d\n', i);
    end

end
GERecon('Archive.Close', archive);
kspace=single(kspace);
tmp=[1:2:90 91:size(kspace,4)];
kspace=kspace(:,:,:,tmp);

%% image parameters
NFreq_outres=146;
NFreq_inres=326;
Ncoils=44;
NLocPhz=3;
NImgLin=54;
NSlc=90;
hdr.NSlc=NSlc;
hdr.R=2;
hdr.SMS=2;
hdr.FOVshift=2; % it is 2
PhaseShift=2*pi/hdr.FOVshift;
Ngroup=NSlc/hdr.SMS;
NRep=62;
load('/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series3/prot.mat')

%% load rawdata
raw=kspace(:,:,:,Ngroup+1:end);
ref=kspace(:,:,4:6,1:Ngroup);
ref=permute(ref,[1 3 2 4]);
ref=repmat(ref,[1 1 1 1 NRep]);
clear kspace

raw=permute(raw,[1 3 2 4]);
raw=reshape(raw,[326 57 44 46 62]);
raw=raw(:,:,:,2:end,:);

% ref=raw(:,1:3,:,:,:);
raw=raw(:,4:end,:,:,:);
raw(:,2:2:end,:,:,:)=flipdim(raw(:,2:2:end,:,:,:),1);
ref(:,2:2:end,:,:,:)=flipdim(ref(:,2:2:end,:,:,:),1);
clear tmp

%% load vrgf dat0
tmp=strcat(data_path, 'vrgf.dat');
fid = fopen(tmp, 'r');
data = fread(fid, inf, 'float'); % 'float' specifies the data type
fclose(fid);
v=reshape(data,[326,146]);
clear data
%% Image k-space data

sz_phasecor =[NFreq_outres NLocPhz Ncoils   Ngroup   NRep];
sz_image =   [NFreq_outres  NImgLin Ncoils   Ngroup   NRep];

% image data extraction
vrgf_phasecor =reshape(v'*reshape(ref,[NFreq_inres Ncoils*NLocPhz*Ngroup*NRep]), sz_phasecor);
vrgf_phasecor=permute(vrgf_phasecor,[1 3 2 4 5]);
% imgscanPC dimension:
vrgf_phasecor_ROp_ROn(:,:,1,:,:)= squeeze((vrgf_phasecor(:,:,1,:,:)+vrgf_phasecor(:,:,3,:,:))/2);
vrgf_phasecor_ROp_ROn(:,:,2,:,:)= squeeze( vrgf_phasecor(:,:,2,:,:));

vrgf_image=reshape(v'*reshape(raw,[NFreq_inres Ncoils*NImgLin*Ngroup*NRep]), sz_image);
vrgf_image=permute(vrgf_image,[1 3 2 4 5]);
Kimage= PhaseCorrect_yang(vrgf_phasecor_ROp_ROn,vrgf_image);
Kimage=permute(Kimage,[1 3 2 4 5]);

sli_idx=[1:2:45 2:2:45];
[~,or]=sort(sli_idx);
Kimage_short=Kimage(:,:,:,or,:);
clear Kimage
img_k=Kimage_short;
clear t_trgt Kimage_short

%% extract the sense map from the ref data itself using ESPIRIT

 load('/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series8/surfaceImages.mat')
 img_patref=permute(surfaceImages,[1 2 4 3]);

num_acs = 24;
kernel_size = [6,6]; % [6 6]
eigen_thresh = 0.8;

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
% sens_gre=flipdim(sens_gre,2);
sens_gre=flipdim(sens_gre,1);
%% SMS-SENSE recon
PhaseShiftBase=pi; % 0 or pi
tmp=zeros(146,55,44,45,62);
tmp(:,2:end,:,:,:)=img_k;
img_k=tmp; clear tmp;
AccY = prot.lAccelFactPE;
AccZ = 2;
kspace_cor=zeros(size(img_k,1), size(img_k,2)*AccY, size(img_k,3),size(img_k,4), size(img_k,5));
kspace_cor(:,2:AccY:end,:,:,:)=img_k;
kspace_cor_tmp=zeros(size(img_k,1), size(img_k,2)*AccY+36, size(img_k,3),size(img_k,4), size(img_k,5));
% img_recon=zeros(size(img_k,1), size(img_k,1), NSlc, NRep);
kdata_dwi=zeros(146*2,146,size(img_k,3), size(img_k,4), size(img_k,5));

for iReps =1:NRep
    % zero-pad for the partial Fourier part
    pf = prot.ucPhasePartialFourier;
    PE_raw= size(kspace_cor,2);
    PE = max(ceil(PE_raw/AccY)*AccY,prot.lPhaseEncodingLines);
    kspace_cor_tmp(:,:,:,:,iReps) = mrir_zeropad(kspace_cor(:,:,:,:,iReps),[0 PE-PE_raw 0 0 0 0 0 0 0 0 0],'pre');
    [kspace_cor_tmp(:,:,:,:,iReps), ky_idx]=detect_kyLines_QL_v1(kspace_cor_tmp(:,:,:,:,iReps) );
%     delete(gcp('nocreate'))
    tic
    show_mercy = 2;
    [kdata,sens] = recon_SMS_data_XCLv3_parfor_QL_bart(squeeze(kspace_cor_tmp(:,:,:,:,iReps)) , ky_idx, sens_gre,AccY, AccZ,PhaseShiftBase,show_mercy);
    toc
    kdata_dwi(:,:,:,:,iReps)=kdata;
    disp(iReps)
end

cd /rfanfs/pnl-zorro/home/ql087/sms_bart/rawdata/
sens=permute(sens,[1 2 14 3 5:13 4]);
writecfl('sens',sens);
kdata_dwi=permute(kdata_dwi,[1 2 14 3 5:13 4]);
[FE, PE, ~, COIL, DWI, ~, ~, ~, ~, ~, ~, ~, ~, SLICE] = size(kdata_dwi);
for dwi_idx = 1:DWI
    kdata_slice = kdata_dwi(:, :, :, :, dwi_idx, :, :, :, :, :, :, :, :, :);
    writecfl(sprintf('kdata_%d', dwi_idx), kdata_slice);
end

















