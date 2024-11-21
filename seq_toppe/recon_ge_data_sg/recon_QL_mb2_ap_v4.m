% GE SA file recon
% full version
% Test SENSE (it uses parfor so it is fast to test)
% Jul/11/2024
% use 3-line scan as the ref scan
% Qiang Liu

clear;close all;clc;

addpath(genpath('/data/pnl/home/ql087/orchestra-sdk-2.1-1.matlab'))
addpath(genpath('/data/pnl/home/ql087/arrShow-develop'));
addpath(genpath('/data/pnl/home/ql087/Bruker_2022'));
addpath(genpath('/data/pnl/home/ql087/functions_recon'));
addpath(genpath('/data/pnl/home/ql087/Pulseq_Mprage_Recon_Toolbox'));

data_path='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_03_28_ge_sub2/Exam16351/Series2/'; % image
data_path1='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_03_28_ge_sub2/Exam16351/Series5/'; % ref

currentDir = pwd;

if ~strcmp(currentDir, data_path1)
    cd(data_path);
end
disp('Currently in directory ref scan');


%% read GE MB2 data

cd(data_path);
disp('Currently in directory rawdata');

tmp=strcat(data_path,'ScanArchive_LONGWOOD30MR2_20240329_205337267.h5');
pfile = fullfile(tmp);
archive = GERecon('Archive.Load', pfile);

% Scan parameters
xRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
yRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres - 1;
stop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.stop_rcv;
start = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.start_rcv;
nChannels = stop - start + 1;
phs = archive.Passes; % number of phases

pass = 1;
zRes = archive.SlicesPerPass(pass);


% Scan parameters
SA_head =  archive.DownloadData.rdb_hdr_rec;
mydiff_head.MB = SA_head.rdb_hdr_mb_factor;
mydiff_head.total_slice = SA_head.rdb_hdr_nslices;
mydiff_head.pass = SA_head.rdb_hdr_npasses;
mydiff_head.diffdir= SA_head.rdb_hdr_numdifdirs;
mydiff_head.logic_slice = mydiff_head.total_slice / mydiff_head.pass;
mydiff_head.b0pass= mydiff_head.pass-mydiff_head.diffdir;
mydiff_head.tensor= SA_head.rdb_hdr_user11;
mydiff_head.ms= SA_head.rdb_hdr_ileaves;
mydiff_head.ExamNumber= archive.ExamNumber;
mydiff_head.pcref_start= SA_head.rdb_hdr_pcref_start;
mydiff_head.pcref_stop= SA_head.rdb_hdr_pcref_stop;
mydiff_head.noise_cal = '';
mydiff_head.asset_cal = '';
mydiff_head.raw_cal = '';
mydiff_head.slthick= archive.DownloadData.rdb_hdr_image.slthick;
mydiff_head.slspace= archive.DownloadData.rdb_hdr_image.scanspacing;
mydiff_head.fov= archive.DownloadData.rdb_hdr_image.dfov;
mydiff_head.dx= archive.DownloadData.rdb_hdr_image.dim_X;
mydiff_head.dy= archive.DownloadData.rdb_hdr_image.dim_Y;
mydiff_head.day= SA_head.rdb_hdr_da_yres -1;
mydiff_head.ctr_R= archive.DownloadData.rdb_hdr_image.ctr_R;
mydiff_head.ctr_A= archive.DownloadData.rdb_hdr_image.ctr_A;
mydiff_head.ctr_S= archive.DownloadData.rdb_hdr_image.ctr_S;
mydiff_head.tableposition= archive.DownloadData.rdb_hdr_series.tablePosition;
mydiff_head.im_size= SA_head.rdb_hdr_im_size;
mydiff_head.RawHeader = SA_head;
mydiff_head.se_no = archive.SeriesNumber;
mydiff_head.se_desc= archive.DownloadData.rdb_hdr_series.se_desc;
mydiff_head.bval = SA_head.rdb_hdr_bvalstab(1);
mydiff_head.sliceorder1 = archive.DownloadData.rdb_hdr_series.start_ras;
mydiff_head.sliceorder2 = archive.DownloadData.rdb_hdr_series.end_ras;
mydiff_head.pepolar= archive.DownloadData.rdb_hdr_image.ihpepolar;
mydiff_head.phasefov=SA_head.rdb_hdr_phase_scale;

mydiff_head.geo_slice=mydiff_head.logic_slice/mydiff_head.MB;
mydiff_head.ms=1;

ph=1;
count=0;
for i=1:archive.ControlCount

    currentControl = GERecon('Archive.Next', archive);
    count=count +1;
    if i==1
        [xres ch yres]=size(currentControl.Data);
        kspace_temp= single(zeros(xres, ch, yres*mydiff_head.ms,  mydiff_head.geo_slice, mydiff_head.pass));
        kspace=complex(kspace_temp);
    end
    if((currentControl.opcode==6) || (currentControl.opcode==14))
        %  count=count+1
        if ((currentControl.opcode==14) && (currentControl.frameType == 2))
            %currentControl.bValueIndex+1
            curr_dir= currentControl.diffDirIndex+1+ mydiff_head.b0pass;
            kspace(:,:,currentControl.viewNum:currentControl.viewSkip:1,currentControl.sliceNum+1,curr_dir)=currentControl.Data;
        else
            kspace(:,:,currentControl.viewNum:currentControl.viewSkip:1,currentControl.sliceNum+1,ph)=currentControl.Data;
        end
        if((currentControl.sliceNum+1== mydiff_head.geo_slice) && ((abs(currentControl.viewSkip)*currentControl.numViews) == (currentControl.viewNum-1+abs(currentControl.viewSkip))) )
            ph=ph+1;
        end
    end
end
kspace = permute(kspace,[1 3 2 4 5]);
clear kspace_temp

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
% raw=kspace(:,:,:,Ngroup+1:end);
ref=kspace(:,:,4:6,1:Ngroup);
ref=permute(ref,[1 3 2 4]);
ref=repmat(ref,[1 1 1 1 NRep]);

% raw=permute(raw,[1 3 2 4]);
% raw=reshape(raw,[326 57 44 46 62]);
% raw=raw(:,:,:,2:end,:);

raw=kspace(:,end:-1:1,:,:,2:end);
ref=kspace(:,end:-1:1,:,:,1);
ref=ref(:,4:6,:,:);
ref=repmat(ref,[1 1 1 1 NRep]);
raw=raw(:,4:end,:,:,:);

raw(:,2:2:end,:,:,:)=flipdim(raw(:,2:2:end,:,:,:),1);
ref(:,2:2:end,:,:,:)=flipdim(ref(:,2:2:end,:,:,:),1);

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

load('/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_03_28_ge_sub2/Exam16351/Series2/surfaceImages.mat')
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
img_recon=zeros(size(img_k,1), size(img_k,1), NSlc, NRep);
for iReps =1:NRep
    % zero-pad for the partial Fourier part
    pf = prot.ucPhasePartialFourier;
    PE_raw= size(kspace_cor,2);
    PE = max(ceil(PE_raw/AccY)*AccY,prot.lPhaseEncodingLines);
    kspace_cor_tmp(:,:,:,:,iReps) = mrir_zeropad(kspace_cor(:,:,:,:,iReps),[0 PE-PE_raw 0 0 0 0 0 0 0 0 0],'pre');

    % detect ky lines of each shot
    [kspace_cor_tmp(:,:,:,:,iReps), ky_idx]=detect_kyLines_QL_v1(kspace_cor_tmp(:,:,:,:,iReps) );

    delete(gcp('nocreate'))
    tic
    show_mercy = 2;
    [img_recon(:,:,:,iReps) ] = recon_SMS_data_XCLv3_parfor(squeeze(kspace_cor_tmp(:,:,:,:,iReps)) , ky_idx, sens_gre,AccY, AccZ,PhaseShiftBase,show_mercy);
    toc

    disp(iReps)
end

