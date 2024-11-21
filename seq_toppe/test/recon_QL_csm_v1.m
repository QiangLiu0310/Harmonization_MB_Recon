% GE SA file recon
% Qiang Liu 
% Jul/11/2024
 clear;close all;clc;
 addpath(genpath('/rfanfs/pnl-zorro/home/ql087/qiang_gSlider_data/lq/functions_recon_nomapVBVD'))
 data_path='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series5/';
%% image parameters
NFreq_outres=146;
NFreq_inres=330;
Ncoils=44;
NLocPhz=3;
NImgLin=73;
NSlc=90;
SMS=2;
Ngroup=NSlc*2;%NSlc/SMS
NRep=1;
%% load rawdata
tmp=strcat(data_path, 'kspace.mat');
raw=load(tmp);
raw=raw.kspace;
raw=raw(:,:,1:2:end,:);
raw(:,:,2:2:end,:)=flipdim(raw(:,:,2:2:end,:),1);
raw=permute(raw,[1 3 2 4]);
tmp=strcat(data_path, 'kspace_ref.mat');
ref=load(tmp);
ref=ref.kspace;
ref(:,:,3:4:end,:)=flipdim(ref(:,:,3:4:end,:),1);
ref=permute(ref,[1 3 2 4]);
ref=ref(:,73:2:77,:,:);
clear kspace tmp

%% flip every second line

%% load vrgf dat0
tmp=strcat(data_path, 'vrgf.dat');
fid = fopen(tmp, 'r');
data = fread(fid, inf, 'float'); % 'float' specifies the data type
fclose(fid);
v=reshape(data,[330,146]);
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
vrgf_image=permute(vrgf_image,[1 3 2 4]);
Kimage= PhaseCorrect_yang(vrgf_phasecor_ROp_ROn,vrgf_image);

tmp=permute(Kimage,[1 3 2 4]);
tmp1=zeros(146,146,44,90);
tmp1(:,1:2:end,:,:)=tmp(:,1:1:end,:,1:90);
tmp1(:,2:2:end,:,:)= -tmp(:,1:1:end,:,91:180);
as(squeeze(sos(fft2cc(tmp1),3)))
clear tmp
tmp1=single(tmp1);
save('csm_k.mat','tmp1')
