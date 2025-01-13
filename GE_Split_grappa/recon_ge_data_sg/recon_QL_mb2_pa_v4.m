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

data_path='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series3/'; % image
data_path1='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series5/'; % ref

currentDir = pwd;

if ~strcmp(currentDir, data_path1)
    cd(data_path1);
end
disp('Currently in directory ref scan');

%% read GE 2-shot ref scan
% image parameters
NFreq_outres=146;
NFreq_inres=330;
Ncoils=44;
NLocPhz=3;
NImgLin=73;
NSlc=90;
SMS=2;
Ngroup=NSlc*2;
NRep=1;

%% load rawdata

tmp=strcat(data_path1,'ScanArchive_LONGWOOD30MR2_20240409_212133980.h5'); % the last one
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

kspace = complex(zeros(xRes,  nChannels, yRes, zRes));

for i = 1:archive.ControlCount

    control = GERecon('Archive.Next', archive);

    if isfield(control, 'Data')
        kspace(:,:,1:2:end,i) = squeeze(control.Data);
    else
        fprintf('%d\n', i);
    end

end
GERecon('Archive.Close', archive);
kspace=single(kspace);

raw=kspace(:,:,1:2:end,:); 
clear kspace
raw(:,:,2:2:end,:)=flipdim(raw(:,:,2:2:end,:),1);
raw=permute(raw,[1 3 2 4]);

tmp=strcat(data_path1,'ScanArchive_LONGWOOD30MR2_20240409_212102346.h5'); % the first one
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

kspace = complex(zeros(xRes,  nChannels, yRes, zRes));

for i = 1:archive.ControlCount

    control = GERecon('Archive.Next', archive);

    if isfield(control, 'Data')
        kspace(:,:,1:2:end,i) = squeeze(control.Data);
    else
        fprintf('%d\n', i);
    end

end
GERecon('Archive.Close', archive);
kspace=single(kspace);
kspace=kspace(:,:,:,1:2:end);
ref=kspace;
clear kspace
ref(:,:,3:4:end,:)=flipdim(ref(:,:,3:4:end,:),1);
ref=permute(ref,[1 3 2 4]);
ref=ref(:,73:2:77,:,:);

%% load vrgf dat0
tmp=strcat(data_path1, 'vrgf.dat');
fid = fopen(tmp, 'r');
data = fread(fid, inf, 'float'); % 'float' specifies the data type
fclose(fid);
v=reshape(data,[330,146]);
clear data

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
clear tmp
tmp1=single(tmp1);

coil=tmp1; clear tmp1
sli_idx=[1:2:NSlc 2:2:NSlc];
[~,or]=sort(sli_idx);
coil=coil(end:-1:1,:,:,or); % FE was flipped
% sli_idx=[46:90 1:45 ];
% [~,or]=sort(sli_idx);
% coil=coil(:,:,:,or); 
coil=coil(:,61:84,:,:);
k_trgt = permute(coil,[1 2 4 3]); 

clearvars -except k_trgt data_path data_path1
sz = size(k_trgt);

%% read GE MB2 data

cd(data_path);
disp('Currently in directory rawdata');

tmp=strcat(data_path,'ScanArchive_LONGWOOD30MR2_20240409_211609044.h5');
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
NRep=1;
load('/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series3/prot.mat')

%% load rawdata
raw=kspace(:,:,:,Ngroup+2:Ngroup+46); 
ref=kspace(:,:,4:6,1:Ngroup);
ref=permute(ref,[1 3 2 4]);
% clear kspace
raw=permute(raw,[1 3 2 4]);
raw=raw(:,4:end,:,:);
raw(:,2:2:end,:,:)=flipdim(raw(:,2:2:end,:,:),1);
ref(:,2:2:end,:,:)=flipdim(ref(:,2:2:end,:,:),1);
clear tmp

%% load vrgf dat0
tmp=strcat(data_path1, 'vrgf.dat');
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
vrgf_image=permute(vrgf_image,[1 3 2 4]);
Kimage= PhaseCorrect_yang(vrgf_phasecor_ROp_ROn,vrgf_image);
Kimage=permute(Kimage,[1 3 2 4]);
Kimage=single(Kimage);

sli_idx=[1:2:45 2:2:45];
[~,or]=sort(sli_idx);
Kimage=Kimage(:,:,:,or);
Kimage_short=permute(Kimage(:,:,:,:),[1 2 4 3]);
clear Kimage


%% Grappa kernel

if ~exist('slices','var'); slices = 1:hdr.NSlc; end;
if ~exist('chi','var'); chi = 1e-6; end;
if ~exist('eta','var'); eta = 1; end;

phzabs = exp(j*(PhaseShift*(0:(size(k_trgt,2)-1))/hdr.R - (0.0*pi) - PhaseShift ) );

for cnt=0:(hdr.SMS-1),
    k_trgt(:,:,(cnt*Ngroup)+(1:Ngroup),:) = tmult( k_trgt(:,:,(cnt*Ngroup)+(1:Ngroup),:), diag(phzabs.^cnt), 2);
end;

if ~exist('Npppp','var')
    
    fprintf('%d',mod(1:length(slices),10));fprintf('\n');
    for slc = 1:NSlc;
        fprintf('o');
        [~,~,Np{slc}]=recongrappa_multik(sz([2 1 4]),permute(k_trgt(:,:,slc,:),[2 1 4 3]),[],'kernel','2x5','dks',hdr.R*[1 2],...
            'chi',chi,'eta',eta );
    end;
    fprintf('\n');
    
end

z = [];
for cnt=1:length(Np{slices(1)})
    dk = diff(find( Np{slices(1)}(cnt).pattern == '*' ));
    z(dk) =  1;
end;

if ( (length(z)< 2*hdr.R) || ( z( 2*hdr.R ) == 0 ) )
    disp(['for this script to run properly, GRAPPA parameters in Np need to cover both R and 2R accelerations']);
    error(['need to regenerate Np parameters'])
end;

%% recon
prot.sSliceAcceleration.lMultiBandFactor=hdr.SMS;
prot.sSliceAcceleration.lFOVShiftFactor=hdr.FOVshift;

for slc = 1:Ngroup
    iii = 1:size(k_trgt,2);
    sz = size(k_trgt);
    for cnt=1:hdr.SMS,
        in2{cnt} = zeros( sz([1 2 4]) );
        in2{cnt}(:,(sz(2)-length(iii))/2+iii,:) = squeeze( k_trgt(:,:,slc+(cnt-1)*Ngroup,:) );
    end;
    [~,w{slc}] =  MultisliceGRAPPA_2kernal_leakBlock( in2{1}, in2, [ 5 3 1 hdr.R ], 'full', prot);
    
end

for cRep=1:NRep
    
    k_data_gc0 = Kimage_short;
    k_data_gc  = k_data_gc0(:,:,:,:,cRep);
    %     k_data_gc_deblur = CaipirinhaDeblur_v3_ge_1( k_data_gc(:,:,:,:) );
    k_data_gc_deblur=k_data_gc;
    if exist('runmod')
        if ( runmod == 1 )
            k_data_gc_deblur = k_data_gc_deblur(:,:,:,coils);
        end;
    end;
    
    nky = size(k_data_gc_deblur,2);
    
    fprintf('%d',mod(1:length(slices),10)); fprintf('\n');
    for slc = 1:Ngroup;
        fprintf('.');
        
        in1 = zeros([ size(k_data_gc_deblur,1) nky*hdr.R sz(4)]);
        
        in1(:,hdr.R*(0:size(k_data_gc_deblur,2)-1)+1,:) = squeeze( k_data_gc_deblur(:,:,slc,:) );
        
        tmp = squeeze(MultisliceGRAPPA_2kernal_leakBlock( in1, w{slc}, [ 5 3 1 hdr.R ]) );
        
        phzabs = exp(j*(PhaseShift*(0:(size(tmp,2)-1))/hdr.R - (0.0*pi) - PhaseShift ) );
        slcgrp = slc + [ 0:(hdr.SMS-1) ]*Ngroup;
        
        sz_in1 = size(in1);
        for cnt=1:hdr.SMS,
            curslc = slcgrp(cnt);
            Fa1 = recongrappa_multik(sz_in1([2 1 3]),permute(tmp(:,1:2*hdr.R:end,:,cnt),[2 1 3]),1:2*hdr.R:sz_in1(2),'kernel','2x5','N',Np{curslc});
            Fb1 = recongrappa_multik(sz_in1([2 1 3]),permute(tmp(:,(1+hdr.R):2*hdr.R:end,:,cnt),[2 1 3]),(1+hdr.R):2*hdr.R:sz_in1(2),'kernel','2x5','N',Np{curslc});
            
            Fa2 = permute(tunfold(fif(Fa1),2),[2 1]);
            Fb2 = permute(tunfold(fif(Fb1),2),[2 1]);
            
            Fb3 = phzshift( Fa2, Fb2,{'nofft','nocombo'} );
            Fb4 = ifi( trefold(permute(Fb3,[2 1]),sz_in1([2 1 3]),2) );
            
            Fc1 = zeros(size(Fa1));
            Fc1( 1:2*hdr.R:end, :, : ) = Fa1( 1:2*hdr.R:end, :, : );
            Fc1( (1+hdr.R):2*hdr.R:end, :, : ) = Fb4( (1+hdr.R):2*hdr.R:end, :, : );
            
            Fd1 = recongrappa_multik(size(Fc1),Fc1,[],'kernel','2x5','dks',hdr.R,'N',Np{curslc});
            
            F2klb(:,:,curslc,:) = Fd1;
            F2klb(:,:,curslc,:) =  tmult( F2klb(:,:,curslc,:), diag(conj(phzabs).^(cnt-1)), 1);
        end;
        
        %if verbose, keyboard; end;
    end
    fprintf('\n');
    k_all(:,:,:,:,cRep)=F2klb;
end

%% partial Fourier
aa=k_all;
for cntsli=1:size(aa,3)
    
    for cntcoil=1:size(aa,4)
        for cRep=1:size(aa,5)
            img_aa(:,:,cntcoil,cRep)=flip(reconhd(aa(:,:,cntsli,cntcoil,cRep),size(aa,1),size(aa,2)));
        end
    end
     fprintf('Loop counter: %d\n', cntsli);
    img(:,:,cntsli,:)=sqrt(sum(abs( fif( img_aa ) ).^2,3));
end

sli_idx=[46:90 1:45];
[~,or]=sort(sli_idx);
img=img(:,:,or); 


cd /data/pnlx/home/ql087/data_processing/2024_recon_reproducibility/sub3/ge/scan1/split_grappa/
save('pa_scan1_sg.mat','img')

