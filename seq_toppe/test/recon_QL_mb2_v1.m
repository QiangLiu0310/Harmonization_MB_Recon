% GE SA file recon
% Qiang Liu
% PA data
% split-GRAPPA
% Jul/11/2024
clear;close all;clc;

addpath(genpath('/data/pnl/home/ql087/arrShow-develop'));
addpath(genpath('/data/pnl/home/ql087/Bruker_2022'));
addpath(genpath('/data/pnl/home/ql087/functions_recon'));
data_path='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series3/';
data_path1='/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series5/';
load('/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_04_09_ge_sub3/Exam16500/Series3/prot.mat')
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

%{
%% load rawdata
tmp=strcat(data_path, 'kspace.mat');
raw=load(tmp);
raw=raw.kspace;
raw=raw(:,:,:,Ngroup+2:Ngroup+46);
% raw(:,:,2:2:end,:)=flipdim(raw(:,:,2:2:end,:),1);
raw=permute(raw,[1 3 2 4]);
ref=raw(:,1:3,:,:);
raw=raw(:,4:end,:,:);
raw(:,2:2:end,:,:)=flipdim(raw(:,2:2:end,:,:),1);
ref(:,2:2:end,:,:)=flipdim(ref(:,2:2:end,:,:),1);
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
vrgf_image=permute(vrgf_image,[1 3 2 4]);
Kimage= PhaseCorrect_yang(vrgf_phasecor_ROp_ROn,vrgf_image);
Kimage=permute(Kimage,[1 3 2 4]);
Kimage=single(Kimage);
save('kspace_mb2.mat','Kimage','-v7.3')
%}

tmp=strcat(data_path, 'kspace_mb2.mat');
raw=load(tmp);
raw=raw.Kimage;
sli_idx=[1:2:45 2:2:45];
[~,or]=sort(sli_idx);
raw=raw(:,:,:,or);
Kimage_short=permute(raw(:,:,:,:),[1 2 4 3]);


%% grappa kernel
tmp=strcat(data_path1, 'csm_k.mat');
coil=load(tmp);
coil=coil.tmp1;
sli_idx=[1:2:NSlc 2:2:NSlc  ];
[~,or]=sort(sli_idx);
coil=coil(end:-1:1,:,:,or); % FE was flipped
sli_idx=[46:90 1:45 ];
[~,or]=sort(sli_idx);
coil=coil(:,:,:,or); 
coil=coil(:,61:84,:,:);
k_trgt = permute(coil,[1 2 4 3]); 
clear or sli_ idx coil raw
sz = size(k_trgt);

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




