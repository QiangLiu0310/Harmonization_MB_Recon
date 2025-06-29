function [k_trgt, Kimage_short] = recon_sms_std_xa30_6( varargin, prot )
twix = varargin;
if iscell(twix)
    k = twix{length(twix)};
else
    k = twix;
end

NMultiplex = 1;
hdr = parse_measdat_hdr(k.image.filename);
Ncoils=size(prot.sCoilElementID_tElement,2);
NRefLin = hdr.NRefLin;
hdr.NSlc=prot.sSliceArray_lSize;
NSlc = hdr.NSlc; 
NRep = k.image.NRep;
hdr.R=prot.lAccelFactPE; 
NImgLin= length(k.image.Lin(1):hdr.R:k.image.NLin);
NLocPhz=3;
Ngroup=hdr.NSlc/hdr.SMS;
PhaseShift=2*pi/hdr.FOVshift;
slic_indexS=[2:2:NSlc, 1:2:NSlc]; % QL
[~,order_acs]=sort(slic_indexS);
starP=NImgLin*NMultiplex*NSlc+1;
endP =NImgLin*NMultiplex*NSlc+1+NImgLin*NMultiplex*(Ngroup-1);
intervalP=NImgLin*NMultiplex;
slicePos =k.image.slicePos(1:3,starP:intervalP:endP);
[~,order_MB]=sort(sum(slicePos,1));
evp.NSlcMeas=NSlc; 

%% ramp sampling correction
v = extract_vrgf(k.image.filename,prot);      
NFreq_inres=k.phasecor.sqzSize(1);
NFreq_outres=size(v,2);

if hdr.R>2
    NSeg=hdr.R; %reference data to be acquired segmented
else
    NSeg=1;
end

sz_refscanPC=[NFreq_outres  Ncoils  NLocPhz          NMultiplex*NSlc  NSeg]; 
sz_refscan  =[NFreq_outres  Ncoils  NRefLin/NSeg     NMultiplex*NSlc   NSeg];

sz_phasecor =[NFreq_outres  Ncoils  NLocPhz        NMultiplex   Ngroup   NRep];
sz_image =   [NFreq_outres  Ncoils  NImgLin        NMultiplex   Ngroup   NRep];


%% GRAPPA calibration data extraction
refscanPC_raw_tmp  = k.refscanPC.unsorted();
refscanPC_raw = v'*reshape(refscanPC_raw_tmp,[NFreq_inres Ncoils*NLocPhz*NMultiplex*NSlc*NSeg]);
refscanPC_raw=reshape(refscanPC_raw,[NFreq_outres Ncoils  size(refscanPC_raw,2)/Ncoils]);
refscanPC_ROp_ROn=zeros(NFreq_outres,Ncoils,2,NSlc);
refscanPC_ROp_ROn(:,:,1,:)=squeeze((refscanPC_raw(:,:,1:3:end)+refscanPC_raw(:,:,3:3:end))/2);
refscanPC_ROp_ROn(:,:,2,:)=refscanPC_raw(:,:,2:3:end);
%refscan: ACS data that can be used to train GRAPPA coefficients
refscan_raw_tmp = k.refscan.unsorted();
refscan_raw = reshape(v'*reshape(refscan_raw_tmp,[NFreq_inres Ncoils*NRefLin*NMultiplex*NSlc]), sz_refscan);

refscanPC  = squeeze(refscanPC_ROp_ROn(:,:,:,order_acs,:));
refscan    = squeeze(refscan_raw(:,:,:,order_acs,:));
refscan_pc_tmp     = PhaseCorrect_yang(refscanPC,refscan);

for i=1:NSeg
    refscan_pc(:,:,i:NSeg:NRefLin,:)=refscan_pc_tmp(:,:,:,:,i);
end

C2=permute(refscan_pc,[1 3 4 2]);
k_trgt = C2;

start_index = floor((32 - 24) / 2) + 1; 
end_index = start_index + 23;
k_trgt = k_trgt(:, start_index:end_index, :, :);


phzabs = exp(j*(PhaseShift*(0:(size(k_trgt,2)-1))/hdr.R - (0.0*pi) - PhaseShift ) );

% for cnt=0:(hdr.SMS-1),
%     k_trgt(:,:,(cnt*Ngroup)+(1:Ngroup),:) = tmult( k_trgt(:,:,(cnt*Ngroup)+(1:Ngroup),:), diag(phzabs.^cnt), 2);
% end;

%% image data extraction
phasecor_raw_tmp = k.phasecor.unsorted();
phasecor_raw =phasecor_raw_tmp(:,:,NLocPhz*NSlc+1:end);
vrgf_phasecor =reshape(v'*reshape(phasecor_raw,[NFreq_inres Ncoils*NLocPhz*Ngroup*NRep]), sz_phasecor);

% imgscanPC dimension:
vrgf_phasecor_ROp_ROn(:,:,1,:,:,:)= squeeze((vrgf_phasecor(:,:,1,:,:,:)+vrgf_phasecor(:,:,3,:,:,:))/2);
vrgf_phasecor_ROp_ROn(:,:,2,:,:,:)= squeeze( vrgf_phasecor(:,:,2,:,:,:));

%imgscan: ACS data that can be used to train GRAPPA coefficients
image_raw_tmp = k.image.unsorted();
image_raw=image_raw_tmp(:,:,NImgLin*NMultiplex*NSlc+1:end);
vrgf_image=reshape(v'*reshape(image_raw,[NFreq_inres Ncoils*NImgLin*Ngroup*NRep]), sz_image);
Kimage= PhaseCorrect_yang(vrgf_phasecor_ROp_ROn,vrgf_image);
Kimage_short  =squeeze(Kimage(:,:,:,1,order_MB,:));
