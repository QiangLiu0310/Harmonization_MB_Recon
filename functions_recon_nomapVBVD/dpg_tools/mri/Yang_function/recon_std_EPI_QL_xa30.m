function [Kimage_first_full] = recon_std_EPI_QL_xa30( varargin , prot)
% twix_  = varargin{1};
% twix=twix_{1,2};

twix = varargin;
if iscell(twix)
    k = twix{length(twix)};
else
    k = twix;
end


twix=k; % debug QL
hdr = parse_measdat_hdr(k.image.filename);

% Ncoils = hdr.Ncoils;
% NSlc = hdr.NSlc;
% R = hdr.R;
% slic_indexS = hdr.slcindx;

R = 2; 
NSlc=40;
Ncoils=32; 
NRefLin = hdr.NRefLin;
NEcho =1;
NPhaCorrLin=3;
NRep = twix.image.NRep;
slic_indexS=[2:2:NSlc, 1:2:NSlc]; % QL
[~,reordr_S]=sort(slic_indexS);

%% % for higher efficiency, high-spatital frequency EPI data is sampled 'on the
% ramps', meaning before the gradient trajectory is stable.  Thus, it
% needs to be 'gridded' onto a constant-velocity-gradient grid prior to
% processing. The VRGF approach does this, e.g. converting 256 sampled-grid
% points into 128 Cartesian-grid points.
%
v = extract_vrgf(twix.image.filename, prot);         % build a VRGF conversion

NFreq_inres=twix.phasecor.sqzSize(1);
NFreq_outres=size(v,2);
NImgLin = size(twix.image.unsorted(),3)/(NEcho*NSlc*NRep);

if R>2
    NSeg=R; %reference data to be acquired segmented
else
    NSeg=1;
end

sz_refscanPC=[NFreq_outres  Ncoils  NPhaCorrLin    NEcho   NSlc   NSeg];
sz_refscan  =[NFreq_outres  Ncoils  NRefLin/NSeg   NEcho   NSlc   NSeg];

sz_phasecor =[NFreq_outres  Ncoils  NPhaCorrLin    NEcho   NSlc   NRep];
sz_image =   [NFreq_outres  Ncoils  NImgLin        NEcho   NSlc   NRep];



%% GRAPPA Reference data extraction 

if R~=1
    
    %refscanPC: ghost-correction data for k.refscan
    refscanPC_raw = tmult(twix.refscanPC.unsorted(),v',1);
    vrgf_refscanPC = reshape(refscanPC_raw, sz_refscanPC);
    % refscanPC dimension:
    % 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','2(first and second Echo)','NSlc
    vrgf_refscanPC_ROp_ROn(:,:,1,:,:,:)= squeeze((vrgf_refscanPC(:,:,1,:,:,:)+vrgf_refscanPC(:,:,3,:,:,:))/2);%ROp
    vrgf_refscanPC_ROp_ROn(:,:,2,:,:,:)= squeeze(vrgf_refscanPC(:,:,2,:,:,:));%ROn
    
    
    %refscan: ACS data that can be used to train GRAPPA coefficients
    vrgf_refscan_raw= tmult(twix.refscan.unsorted(),v',1);
    vrgf_refscan = reshape(vrgf_refscan_raw, sz_refscan);
    
end

%% image data extraction
phasecor_raw   = tmult(twix.phasecor.unsorted(),v',1);
vrgf_phasecor  = reshape(phasecor_raw, sz_phasecor);
% imgscanPC dimension:
% 'NFreq_outres', 'Ncoils', '2(RO+ and RO-)','2(first and second Echo)','NSlc
vrgf_phasecor_ROp_ROn(:,:,1,:,:,:) = squeeze((vrgf_phasecor(:,:,1,:,:,:)+vrgf_phasecor(:,:,3,:,:,:))/2);%ROp
vrgf_phasecor_ROp_ROn(:,:,2,:,:,:) = squeeze(vrgf_phasecor(:,:,2,:,:,:));%ROn

%imgscan: ACS data that can be used to train GRAPPA coefficients
vrgf_image_raw = tmult(twix.image.unsorted(),v',1);
vrgf_image     = reshape(vrgf_image_raw, sz_image);

%%  the next step is to generate the grappa parameters


%% ########-----data order from mp2echo-----########
%%                            ACSscan     
%       refscanPC:   refscanPC#2   refscanPC#1   
%        refscan :     refscan#1     refscan#2    
%-->--->--->-->--->--->--sequential data storage-->--->--->-->--->
%
%%                           Imagescan
%        phasecor:    phasecor#2    phasecor#1    
%           image:       image#1       image#2       
%-->--->--->-->--->--->--sequential data storage-->--->--->-->--->
% The difference between PAT2 and PAT3 could be related to the fact that for PAT3 the PAT reference data are
% usually acquired segmented while for PAT2 this is done in a single shot
%
%% ########---------------------------------########

if R~=1
    vrgf_refscanPC_first = squeeze(vrgf_refscanPC_ROp_ROn(:,:,:,reordr_S,:));
    vrgf_refscan_first = squeeze(vrgf_refscan(:,:,:,1,reordr_S,:));
end

vrgf_phasecor_first = squeeze(vrgf_phasecor_ROp_ROn(:,:,:,reordr_S,:));
vrgf_image_first    = squeeze(vrgf_image(:,:,:,1,reordr_S,:));

if R~=1
    refscan_first_temp = PhaseCorrect_yang(vrgf_refscanPC_first,vrgf_refscan_first);
    for i=1:NSeg
        refscan_first(:,:,i:NSeg:NRefLin,:)=refscan_first_temp(:,:,:,:,i);
    end
end

sz_k =[NFreq_outres  Ncoils  NImgLin     NSlc   NRep];
Kimage_first = reshape(PhaseCorrect_yang(vrgf_phasecor_first,vrgf_image_first),sz_k);



% if R~=1
%     
%     if ~exist('chi','var'); chi = 1e-6; end;
%     if ~exist('eta','var'); eta = 1; end;
%     %  the next step is to generate the grappa parameters
%     
%     if R~=1
%         for slc = 1:NSlc
%             %[~,~,Ng_first{slc}] = recongrappa_multik_Yang(refscan_first(:,:,:,slc), vec(1:NRefLin),'kernel','2x5','dks',R*[1;2]);
%             [~,~,Ng_first{slc}]=recongrappa_multik([hdr.NRefLin size(v,2) hdr.Ncoils],permute(refscan_first(:,:,:,slc),[3 1 2]),[],'kernel','2x5','dks',hdr.R*[1 2],...
%                 'chi',chi,'eta',eta );
%         end
%     end
%     
% end

% pad the data?
Kimage_first_full=zeros(NFreq_outres, Ncoils, NImgLin*R,NSlc,NRep);
Kimage_first_full(:,:,R:R:end,:,:) = Kimage_first;    





    
end

