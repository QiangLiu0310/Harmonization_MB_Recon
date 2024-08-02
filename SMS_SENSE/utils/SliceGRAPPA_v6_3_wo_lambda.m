function [reconComb , reconSep, gfact, weights1_deblur, weights2_deblur] = SliceGRAPPA_v6_3_wo_lambda(imfold,imtrg,nKy,nKx,PhaseShiftBase)
%function [reconComb , reconSep, weights1, weights2] = SliceGRAPPA_v6_3(imfold,imtrg,nKy,nKx,PhaseShiftBase,lamda)
%% ________________________________________________________________________
%% Split-Slice GRAPPA
%% ________________________________________________________________________
%%
%%
%% INPUT: 
%% 
%% imfold               folded image (Nc x Ny x Nx)
%% imtrg                ull image as target data (Nc x Ny x Nx x MB)
%% nKy                  kernel size along ky  (odd number)
%% nKx                  kernel size along kx  (odd number)
%% PhaseShiftBase       CaipiShift Phase
%% lamda                weight regulatzation parameter
%%
%%
%% Outputs:
%%
%% reconSep             the slice-seperated recon dataset in image space
%% reconComb            the sos of recon dataset in image space
%% weights1_deblur      the deblured weights of odd lines in image space 
%% weights2_deblur      the deblured weights of even lines in image space 
%%
%% add leakage block (sp-SG) kernel calc to reduce leakage artifact
%% Sept 13 2016 by Haifeng Wang
%%
%% try to use none deblur weights
%% Sept 14 2016 by Haifeng Wang
%%
%% try to use deblur weights
%% Sept 15 2016 by Haifeng Wang
%%
%% try to seperatly deblur using different CAIPISHIFT for odd and even
%% kernels ?????
%% Sept 24 2016 by Haifeng Wang
%% ________________________________________________________________________

[nc,    nyr,   nx   ]   = size(imfold);
[nctrg, nytrg, nxtrg, MB]   = size(imtrg);
ncsrc = nctrg;
nysrc = nytrg;
nxsrc = nxtrg;
ny = nyr;

if nargin < 7;
    NoiseCov = eye(nc);
%     disp('No noise correlations has been passed ......')
%     disp('Assuming perfect coils without correlations......')
end

% Fourier transform in k-space
trg = ifftshift(ifftshift(fft(fft(fftshift(fftshift(imtrg,2),3),[],2),[],3),2),3);

% Two kernels
[kernel1, kernel2] = calcKernel    (trg,nKy,nKx);

SliceGroup = [0:MB-1];
for slcInd = 1:MB
    % Two weights
    weights1(:,:,:,:,slcInd) = getWeights    (kernel1{slcInd},ny,nx);
    weights2(:,:,:,:,slcInd) = getWeights    (kernel2{slcInd},ny,nx);
    % Two recon
    recon1(:,:,:,slcInd)   = applyWeights  (imfold,weights1(:,:,:,:,slcInd));
    recon2(:,:,:,slcInd)   = applyWeights  (imfold,weights2(:,:,:,:,slcInd));
    % Combine two recon
    k_recon1 = mrir_fDFT_freqencode(mrir_fDFT_phasencode(recon1(:,:,:,slcInd)));
    k_recon2 = mrir_fDFT_freqencode(mrir_fDFT_phasencode(recon2(:,:,:,slcInd)));    
    k_recon(:,1:2:ny,:) = k_recon1(:,1:2:ny,:);
    k_recon(:,2:2:ny,:) = k_recon2(:,2:2:ny,:);
    recon(:,:,:,slcInd) = mrir_iDFT_freqencode(mrir_iDFT_phasencode(k_recon(:,:,:)));
    
    %recon_deblur(:,:,:,slcInd) = mrir_iDFT_phasencode(CaipiShift_K(mrir_fDFT_phasencode(recon(:,:,:,slcInd)),SliceGroup(slcInd),-PhaseShiftBase));
    %imtrg_deblur(:,:,:,slcInd) = mrir_iDFT_phasencode(CaipiShift_K(mrir_fDFT_phasencode(imtrg(:,:,:,slcInd)),SliceGroup(slcInd),-PhaseShiftBase));
    %weights1_deblur(:,:,:,:,slcInd) = permute(mrir_iDFT_phasencode(CaipiShift_K(mrir_fDFT_phasencode(permute(weights1(:,:,:,:,slcInd),[1 3 2 4 5])),SliceGroup(slcInd),-PhaseShiftBase)),[1 3 2 4 5]);%[nc,nc,srcy,srcx,MB]
    %weights2_deblur(:,:,:,:,slcInd) = permute(mrir_iDFT_phasencode(CaipiShift_K(mrir_fDFT_phasencode(permute(weights2(:,:,:,:,slcInd),[1 3 2 4 5])),SliceGroup(slcInd),-PhaseShiftBase)),[1 3 2 4 5]);%[nc,nc,srcy,srcx,MB]
    
    recon_deblur(:,:,:,slcInd) = mrir_iDFT_phasencode(CaipirinhaShift_K_v2(mrir_fDFT_phasencode(recon(:,:,:,slcInd)),SliceGroup(slcInd),-PhaseShiftBase));
    imtrg_deblur(:,:,:,slcInd) = mrir_iDFT_phasencode(CaipirinhaShift_K_v2(mrir_fDFT_phasencode(imtrg(:,:,:,slcInd)),SliceGroup(slcInd),-PhaseShiftBase));
    weights1_deblur(:,:,:,:,slcInd) = permute(mrir_iDFT_phasencode(CaipirinhaShift_K_v2(mrir_fDFT_phasencode(permute(weights1(:,:,:,:,slcInd),[1 3 2 4 5])),SliceGroup(slcInd),-PhaseShiftBase)),[1 3 2 4 5]);%[nc,nc,srcy,srcx,MB]
    weights2_deblur(:,:,:,:,slcInd) = permute(mrir_iDFT_phasencode(CaipirinhaShift_K_v2(mrir_fDFT_phasencode(permute(weights2(:,:,:,:,slcInd),[1 3 2 4 5])),SliceGroup(slcInd),-PhaseShiftBase)),[1 3 2 4 5]);%[nc,nc,srcy,srcx,MB]
       
    
    % Shift 1 for even kernel in k-space ????
    %weights2_deblur(:,:,:,:,slcInd) = permute(mrir_iDFT_phasencode(CaipiShift_K2(mrir_fDFT_phasencode(permute(weights2(:,:,:,:,slcInd),[1 3 2 4 5])),SliceGroup(slcInd),-PhaseShiftBase)),[1 3 2 4 5]);%[nc,nc,srcy,srcx,MB]
    %weights1_deblur(:,:,:,:,slcInd) = permute(mrir_iDFT_phasencode(CaipiShift_K(mrir_fDFT_phasencode(permute(weights1(:,:,:,:,slcInd),[1 3 2 4 5])),SliceGroup(slcInd),PhaseShiftBase)),[1 3 2 4 5]);%[nc,nc,srcy,srcx,MB]
    %weights2_deblur(:,:,:,:,slcInd) = permute(mrir_iDFT_phasencode(CaipiShift_K(mrir_fDFT_phasencode(permute(weights2(:,:,:,:,slcInd),[1 3 2 4 5])),SliceGroup(slcInd),PhaseShiftBase)),[1 3 2 4 5]);%[nc,nc,srcy,srcx,MB]
    
    % test g-factor calculation
    %gfact1(:,:,slcInd)   = calcGfact     (weights1_deblur(:,:,:,:,slcInd),imtrg_deblur(:,:,:,slcInd),NoiseCov);
    %gfact2(:,:,slcInd)   = calcGfact     (weights2_deblur(:,:,:,:,slcInd),imtrg_deblur(:,:,:,slcInd),NoiseCov);
    %gfact(:,:,slcInd)   = calcGfact2kernel(weights1_deblur(:,:,:,:,slcInd),weights2_deblur(:,:,:,:,slcInd),imtrg_deblur(:,:,:,slcInd),NoiseCov);
    gfact(:,:,slcInd)   = 0;
    %figure;imagesc(1./gfact1(:,:,slcInd),[0 1]);title('Slice GRAPPA(odd)');colorbar;
    %figure;imagesc(1./gfact2(:,:,slcInd),[0 1]);title('Slice GRAPPA(even)');colorbar;
    %figure;imagesc(1./gfact(:,:,slcInd),[0 1]);title('Slice GRAPPA');colorbar;
    
end

reconSep = squeeze(recon_deblur);

reconComb = squeeze(sqrt(sum(abs(recon_deblur).^2,1)));

end


%%  Helper functions



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Calculation of GRAPPA weights in kspace from ACS data
%  see Griswold et al
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ws_kernel1,ws_kernel2] = calcKernel(acs_trg,srcy,srcx)

[nc_src,ny_src,nx_src,MB]=size(acs_trg);

nc    = nc_src;
nyacs = ny_src;
nxacs = nx_src;

% tic
centerKyKernelMinus1 = floor(srcy/2);
centerKxKernelMinus1 = floor(srcx/2);
WyMinus1 = (srcy-1);
WxMinus1 = (srcx-1);

src1=zeros(nc*srcy*srcx,ceil(0.5*(nyacs-WyMinus1))*(nxacs-WxMinus1)*MB);
src2=zeros(nc*srcy*srcx,floor(0.5*(nyacs-WyMinus1))*(nxacs-WxMinus1)*MB);
cnt1 = 0;
cnt2 = 0;
clear trg
for slcInd = 1:MB
    slcInd  ;   
    trg1{slcInd} = reshape(acs_trg(:,[1:2:nyacs-WyMinus1]+centerKyKernelMinus1,[1:nxacs-WxMinus1]+centerKxKernelMinus1,slcInd),nc,[]);
    trg2{slcInd} = reshape(acs_trg(:,[2:2:nyacs-WyMinus1]+centerKyKernelMinus1,[1:nxacs-WxMinus1]+centerKxKernelMinus1,slcInd),nc,[]);
    for xind=1:nxacs-WxMinus1
        % odd line
        for yind=1:2:nyacs-WyMinus1
            cnt1=cnt1+1;
            src1(:,cnt1)=reshape(acs_trg(:,yind:yind+WyMinus1,xind:xind+WxMinus1,slcInd),[],1);
        end
        % even line
        for yind=2:2:nyacs-WyMinus1
            cnt2=cnt2+1;
            src2(:,cnt2)=reshape(acs_trg(:,yind:yind+WyMinus1,xind:xind+WxMinus1,slcInd),[],1);
        end
    end
end
% toc

src1_inv = pinv(src1);
src2_inv = pinv(src2);
% src1_inv = pinv_reg(src1,lamda);
% src2_inv = pinv_reg(src2,lamda);
% src1_inv = pseudoinverse(src1,[],'lsqr');
% src2_inv = pseudoinverse(src2,[],'lsqr');

sT1= size(trg1{1});
sT2= size(trg2{1});
for slcInd = 1:MB
    current_trg1 = cat(2,zeros(sT1(1),sT1(2)*(slcInd-1)),trg1{slcInd}, zeros(sT1(1),sT1(2)*(MB -slcInd)));
    current_trg2 = cat(2,zeros(sT2(1),sT2(2)*(slcInd-1)),trg2{slcInd}, zeros(sT2(1),sT2(2)*(MB -slcInd)));
    % odd line
    ws_tmp1 = current_trg1*src1_inv;
    ws_tmp1 = reshape(ws_tmp1,[nc,nc,srcy,srcx]);                                 %Reshape weight set
    ws_kernel1{slcInd} = flipdim(flipdim(ws_tmp1,3),4);                           %flip source points in ky and kx for the convolution
    % even line
    ws_tmp2 = current_trg2*src2_inv;
    ws_tmp2 = reshape(ws_tmp2,[nc,nc,srcy,srcx]);                                 %Reshape weight set
    ws_kernel2{slcInd} = flipdim(flipdim(ws_tmp2,3),4);                           %flip source points in ky and kx for the convolution    
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                   
%                                                                    
%  Calculate grappa weights for reconstruction in image space                                                                  
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ws_img] = getWeights(ws_kernel,ny,nx)

[nc,temp,nky,nkx] = size(ws_kernel);

ws_k = zeros(nc,nc,ny,nx);
ws_k(:,:,ceil((ny-nky)/2)+1:ceil((ny+nky)/2),ceil((nx-nkx)/2+1):ceil((nx+nkx)/2)) = ws_kernel;  %put reconstruction kernel in the center of matrix

tmp0 = ifftshift(ws_k,3);           % shift in phase
tmp1 = ifftshift(tmp0,4);           % shift in read
tmp0 = ifft(tmp1,[],3);             % ifft in phase
tmp1 = ifft(tmp0,[],4);             % ifft in read
tmp0 = ifftshift(tmp1,3);           % shift in phase
tmp1 = ifftshift(tmp0,4);           % shift in read
ws_img = ny*nx*tmp1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Weights application in image space
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [recon] = applyWeights(img_red,ws_img)


nc = size(ws_img,1);
ny = size(ws_img,3);
nx = size(ws_img,4);

nyr = size(img_red,2);

if (ny > nyr);
  
    af = round(ny/nyr);
    
    disp('Assuming the data is passed without zeros at not sampled lines ......')
    disp(['Acceleration factor is af = ' num2str(af) '.....'] );
       
    sig_red = fftshift(fftshift(fft(fft(fftshift(fftshift(img_red,2),3),[],2),[],3),2),3);

    sig_new = zeros(nc,ny,nx);
    sig_new(:,1:af:end,:) = sig_red;
    sig_red = sig_new;
    
    img_red = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sig_red,2),3),[],2),[],3),2),3);
        
end

recon = zeros(nc,ny,nx);

for k = 1:nc,
    recon(k,:,:) = sum(squeeze(ws_img(k,:,:,:)).*img_red,1);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gfactor calculation according to Breuer et al
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [g] = calcGfact(ws_img,imgref,R)


nc = size(ws_img,1);
ny = size(ws_img,3);
nx = size(ws_img,4);

ws_img = ws_img;

if nargin < 3;
    R = eye(nc);
%     disp('No noise correlations has been passed ......')
%     disp('Assuming perfect coils without correlations......')
end

disp('Calculating G-factor......')

sigref = fftshift(fftshift(fft(fft(fftshift(fftshift(imgref,2),3),[],2),[],3),2),3);

[nc, nyref, nxref] = size(sigref);

% % filter kspace
% sigref = bsxfun(@times,sigref,reshape(tukeywin(nyref,1),[1 nyref]));
% sigref = bsxfun(@times,sigref,reshape(tukeywin(nxref,1),[1 1 nxref]));
% % 
% sigreff = zeros(nc,ny,nx);
% yidx = floor((ny-nyref)/2) + [1:nyref];
% xidx = floor((nx-nxref)/2) + [1:nxref];
% sigreff(:,yidx,xidx) = sigref;
% 
% imgref = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sigreff,2),3),[],2),[],3),2),3);

g = zeros(ny,nx);


Z = eye(nc);

for y = 1:ny,
    for x = 1:nx,
        W = ws_img(:,:,y,x);            % Weights in image space
        tmp = imgref(:,y,x);
        n = tmp'./sqrt(sum(abs(tmp).^2,1));
        g(y,x) = sqrt(abs((n*W)*R*(n*W)'))./sqrt(abs((n*Z)*R*(n*Z)'));       % This is the generalized g-factor formulation
    end
    
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gfactor calculation of two kernels
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g] = calcGfact2kernel(ws_img1,ws_img2,imgref,R)


nc = size(ws_img1,1);
ny = size(ws_img1,3);
nx = size(ws_img1,4);

ws_img_Slice_new_odd = ws_img1;
ws_img_Slice_new_even = ws_img2;

ws_img1a = ws_img_Slice_new_odd;
ws_img1b = ws_img_Slice_new_even;
% shift FOV/2 ALONG y
ws_img1a2 = ifftshift(ws_img_Slice_new_odd,3);
ws_img1b2 = ifftshift(ws_img_Slice_new_even,3);

if nargin < 3;
    R = eye(nc);
%     disp('No noise correlations has been passed ......')
%     disp('Assuming perfect coils without correlations......')
end

disp('Calculating G-factor......')

sigref = fftshift(fftshift(fft(fft(fftshift(fftshift(imgref,2),3),[],2),[],3),2),3);

[nc, nyref, nxref] = size(sigref);

% % filter kspace
% sigref = bsxfun(@times,sigref,reshape(tukeywin(nyref,1),[1 nyref]));
% sigref = bsxfun(@times,sigref,reshape(tukeywin(nxref,1),[1 1 nxref]));
% % 
% sigreff = zeros(nc,ny,nx);
% yidx = floor((ny-nyref)/2) + [1:nyref];
% xidx = floor((nx-nxref)/2) + [1:nxref];
% sigreff(:,yidx,xidx) = sigref;
% 
% imgref = ifftshift(ifftshift(ifft(ifft(ifftshift(ifftshift(sigreff,2),3),[],2),[],3),2),3);

g = zeros(ny,nx);


Z = eye(nc);
RR = blkdiag(R,R);
for y = 1:ny,
    for x = 1:nx,
        % Calculate the average slice weight matrices between weights of odd and even lines
        W1a = ws_img1a(:,:,y,x);            % Weights in image space
        W1b = ws_img1b(:,:,y,x);            % Weights in image space
        W11 = (W1a+W1b)./2;
        % Calculate the difference slice weight matrices between weights of odd and even lines
        W1a2 = ws_img1a2(:,:,y,x);            % Weights in image space
        W1b2 = ws_img1b2(:,:,y,x);            % Weights in image space        
        W12 = (W1a2-W1b2)./2;
        % Augment slice weight matrices
        WW1 = [W11,W12];% Weights in image space
        tmp = imgref(:,y,x);
        n = tmp'./sqrt(sum(abs(tmp).^2,1));
        g(y,x) = sqrt(abs((n*(WW1))*RR*(n*(WW1))'))./(sqrt(abs((n*Z)*R*(n*Z)'))+eps);
    end
    
end


end



function X = pinv_reg(A,lambda)
%PINV   Pseudoinverse.
%   X = PINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS.
%
%   PINV(A,TOL) uses the tolerance TOL instead of the default.
%
%   See also RANK.

%   Copyright 1984-2001 The MathWorks, Inc. 
%   $Revision: 5.11 $  $Date: 2001/04/15 12:01:37 $
[m,n] = size(A);

if n > m
   X = pinv_reg(A',lambda)';
else
    
   AA = A'*A;

   S = svd(AA,0);
   %S = svd(AA,'econ');
   
   S = sqrt(max(abs(S)));
   X = (AA+eye(size(AA)).*lambda.*S.^2)\A'; 
  
end

end
