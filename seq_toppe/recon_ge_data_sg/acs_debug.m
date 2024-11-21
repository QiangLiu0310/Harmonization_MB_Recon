% try Entropy NGC
clear;close all;clc;
load('/data/pnlx/home/ql087/data_bwh/product_scan_rescan/2024_09_14_bwh_ge_sub5/Exam18545/Series5/vrgf_image_new.mat')
tmp=permute(vrgf_image,[1 3 2 4]);

product=zeros(146,146,44,90);
product(:,1:2:end,:,:)=tmp(:,1:1:end,:,1:90);
product(:,2:2:end,:,:)=-tmp(:,1:1:end,:,91:180);

product=permute(product,[2 1 3 4]);
xKy = fftshift(ifft(ifftshift(product,2),[],2),2);

info.NShot=2;
phasepara=zeros(2,size(xKy,4),size(xKy,5));%notice the size: size(phasepara,1)= the number of phase correction order, size: [2 slice 1] -- QL
CorxKy=zeros(size(xKy));
fit_order=2;% 2:linear fitting ;3:2nd order fitting
%% entropy-based 1D linear phase correction

for iSlice = 1:size(xKy,4)
    [CorxKy(:,:,:,iSlice), phasepara(:,iSlice, 1)] = OIEntropyBasedCor(squeeze(xKy(:,:,:,iSlice)),info.NShot,fit_order);
    rep1phaseError=phasepara(:,iSlice,1); % the time point is used for initial searching point.
    [CorxKy(:,:,:,iSlice), phasepara(:,iSlice)] = OIEntropyBasedCor(squeeze(xKy(:,:,:,iSlice)),info.NShot,fit_order,rep1phaseError);% use first rep as start point
    disp(num2str(iSlice))
end

CorKxy=squeeze(fftshift(fft(ifftshift(CorxKy,2),[],2),2));  








