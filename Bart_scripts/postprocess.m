
%% bart
% bart pics -i100 -e -RW:3:0:1e-06 kdata sens recon_sms
% read the data from bart
clear;close all;clc;
N=[146 146];
AccZ=2;
AccY=2;
num_slc=90;
ndwi=62;
img_bart=zeros(146,146,num_slc,ndwi);

for i_dwi=1:ndwi
% res_ap=squeeze(readcfl('/rfanfs/pnl-zorro/home/ql087/sms_bart/rawdata/recon_sms_1'));

filename = sprintf('/rfanfs/pnl-zorro/home/ql087/sms_bart/rawdata1/recon_sms_%d', i_dwi);
% filename = sprintf('/rfanfs/pnl-zorro/home/ql087/sms_bart/rawdata/recon_sms');
res_ap = squeeze(readcfl(filename));
img_sense = zeros([N.*[AccZ,1],num_slc/2]);
mtx=[N.*[AccZ,1],num_slc/2];
img_sense= reshape(res_ap, mtx);
img_sense=reshape(permute(reshape(img_sense,[N(1),AccZ,N(2),num_slc/2]),[1,3,4,2]),[N(1),N(2),num_slc]);

% shift back
% peshift=ceil(2*pi/ (2/3*pi)); % siemens data
peshift=ceil(2*pi/ (-pi)); % ge data

for ii_temp=num_slc/2+1:num_slc
    img_sense(:,:,ii_temp) = circshift(img_sense(:,:,ii_temp),ceil(N(2)/(AccY*peshift)),2);
end

img_bart(:,:,:,i_dwi)=img_sense;
clear img_sense
end


