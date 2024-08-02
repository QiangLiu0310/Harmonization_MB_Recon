function [kspace,ky_idx,kz_idx]=detect_kykzLines(K_epi)
 
K_epi= squeeze(K_epi);
% remove readout oversampling
K_epi = ifftc(K_epi,1);
kspace = fftc(K_epi(1+end/4:3*end/4,:,:,:),1);

% detect ky lines
msk = abs(kspace(:,:,1,1)).^(1/3) > 0.001; % Mask to set unacquired points to zero.
ky_idx = find(msk(1,:)~=0);

clear msk

% detect kz lines
msk = abs(kspace(:,ky_idx(1),:,1)).^(1/3) > 0.001; % Mask to set unacquired points to zero.
kz_idx = find(msk(1,:)~=0);
end