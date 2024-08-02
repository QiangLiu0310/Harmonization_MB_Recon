function [img_recon] = recon_3DSense(kspace, ky_idx,kz_idx,sens_map)


[N(1), N(2), N(3), num_chan]= size(sens_map);  % x*y*z*coil


img_recon = zeros(N) ;


lsqr_iter = 300;
lsqr_tol = 1e-3;
mask = zeros(N(1),N(2),N(3), num_chan);

mask(:,ky_idx,kz_idx,:) =1;
kspace = kspace.*mask;

% SMS-sense using LSQR
param = [];
param.N = N ;
param.num_chan = num_chan;
param.lambda = 1e-3;
param.sens = sens_map;
param.m3d = kspace~=0;

% 3D-SENSE
tic
res = lsqr(@apply_sense3D_CLv1, cat(1, kspace(:), zeross([prod(N),1])), lsqr_tol, lsqr_iter, [], [], [], param);
toc

img_recon = reshape(res, N);



end