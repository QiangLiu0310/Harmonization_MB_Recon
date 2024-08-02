function [Img_PAT, K_PAT, whtcc]= ref_prewht(ref,noise,N,nc)

% gre = squeeze(ref);

gre = permute(squeeze(ref), [1 2 4 3]);
gre = fft2call(gre(1+end/4:3*end/4,:,:,:,:));
% gre = fft2call(gre);   % go to k-space
[n(1), n(2), ~, ~] = size(gre);
% zero pad patref k-space
rref = ifft2call(padarray(gre, (N-n)/2));


% load noise matrix
noise_matrix = squeeze(noise);
nmat = permute(noise_matrix, [2, 1, 3, 4]);
nmat = reshape(nmat, size(nmat, 1), prod(size(nmat))/size(nmat, 1));


[rx, ry, rz, rc] = size(rref);
% SVD to get coil compression matrix.
rref = reshape(rref, rx * ry * rz, rc).';
[u, s, ~] = svd(rref, 'econ');
u = u(:, 1:nc);
s = diag(s);

% Coil compressing noise.
nmat = u' * nmat;

% Estimating whitening matrix.
covm = (nmat * nmat')/(size(nmat, 2) - 1);
whmt = inv(chol(covm, 'lower'));
whmt=whmt./max(max(abs(whmt)));

% Joint coil compression and whitening matrix.
whtcc = whmt * u';

% Whiten and coil compress reference data.
rref = whtcc * rref;
% rref = u'*rref;
rref = rref.';
Img_PAT = reshape(rref, rx, ry, rz, nc);
Img_PAT = permute( Img_PAT, [1 2 4 3]);

K_PAT = fft2call(Img_PAT);



end