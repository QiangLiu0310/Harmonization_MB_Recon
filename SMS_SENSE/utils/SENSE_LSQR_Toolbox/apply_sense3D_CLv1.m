function [ res, tflag ] = apply_sense3D_CLv1( in, params, tflag )
%APPLY_COLLAPSE Summary of this function goes here
%   Detailed explanation goes here

if strcmp(tflag,'transp');
    
    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels
    b = in(1+params.num_chan*prod(params.N):end);
    
    kspace_coils = reshape(in(1:params.num_chan*prod(params.N)), [params.N, params.num_chan]);
    
    img_coils = ifft3c( kspace_coils .* params.m3d );
    
    Res = sum(img_coils .* conj(params.sens), 4);
    
    res = Res(:) + sqrt(params.lambda) * b;
    

else
    
   
        % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT
    
%     img_coils = repmat(reshape(in, params.N), [1,1,1,params.num_chan]) .* params.sens;
    img_coils = bsxfun(@times, reshape(in, params.N), params.sens);

    kspace_coils = fft3c(img_coils) .* params.m3d;
    
    res = cat(1, kspace_coils(:), sqrt(params.lambda)*in);
    

    
end

end


