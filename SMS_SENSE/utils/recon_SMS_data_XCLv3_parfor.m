function [img_recon] = recon_SMS_data_XCLv3_parfor(kspace_cor, ky_idx,sens_gre,AccY, AccZ,PhaseShiftBase,show_mercy)

if (PhaseShiftBase ~=0)
    peshift = ceil(2*pi/PhaseShiftBase);
else
    peshift = 0;
end

[N(1), N(2), num_chan, num_slc]= size(sens_gre);

img_recon=zeros(N(1),N(2),1*AccZ,num_slc/AccZ);
num_slc_raw=num_slc/AccZ;

lsqr_iter = 300;
lsqr_tol = 1e-3; % 1e-3
mask = zeros(N(1),N(2),num_chan);
mask(:,ky_idx,:) =1;


%% determine how many  cores u gonna use
delete(gcp('nocreate'))
c = parcluster('local'); 
total_cores = c.NumWorkers;
be_mercy=show_mercy;  
for ii=1:num_slc_raw
    use_cores=ceil(num_slc_raw/ii);
    if use_cores<=total_cores-be_mercy
        break
    end
end
parpool(use_cores)

parfor ii_slc=1:size(kspace_cor,4) % parfor
    slc_select = ii_slc;
    if AccZ>1
        for ii_temp=2:AccZ
            slc_select=[slc_select,ii_slc+num_slc_raw*(ii_temp-1)];
        end
    end
    img_sense = zeros([N.*[AccZ,1],1]);

    % sens shift
    sens_gre_shift =sens_gre(:,:,:,slc_select);
    sensitivity = zeros(size(sens_gre_shift));
    if (PhaseShiftBase ~=0)
        pes_index  = mod((1:AccZ)-1,2) ;
        peshift=ceil(2*pi/pi); % there is no -
        for jj=1:AccZ
            sensitivity(:,:,:,jj) = circshift(sens_gre_shift(:,:,:,jj),ceil(N(2)/(AccY*peshift))*pes_index(jj),2);
        end

    else
        sensitivity  = sens_gre_shift;
    end

    % AP data
    kspace_ap_collaps = kspace_cor(:,:,:,slc_select(1));
    img_ap_collaps = ifft2call (kspace_ap_collaps );
    img_ap_collaps = repmat(img_ap_collaps, [AccZ,1,1,1,1])/AccZ;
    kspace_ap = fft2call(img_ap_collaps);
    kspace_ap = kspace_ap.* repmat(mask,[AccZ,1,1,1,1]);

    if AccZ>1
        for ii_temp=2:AccZ
            kspace_ap(ii_temp:AccZ:end,:,:,:,:) = 0;
        end
    end

    % SMS-sense using LSQR
    param = [];
    param.N = N .* [AccZ,1];
    param.num_chan = num_chan;
    param.lambda =1e-3; % 1e-3
    param.sens = permute(reshape(permute(sensitivity,[2 3 1 4]),[N(2), num_chan, N(1)*AccZ]),[3 1 2]);
    param.m2d = kspace_ap~=0;
    % sms-sense for AP shot
    res_ap = lsqr(@apply_sense_tikc, cat(1, kspace_ap(:), zeross([prod(N)*AccZ,1])), lsqr_tol, lsqr_iter, [], [], [], param);
    img_sense(:,:,1) = reshape(res_ap, N.*[AccZ,1]);

    img_sense=reshape(permute(reshape(img_sense,[N(1),AccZ,N(2),1]),[1,3,4,2]),[N(1),N(2),1*AccZ]);

    % shift back
    % peshift=ceil(2*pi/ (-pi)); %  (2/3*pi)
peshift=ceil(  (2/3*pi)); %  (2/3*pi)
    pes_index  = mod((1:AccZ)-1,2) ;
    for ii_temp=1:length(pes_index)
        img_sense(:,:,ii_temp) = circshift(img_sense(:,:,ii_temp),ceil(N(2)/(AccY*peshift))*pes_index(ii_temp),2);
    end
    img_recon(:,:,:,ii_slc)=img_sense;

end

img_recon=reshape(permute(reshape(img_recon,[N(1),N(2),1,AccZ,num_slc/AccZ]),[1,2,3,5,4]),[N(1),N(2),1,num_slc]);
img_recon = permute(img_recon, [1 2 4 3]);

end