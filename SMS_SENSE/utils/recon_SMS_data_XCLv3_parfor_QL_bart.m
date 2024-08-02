function [kdata,sens] = recon_SMS_data_XCLv3_parfor_QL_bart(kspace_cor, ky_idx,sens_gre,AccY, AccZ,PhaseShiftBase,show_mercy)


[N(1), N(2), num_chan, num_slc]= size(sens_gre);

num_slc_raw=num_slc/AccZ;
mask = zeros(N(1),N(2),num_chan);
mask(:,ky_idx,:) =1;
kdata=zeros(N(1)*2,N(2),num_chan, num_slc_raw);
sens= zeros(N(1)*2,N(2),num_chan, num_slc_raw);
sens_tmp=zeros(N(1)*2,N(2),num_chan);
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
    kdata(:,:,:,ii_slc)=kspace_ap;
    sens_tmp = permute(reshape(permute(sensitivity,[2 3 1 4]),[N(2), num_chan, N(1)*AccZ]),[3 1 2]);
    sens(:,:,:,ii_slc)=sens_tmp;

end
end