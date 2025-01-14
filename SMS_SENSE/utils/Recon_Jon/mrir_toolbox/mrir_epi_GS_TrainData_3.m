function [epi_raw_output, linear_fit_coeff_output, G, epi_img] = mrir_epi_GS_TrainData_old(meas,KernalX,KernalY)

if isfield(meas, 'patrefscan')
    [epi_raw, acs_raw] = mrir_array_GRAPPA_prune(meas.data, meas.patrefscan, meas.evp);
else
    epi_raw = meas.data;
    acs_raw = meas.patrefscan;
end
meas.data = [];
meas.patrefscan = [];

s = size(epi_raw); 
epi_raw = mrir_iDFT_freqencode(epi_raw);

%gridding
epi_raw  = mrir_regrid_trapezoid(epi_raw, meas.prot);
epi_raw_output= mrir_fDFT_freqencode(sum(epi_raw, 8)); %compress partition
clear PhaseUnwindVec

meas.data_phascor1d = mrir_regrid_trapezoid(meas.data_phascor1d, meas.prot);
clear PhaseUnwindVecPh

%ghost
linear_fit_coeff = mrir_artifact_ghost_compute(meas.data_phascor1d);
linear_fit_coeff_output = mean(reshape(linear_fit_coeff,2,[],size(meas.data_phascor1d,7)),3);

epi_raw = mrir_artifact_ghost_correct(epi_raw, linear_fit_coeff);
epi_raw = sum(epi_raw, 8); %compress partition

epi_raw = mrir_fDFT_freqencode(epi_raw);

if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
    epi_raw = mrir_image_slice_deinterleave(epi_raw);
    epi_raw_output = mrir_image_slice_deinterleave(epi_raw_output);
    linear_fit_coeff_output = reshape(mrir_image_slice_deinterleave(reshape(linear_fit_coeff_output,[2,1,s(3),1,1,1,1,1,1,s(10)])),2,[]);
end

if (~isempty(acs_raw)) % if accelerated acquisition
    epi_acs = mrir_iDFT_freqencode(acs_raw);
    s = size(acs_raw); 
    s2 = size(meas.patrefscan_phascor);
    %gridding
    epi_acs = mrir_regrid_trapezoid(epi_acs , meas.prot);

    meas.patrefscan_phascor = mrir_regrid_trapezoid(meas.patrefscan_phascor, meas.prot);
    
    %ghost
    linear_fit_coeff_acs = mrir_artifact_ghost_compute(meas.patrefscan_phascor);
    epi_acs = mrir_artifact_ghost_correct(epi_acs, linear_fit_coeff_acs);
    epi_acs = sum(epi_acs, 8); %compress partition
    
    epi_acs = mrir_fDFT_freqencode(epi_acs);
    
    if ( strcmp(meas.prot.ucMultiSliceMode, 'MSM_INTERLEAVED') ),
        epi_acs = mrir_image_slice_deinterleave(epi_acs);
    end;
    
    [epi_img, G] = mrir_epi_GRAPPA(epi_raw, epi_acs, meas.prot, meas.evp, KernalX, KernalY, 1, 1);
    
else
    epi_img = mrir_image_crop(mrir_conventional_2d(epi_raw));
    G = [];
    
    s_pad = s;
    s_pad(2) = s(2)*4/3 *1/4;
    s_pad(8) = 1;
    epi_rawfull = cat(2,zeros(s_pad),epi_raw);
    epi = mrir_conventional_2d(epi_rawfull);
    epi = mrir_array_combine(epi,0);
    epi = int16(epi); %o.w. wont save properly for some cases
    save_avw(abs(epi),'img_orderswapNoPhase','f',[2 2 2 2]);
    unix(['fslview img_orderswapNoPhase.nii.gz &'])
    
end


