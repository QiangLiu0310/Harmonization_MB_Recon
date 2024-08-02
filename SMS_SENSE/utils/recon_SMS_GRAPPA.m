function recon_SMS_GRAPPA(EPI_data, EPI_nav,PAT_k,SMS_ref,SMS_nav,sens_gre,AccZ,evp,prot)


kSize_GRAPPA = prot.kSize_GRAPPA ;
lambda_tik = prot.lambda_tik;
num_acs = prot.num_acs;
PhaseShiftBase = prot.PhaseShiftBase;
start_line = evp.NFirstLin;
iReps = prot.iReps ;

AccY = prot.lAccelFactPE;     % PAT factor
num_slice = size(EPI_data,10)*AccZ;  % number of slices

for zz = 1:num_slice
    DThickness(zz) = prot.sSliceArray(zz).dThickness;
end


for iSlice = 1: size(EPI_data,10)
    StartSlice =iSlice;
    SensMap= sens_gre(:,:,:,StartSlice:(num_slice/AccZ):num_slice);
    kspace_acs= PAT_k(:,:,:,StartSlice:(num_slice/AccZ):num_slice);
    
    disp(['Slice Group # ',num2str(StartSlice)]);
    
    
    %  EPI ghost corrections
    lin_fit_EPInavCor = mrir_artifact_ghost_compute_CorrelationMethod(EPI_nav);
    EPI_data_EPInavCor = mrir_artifact_ghost_correct_CorrelationMethod_v3(EPI_data, lin_fit_EPInavCor);
    data_hybrid = mrir_iDFT_freqencode(EPI_data_EPInavCor);
    EPI_data_EPInavCor=[];
    data_hybrid_tp = mrir_regrid_trapezoid(data_hybrid, prot);
    EPI_data_EPInavCor = mrir_fDFT_freqencode(data_hybrid_tp);
    %     clear data_hybrid data_hybrid_tp
    data_hybrid=[];data_hybrid_tp=[];
    
    
    dims_EPI_nav = size(EPI_nav);
    dims_SMS_nav = size(SMS_nav);
    
    dims_EPI_data = size(EPI_data);
    dims_SMS_ref = size(SMS_ref);
    
    
    lin_fit_EPInavCor_reshape = reshape(lin_fit_EPInavCor,2,dims_EPI_nav(3),dims_EPI_nav(7),dims_EPI_nav(10));
    lin_fit_EPInavCor_SMSref_reshape = mean(repmat(lin_fit_EPInavCor_reshape,[1,1,1,AccZ]),3);
    lin_fit_EPInavCor_SMSref= reshape(lin_fit_EPInavCor_SMSref_reshape,2,dims_SMS_nav(3)*dims_SMS_nav(10));
    lin_fit_EPInavCor=[]; lin_fit_EPInavCor_reshape =[];
    lin_fit_EPInavCor_SMSref_reshape=[];
    dims_EPI_nav=[]; dims_SMS_nav=[]; dims_EPI_data=[]; dims_SMS_ref=[];
    
    % averaged ghost correction by using SMS navigators
    SMS_ref_EPInavCor = mrir_artifact_ghost_correct_CorrelationMethod_v3(SMS_ref, lin_fit_EPInavCor_SMSref);
    smsrefscan_hybrid = mrir_iDFT_freqencode(SMS_ref_EPInavCor);
    SMS_ref_EPInavCor = [];
    smsrefscan_hybrid_tp = mrir_regrid_trapezoid(smsrefscan_hybrid, prot);
    SMS_ref_EPInavCor = mrir_fDFT_freqencode(smsrefscan_hybrid_tp);
    smsrefscan_hybrid=[]; smsrefscan_hybrid_tp=[];
    
    % Slice-GRAPPA reconstruction
    % CaipiShift ACS data
    temp = sum(SMS_ref_EPInavCor(:,:,:,:,:,:,1,:,:,:),8);
    K_SG_acs = temp(:,start_line:AccY:end,:,:,:,:,1,:,:,:);
    temp=[];
    K_SG_acs_mslc = K_SG_acs(:,:,:,:,:,:,1,:,:,start_line:(num_slice/AccZ):num_slice);
    SliceGroup = [0:(AccZ-1)] ;
    K_SG_acs_mslc_shft = zeros(size(K_SG_acs_mslc));
    for ss = 1:AccZ
        K_SG_acs_mslc_shft(:,:,:,:,:,:,1,:,:,ss) = CaipirinhaShift_K_v2(K_SG_acs_mslc(:,:,:,:,:,:,1,:,:,ss),SliceGroup(ss),PhaseShiftBase);
    end
    I_SG_acs_mslc_shft= mrir_iDFT_freqencode(mrir_iDFT_phasencode(K_SG_acs_mslc_shft));
    SliceGroup=[]; K_SG_acs_mslc_shft=[];  K_SG_acs_mslc=[]; K_SG_acs=[];
    
    %  CaipiShift EPI data
    K_SG_EPI = sum(EPI_data_EPInavCor(:,:,:,:,:,:,1,:,:,:),8); % 1st TR
    
    SliceSep = sum(DThickness(StartSlice:(StartSlice+(num_slice/AccZ)-1)));
    %SliceSep = (nSlice/AccZ)*DThickness(1);
    K_SG_EPI_deblur = CaipirinhaDeblur_v4(K_SG_EPI, prot, evp, PhaseShiftBase, SliceSep);
    K_SG_EPI_1slc_deblur =  K_SG_EPI_deblur(:,:,:,:,:,:,1,:,:,StartSlice);
    I_SG_EPI_1slc_deblur = mrir_iDFT_freqencode(mrir_iDFT_phasencode(K_SG_EPI_1slc_deblur));
    K_SG_EPI_deblur=[]; K_SG_EPI_1slc_deblur=[]; K_SG_EPI =[];SliceSep=[];
    
    %  slice-GRAPPA recon
    %     disp(['Slice:', num2str(StartSlice)]);
    img_EPI= permute(squeeze(I_SG_EPI_1slc_deblur),[3 2 1]);
    img_acs = permute( squeeze(I_SG_acs_mslc_shft),[3 2 1 4]);
    [Irecon_SG,Irecon_SG_sep,gfactor,weights1_Slice,weights2_Slice] = SliceGRAPPA_v6_3_wo_lambda(img_EPI,img_acs,3,3,PhaseShiftBase);%% STD
    img_EPI=[]; img_acs=[];
    
    %  correct recon data of slice-GRAPPA with the residul of Navigates of SMS ref and Data
    % prepare the kspace data +mask even/odd lines
    num_line = size(EPI_data,1);
    num_column = size(EPI_data,2);
    num_coil = size(EPI_data,3);
    temp = zeros(num_line,num_column,num_coil,1,1,1,1,1,1,AccZ);
    temp(:,:,:,1,1,1,1,1,1,:) = permute(Irecon_SG_sep(:,:,:,:),[3 2 1 4]);
    Krecon_SG_sep = zeros(num_line,num_column,num_coil,1,1,1,1,1,1,AccZ);
    Krecon_SG_sep  = mrir_fDFT_freqencode(mrir_fDFT_phasencode(temp));
    temp=[];
    
    Krecon_SG_sep_mask = zeros(num_line,num_column,num_coil,1,1,1,1,2,1,AccZ);
    Mask_odd = ones(num_line,num_column,num_coil,1,1,1,1,1,1,AccZ);
    Mask_odd(:,2:2:end,:,:) = 0;
    Mask_even = ones(num_line,num_column,num_coil,1,1,1,1,1,1,AccZ);
    Mask_even(:,1:2:end,:,:) = 0;
    Krecon_SG_sep_mask(:,:,:,1,1,1,1,1,1,:) =Krecon_SG_sep .*Mask_odd;
    Krecon_SG_sep_mask(:,:,:,1,1,1,1,2,1,:) =Krecon_SG_sep .*Mask_even;
    Krecon_SG_sep =[]; Mask_odd=[]; Mask_even=[];
    
    %  add EPI_nav back (inversion)
    EPI_nav1slc = EPI_nav(:,:,:,:,:,:,1,:,:,StartSlice);
    lin_fit_EPInavCor1slc = mrir_artifact_ghost_compute_CorrelationMethod(EPI_nav1slc);
    Krecon_SG_sep_EPInavCor_inv= zeros(num_line,num_column,num_coil,1,1,1,1,2,1,AccZ);
    for ss = 1:AccZ
        Krecon_SG_sep_EPInavCor_inv(:,:,:,1,1,1,1,:,1,ss) = mrir_artifact_ghost_correct_CorrelationMethod_v3(Krecon_SG_sep_mask(:,:,:,1,1,1,1,:,1,ss), -lin_fit_EPInavCor1slc );
    end
    
    % correction using SMS_nav
    SMS_nav1slc= SMS_nav(:,:,:,1,1,1,1,:,1,StartSlice:(num_slice/AccZ):num_slice);
    lin_fit_SMSnavCor1slc = mrir_artifact_ghost_compute_CorrelationMethod(SMS_nav1slc);
    temp_corr = mrir_artifact_ghost_correct_CorrelationMethod_v3(Krecon_SG_sep_EPInavCor_inv, lin_fit_SMSnavCor1slc);
    Krecon_SG_sep_SMSnavCor = sum(temp_corr(:,:,:,1,1,1,1,:,1,:),8);
    Irecon_SG_sep_SMSnavCor = squeeze(mrir_iDFT_freqencode(mrir_iDFT_phasencode(Krecon_SG_sep_SMSnavCor)));
    temp_corr=[]; SMS_nav1slc=[]; lin_fit_SMSnavCor1slc=[]; EPI_nav1slc=[]; lin_fit_EPInavCor1slc=[];
    Krecon_SG_sep_EPInavCor_inv=[];
    temp = squeeze(sos(Irecon_SG_sep_SMSnavCor(:,:,:,2),3,2)-sos(permute(Irecon_SG_sep(:,:,:,2),[3 2 1]),3,2));
    
    %     mosaic(squeeze(sos(Irecon_SG_sep_SMSnavCor(:,:,:,2),3,2)-sos(permute(Irecon_SG_sep(:,:,:,2),[3 2 1]),3,2)),1,1,41,['mean of residul = ',num2str(mean_delta_residul)]);
    temp=[]; mean_delta_residul=[];
    
    %get ghost-free SMS_ref data
    
    lin_fit_SMSnavCor= mrir_artifact_ghost_compute_CorrelationMethod(SMS_nav);
    SMS_ref_SMSnavCor = mrir_artifact_ghost_correct_CorrelationMethod_v3(SMS_ref,lin_fit_SMSnavCor);
    smsrefscan_hybrid = mrir_iDFT_freqencode(SMS_ref_SMSnavCor);
    SMS_ref_SMSnavCor = [];
    smsrefscan_free_hybrid = mrir_regrid_trapezoid(smsrefscan_hybrid, prot);
    SMS_ref_SMSnavCor= mrir_fDFT_freqencode(smsrefscan_free_hybrid);
    smsrefscan_hybrid=[]; smsrefscan_free_hybrid =[];lin_fit_SMSnavCor=[];
    
    
    temp = sum(SMS_ref_SMSnavCor(:,:,:,:,:,:,1,:,:,:),8);
    Krecon_SMS_ref_SMSnavCor = temp(:,2:AccY:end,:,:,:,:,1,:,:,:);
    temp=[];
    Krecon_sms_ref_mslc = Krecon_SMS_ref_SMSnavCor(:,:,:,:,:,:,1,:,:,StartSlice:(num_slice/AccZ):num_slice);
    Irecon_sms_ref_mslc = mrir_iDFT_freqencode(mrir_iDFT_phasencode(Krecon_sms_ref_mslc));
    Irecon_sms_ref_mslc = permute( squeeze(Irecon_sms_ref_mslc),[3 2 1 4]);
    
    
    % in-plane GRAPPA reconstruction
    %  in-plane GRAPPA Recon
    PElines = prot.lPhaseEncodingLines;
    Irecon_EPI_Gfree_fs = zeros(num_line,PElines,num_coil,AccZ);
    Irecon_PAT_acs = zeros(num_line,PElines,num_coil,AccZ);
    
    Irecon_EPI_Gfree_fs_comb = zeros(num_line,PElines,AccZ);
    Irecon_PAT_acs_comb = zeros(num_line,PElines,AccZ);
    
    ksampled_EPI=zeros(num_line,PElines,num_coil,AccZ);
    
    Irecon_SMS_ref_fs = zeros(num_line,PElines,num_coil,AccZ);
    Irecon_SMS_ref_fs_comb = zeros(num_line,PElines,AccZ);
    ksampled_EPI(:,start_line:AccY:end,:,:)=squeeze(Krecon_SG_sep_SMSnavCor);
    
    
    ksampled_SMS=zeros(num_line,PElines,num_coil,AccZ);
    ksampled_SMS(:,start_line:AccY:end,:,:)=squeeze(Krecon_sms_ref_mslc);
    
    weights_Inplane=zeros(num_coil, num_coil,PElines,num_line,AccZ );
    for slc = 1:AccZ
        % in-plane GRAPPA for ky undersampled, ghost-free EPI data (Irecon_SG_sep_SMSnavCor)
        
        [Irecon_EPI_Gfree_fs(:,:,:,slc),img_weights] = grappa_gfactor_2d_vc_CL_v2( ksampled_EPI(:,:,:,slc),kspace_acs(:,:,:,slc), 1, AccY, num_acs, kSize_GRAPPA, lambda_tik ,0, zeros(1,num_coil), zeros(1,num_coil),start_line);
        Irecon_EPI_Gfree_fs_comb(:,:,slc)= sum(conj(SensMap(:,:,:,slc)).*(Irecon_EPI_Gfree_fs(:,:,:,slc)),3) ./ (eps + sum(abs(SensMap(:,:,:,slc)).^2,3));
        if (1)
            [Irecon_SMS_ref_fs(:,:,:,slc),~] =grappa_gfactor_2d_vc_CL_v3(ksampled_SMS(:,:,:,slc),kspace_acs(:,:,:,slc), 1, AccY, num_acs, kSize_GRAPPA, lambda_tik,0, zeros(1,num_coil), zeros(1,num_coil) ,start_line);
            Irecon_SMS_ref_fs_comb(:,:,slc) = sum(conj(SensMap(:,:,:,slc)).*(Irecon_SMS_ref_fs(:,:,:,slc)),3) ./ (eps + sum(abs(SensMap(:,:,:,slc)).^2,3));
        end
        Irecon_PAT_acs (:,:,:,slc) = mrir_iDFT_freqencode(mrir_iDFT_phasencode(kspace_acs(:,:,:,slc))) ;
        Irecon_PAT_acs_comb(:,:,slc) = sum(conj(SensMap(:,:,:,slc)).*(Irecon_PAT_acs (:,:,:,slc)),3) ./ (eps + sum(abs(SensMap(:,:,:,slc)).^2,3));
        
    end
    mosaic(Irecon_SMS_ref_fs_comb,1,3,990,'Recon of SMS Ref');%caxis([0 100]);
    mosaic(Irecon_PAT_acs_comb,1,3,991,'Recon of PAT Ref');%caxis([0 100]);
    mosaic(Irecon_EPI_Gfree_fs_comb,1,3,992,['Recon of EPI TR',num2str(iReps),' slice',num2str(iSlice)]);%caxis([0 50]);
    
end