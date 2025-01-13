Offline code for Multi-band 2 DWI rawdata recon  
Toolboxes are required: 
1. The Berkeley Advanced Reconstruction Toolbox (BART) toolbox https://mrirecon.github.io/bart/
2. GE orchestra-sdk-2.1-1.matlab 
3. Siemens READ_MEAS_PROT  read protocol from VB- and VD-style "meas.dat"
   
In the script, two recon algorithms are included:

1. L1-ESPIRiT:
   Harmonization_MB_Recon/SMS_SENSE/, please run the scripts with the name *_bart and generate the cfl and hdr files as required by BART toolbox.
2. Split-GRAPPA:
   For Siemens, please run Harmonization_MB_Recon/Siemens_Split-grappa/dpg_tools/mri/Yang_function
/recon_sms_test_QL_v0_runthis.m for both AP and PA data;
   For GE, please run Harmonization_MB_Recon/GE_Split_grappa/recon_ge_data_sg
/recon_QL_mb2_ap_v5.m for AP data and Harmonization_MB_Recon/GE_Split_grappa/recon_ge_data_sg
/recon_QL_mb2_pa_v4.m for PA data.

The PA test data is available at: https://huggingface.co/datasets/QiangLiu0310/Harmonization_MB_Recon  

The whole dataset that supports the findings of our paper is available upon reasonable request from the corresponding author (Qiang Liu qliu30@mgh.harvard.edu).


