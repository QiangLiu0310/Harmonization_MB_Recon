Offline code for Multi-band two diffusion MRI rawdata recon  

# Toolboxes are required: 
1. The Berkeley Advanced Reconstruction Toolbox (BART) toolbox https://mrirecon.github.io/bart/
2. GE orchestra-sdk-2.1-1.matlab 
3. Siemens READ_MEAS_PROT  read protocol from VB- and VD-style "meas.dat"
   
# Recon algorithms and code:

1. L1-ESPIRiT:
   SMS_SENSE/*_bart.m for generating the cfl and hdr files as required by BART toolbox, and recon Bart_scripts.
2. Split-GRAPPA:  
   # Siemens  
   Siemens_Split-grappa/dpg_tools/mri/Yang_function/recon_sms_test_QL_v0_runthis.m for both AP and PA data  
   # GE
   GE_Split_grappa/recon_ge_data_sg/recon_QL_mb2_ap_v5.m for AP data;
   GE_Split_grappa/recon_ge_data_sg/recon_QL_mb2_pa_v4.m for PA data.  

# Open-source dataset:

The PA test data is available at: https://huggingface.co/datasets/QiangLiu0310/Harmonization_MB_Recon  

The whole dataset that supports the findings of our paper is available upon reasonable request from the corresponding author (Qiang Liu qliu30@mgh.harvard.edu).


