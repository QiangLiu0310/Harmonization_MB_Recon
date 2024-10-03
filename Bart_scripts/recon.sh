
#!/bin/bash
cd /rfanfs/pnl-zorro/home/ql087/sms_bart/rawdata1/

num_dwi=62

for dwi_idx in $(seq 1 $num_dwi)
do

    kdata_filename="kdata_${dwi_idx}"
    recon_filename="recon_sms_${dwi_idx}"
    bart pics -i50 -e -RW:3:0:1e-06 $kdata_filename sens $recon_filename

done
