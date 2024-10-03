
#!/bin/bash
cd /rfanfs/pnl-zorro/home/ql087/sms_bart/rawdata/


    kdata_filename="kdata"
    recon_filename="recon_sms"
    bart pics -i50 -e -RW:3:0:1e-06 $kdata_filename sens $recon_filename


