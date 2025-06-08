# Multi-band Reconstruction for Recon Reproducibility Study

This repository contains code and workflows for evaluating multi-band (MB) reconstruction methods in the context of diffusion MRI reconstruction reproducibility.

## 🛠 Recon Methods

1. **SMS_SENSE**  
2. **SMS_SENSE with BART**  
3. **Split Slice-GRAPPA**

## 🧾 Reconstruction Steps

Split Slice-GRAPPA
1. GRAPPA ACS data extraction, EPI ramp sampling regrid, EPI Nyquist ghost correction: odd-even readout lines phase correction with linear phase correction (LPC), GRAPPA kernel training
2. image data extraction, EPI ramp sampling regrid, Nyquist ghost correction with LPC, GRAPPA recon
3. Partial Fourier recon with homodyne algorithm (reconhd)
4. Coil combine with sum of square along the coil dimension

L1-ESPIRiT
1. extract the reference data for estimating the coil sensitivity maps (Siemens: EPI, GE: low resolution gre image).
2. Siemens: ref data: EPI ramp sampling regrid, Nyquist ghost removal, FOV shift according to the MB2 data
3. GE: ref data FOV chop, FOV shift according to the MB2 data
4. use writecfl function to generate cfl files for BART

## 📂 Test Data

Available on Hugging Face:  
🔗 [https://huggingface.co/datasets/QiangLiu0310/Harmonization_MB_Recon](https://huggingface.co/datasets/QiangLiu0310/Harmonization_MB_Recon)

## 👤 Author

**Qiang Liu**  
*Email:* qliu30@mgh.harvard.edu  
*Affiliation:* Brigham and Women's Hospital

---

*For any questions or contributions, feel free to open an issue or submit a pull request.*
