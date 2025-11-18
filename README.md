# Multi-band Reconstruction for Recon Reproducibility Study

This repository contains code and workflows for evaluating multi-band (MB) reconstruction methods in the context of diffusion MRI reconstruction reproducibility.

## 🛠 Reconstruction Methods

1. SMS-SENSE (Custom MATLAB Implementation)

2. SMS-SENSE using BART (L1-ESPIRiT)

3. Split Slice-GRAPPA (Custom MATLAB Implementation)

⸻

## 🧾 Reconstruction Pipelines

🔹 Split Slice-GRAPPA

Step 1: Calibration (Reference Data)
	•	a. GRAPPA ACS data extraction
	•	b. EPI ramp sampling regridding
	•	c. Nyquist ghost correction (odd-even readout phase correction using Linear Phase Correction, LPC)
	•	d. GRAPPA kernel training

Step 2: Imaging Data Reconstruction
	•	a. EPI ramp sampling regridding
	•	b. Nyquist ghost correction (LPC)
	•	c. Apply trained GRAPPA kernels

Step 3: Partial Fourier Reconstruction
	•	Homodyne reconstruction using a reconhd function

Step 4: Coil Combination
	•	Sum-of-squares method along the coil dimension

⸻

🔹 L1-ESPIRiT (BART)

Step 1: Reference Data Extraction for Sensitivity Map Estimation
	•	Siemens: EPI-based reference scan
	•	GE: Low-resolution GRE image

Step 2: Preprocessing
	•	Siemens:
	•	EPI ramp sampling regridding
	•	Nyquist ghost correction
	•	FOV shift to match MB2 acquisition
	•	GE:
	•	Image-domain FOV cropping of GRE image
	•	FOV shift to match MB2 acquisition

Step 3: Coil Sensitivity Estimation
	•	Generate .cfl and .hdr files using writecfl for BART
	•	Use BART’s ecalib and pics tools for reconstruction

⸻

## 📂 Test Data

Test data and reference reconstructions are provided for both Siemens and GE datasets. Please refer to the data/ folder and the example scripts in scripts/ for reproducibility.

Available on Hugging Face:  
🔗 [https://huggingface.co/datasets/QiangLiu0310/Harmonization_MB_Recon](https://huggingface.co/datasets/QiangLiu0310/Harmonization_MB_Recon)

## 👤 Author

**Qiang Liu**  
*Email:* qliu30@mgh.harvard.edu  
*Affiliation:* Brigham and Women's Hospital

## Recently Published *Medical Physics* Paper

```bibtex
@article{article,
  author = {Liu, Qiang and Zhu, Ante and Wang, Xiaoqing and Erdogmus, Deniz and Westin, Carl-Fredrik and O'Donnell, Lauren and Bilgic, Berkin and Ning, Lipeng and Rathi, Yogesh},
  year = {2025},
  month = {10},
  pages = {},
  title = {Effect of a consistent reconstruction algorithm on inter‐scanner reproducibility in diffusion MRI},
  volume = {52},
  journal = {Medical Physics},
  doi = {10.1002/mp.70096}
}

