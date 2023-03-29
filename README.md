# perovskite_dataset


**[A]- Decomposition Energy:**
==============================

**(1) Purpose:**

Scripts used to calculate decomposition energy for all halide perovskites, applicable to both pure and alloyed compositions.

Python 3.9, pandas, numpy, sympy and matplotlib required.

Reference energies for all possible decomposed phases (AX, BX2, and A/B/X species) are stored in ref_data.xlsx, which will be read every time the script is run.

[i] DecompositionEnergyCalculator.py
General decomposition energy calculator for all halide perovskite compositions, used to generate all values presented in the paper.

[ii] DecompositionEnergy_Xmix.py
Special version of decomposition energy calculator for X-site mixed halide perovskites, which generally have multiple possibilities for decomposed phases. This script finds the most likely (ABX3 -> AX + BX2) decomposition reaction and calculates the decomposition energy based on the appropriate phases. Process is described in the SI under "Decomposition Energy correction for X site mixed perovskites".


**(2) input files:**

The input files to use for the decomposition calculator is a spreadsheet. 
The format of the input file is a 14-dimensional composition vector (columns 1 to 14) and total DFT energy per ABX3 functional unit as the 15th column.

Example input files include:

[i] decomp_calcs.xlsx -> used for DecompositionEnergyCalculator.py

[ii] Decomp_HSErel_SOC_calcs.xlsx -> used for DecompositionEnergy_Xmix.py


**(3) output files:**

Both scripts will write out a spreadsheet as output, with the last column showing the calculated decomposition energy.

**[B] SLME calculations:**
==========================

The SLME calculation is based on the work of L. Yu, A. Zunger, Phys. Rev. Lett. 108, 068701 (2012). https://doi.org/10.1103/PhysRevLett.108.068701

The source code can be found at: https://github.com/ldwillia/SL3ME

**(1) Scripts and dependency:**

[i] SLME_shift_fromdata.py -> MAIN running entrance

This script performs the entire task of calculating the SLME for any given compound. It is a straightforward task for the PBE computed absorption spectrum, but for the HSE functionals, we do not calculate the absorption spectrum but rather obtain it by shifting the PBE computed spectrum by the difference between the PBE and HSE band gaps.

[ii] SL3ME.py -> SLME calculation function dependency

This is the original code developed and released by the authors of this publication: L. Yu, A. Zunger, Phys. Rev. Lett. 108, 068701 (2012). https://doi.org/10.1103/PhysRevLett.108.068701 
We import the functions from this script as a module.

[iii] am1.5G.dat

This is the "Global Tilt" spectra used as reference spectrum in SL3ME.py.

**(2) input files:**

[i] data.xlsx

This is the input file for calculating SLME for any compound in our dataset. The first column is the perovskite index which refers to its unique chemical composition. The second column is the chemical formula, the third column is the PBE band gap, and the forth column is the band gap from the target functional (PBE or any of the HSE functionals).
Here, we use results from a particular HSE functional as an example input.

[ii] strut_loptics_550 folder

This folder includes all the PBE computed absorption spectra of halide perovskite compounds in our dataset. The file label is the same as the corresponding perovskite index in the spreadsheet perovs_data_final.xlsx.




**[C] Pearson Correlation Calculations:**
==========================

**(1) Scripts and Dependency**

[i] Folders:

All folders are named after the functionals corresponding to the DFT data.

There are 5 files inside:

*_data.csv (e.g. PBE_data.csv) -> File containing PBE computed properties and descriptors for calculating Pearson coefficients of linear correlation.

Func_corr.py (e.g. PBE_corr.py) -> First step script, generates Pearson Correlation values.

test.csv -> output of Pearson correlation values generated by Func_corr.py.

Func_heatmap.py. -> plots correlation values as a heatmap, using test.csv as input.

Corr.xlsx -> label file showing names of all descriptors.

**(2) Input Formats**
The input format is *_data.csv. For every functional, there is an example file provided.

For PBE and HSErel functionals, the format of the .csv file is:

Formula, type of mixing, 4 columns of properties, 14 columns of species fraction descriptors, 36 columns of elemental property descriptors.

For HSErel-SOC and HSE-PBE-SOC functionals, the format is:

Formula, type of mixing, 3 columns of properties, 14 columns of species fraction descriptors, 36 columns of elemental property descriptors.

For more details about the descriptors, please refer to our paper.





![image](https://user-images.githubusercontent.com/32602669/228574400-9792d8ef-5227-4a74-927b-fd9ad574e273.png)
