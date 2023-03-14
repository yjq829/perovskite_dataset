# perovskite_dataset

**[A]- Decomposition Energy:**
==============================
**(1) Purpose:

Scripts used to calculate decomposition energy for all halide perovskits.

Using Python 3.9, pands, numpy, sympy and matplotlib required.

Refenece energy for all possible decomposed phase is restored in ref_data.xlsx, which will be loaded each running.

[i] DecompositionEnergyCalculator.py

general decomposition energy calculator for all halide perovskites. Used in all decomposition energy calculation shown in the papre

[ii] DecompositionEnergy_Xmix.py

Special verision of decomposioin energy calculator for X-mixed halide perovskites. 
Since the X-mixed halide perovskites have multipul decomposed phase possiblity, this code is used to find the most possible decomposition reaction.
Used in section: "Decomposition Enery correction for X site mixed perovskites"


**(2) input files:
The input files to use the decomposition calculator is an spreadsheet. 
The format of input file is a 14 dimentional compositon vector as 1-14th column and total energy p.f.u calcualted from DFT in 15th column.

The example input files can be found as 

[i] decom_calcs.xlsx -> used for DecompositionEnergyCalculator.py

[ii] Decomp_HSErel_SOC_calcs.xlsx -> used for DecompositionEnergy_Xmix.py

**(3) output files:
Both script will write an spreadsheet as output. the last colum will show the calucalted decomposition energy.

**[B] SLME calculations:**
==========================

The SLME cacluation is based on the work of L. Yu, A. Zunger, Phys. Rev. Lett. 108, 068701 (2012). https://doi.org/10.1103/PhysRevLett.108.068701

The source code can be found at: https://github.com/ldwillia/SL3ME

**(1) Scripts and Dependency:

[i] SLME_shift_fromdata.py -> MAIN running entrance

This script calculating SLME for halide perovskites. For HSE funtionals, it's too expensive to calculate absorption spectrum. So we shifted PBE absorption spectrom with band gap differenct to approxiamte HSE absorption spectrum.

[ii] SL3ME.py -> SLME calculation function dependency
This SLME cacluation is published by the work of L. Yu, A. Zunger, Phys. Rev. Lett. 108, 068701 (2012). https://doi.org/10.1103/PhysRevLett.108.068701 
We import the functions from this modules

[iii] am1.5G.dat
This is the "Global Tilt" spectra used as reference spectrum for SL3ME.py module. It's also assigned by the authros of SL3ME.

**(2) input files:
[i] data.xlsx
This is the input file for caculation SLME of our dataset. The first column is the poervskite index refer to the formula in the dataset. The second column is the formula. The third column is the PBE band gap collected from DFT. The forth column is the band gap from target functional.
Here we use results from HSE as an example input here.





