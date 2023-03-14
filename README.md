# perovskite_dataset
[A]- Decomposition Energy:

(1) Purpose:

Scripts used to calculate decomposition energy for all halide perovskits.
Using Python 3.9, pands, numpy, sympy and matplotlib required.

[i] DecompositionEnergyCalculator.py
general decomposition energy calculator for all halide perovskites. Used in all decomposition energy calculation shown in the papre

[ii] DecompositionEnergy_Xmix.py
Special verision of decomposioin energy calculator for X-mixed halide perovskites. 
Since the X-mixed halide perovskites have multipul decomposed phase possiblity, this code is used to find the most possible decomposition reaction.
Used in section: "Decomposition Enery correction for X site mixed perovskites"


(2) input files:
The input files to use the decomposition calculator is an spreadsheet. 
The format of input file is a 14 dimentional compositon vector as 1-14th column and total energy p.f.u calcualted from DFT in 15th column.
The example input files can be found as 
[i] decom_calcs.xlsx -> used for DecompositionEnergyCalculator.py

[ii] Decomp_HSErel_SOC_calcs.xlsx -> used for DecompositionEnergy_Xmix.py

(3) output files:
Both script will write an spreadsheet as output. the last colum will show the calucalted decomposition energy.

[B]
