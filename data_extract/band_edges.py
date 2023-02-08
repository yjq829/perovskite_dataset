import sys
import numpy as np
from numpy import array as npa
 
import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, Procar
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.electronic_structure.core import Spin, Orbital
 


bands = Vasprun("vasprun.xml")
bs = bands.get_band_structure()

cbm_bs = bs.get_cbm()
cbm = cbm_bs.get('energy')
text_file = open("cbm_bs.txt", "w")
text_file.write("CBM: %s" % cbm)
text_file.close()

vbm_bs = bs.get_vbm()
vbm = vbm_bs.get('energy')
text_file = open("vbm_bs.txt", "w")
text_file.write("VBM: %s" % vbm)
text_file.close()

dos = bands.complete_dos.get_cbm_vbm()
dos_cbm = dos[0]
text_file = open("cbm_dos.txt", "w")
text_file.write("DOS CBM: %s" % dos_cbm)
text_file.close()
dos_vbm = dos[1]
text_file = open("vbm_dos.txt", "w")
text_file.write("DOS VBM: %s" % dos_vbm)
text_file.close()

#gap_bs = cbm - vbm
gap_dos = dos_cbm - dos_vbm

gap = bs.get_band_gap().get('energy')
dir_gap = bs.get_direct_band_gap()

text_file = open("gap.txt", "w")
text_file.write("BS_gap: %s" % gap + '\n')
text_file.write("Direct_gap: %s" % dir_gap + '\n')
text_file.write("DOS_gap: %s" % gap_dos + '\n')
text_file.close()


aa = bands.final_structure.lattice.a
bb = bands.final_structure.lattice.b
cc = bands.final_structure.lattice.c
alpha = bands.final_structure.lattice.alpha
beta = bands.final_structure.lattice.beta
gamma = bands.final_structure.lattice.gamma

text_file = open("aa.txt", "w")
text_file.write("a: %s" % aa + '\n')
text_file.close()
text_file = open("bb.txt", "w")
text_file.write("b: %s" % bb + '\n')
text_file.close()
text_file = open("cc.txt", "w")
text_file.write("c: %s" % cc + '\n')
text_file.close()
text_file = open("alpha.txt", "w")
text_file.write("alpha: %s" % alpha + '\n')
text_file.close()
text_file = open("beta.txt", "w")
text_file.write("beta: %s" % beta + '\n')
text_file.close()
text_file = open("gamma.txt", "w")
text_file.write("gamma: %s" % gamma + '\n')
text_file.close()

vol = bands.final_structure.lattice.volume

text_file = open("vol.txt", "w")
text_file.write("volume: %s" % vol + '\n')
text_file.close()


toten = bands.final_energy

text_file = open("toten.txt", "w")
text_file.write("total energy: %s" % toten + '\n')
text_file.close()




