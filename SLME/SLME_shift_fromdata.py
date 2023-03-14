##########################################################
#
# Python Script for Shifting adsorption spectrum with Bg_HSE - Bg_PBE
#
##########################################################

import numpy as np
import copy
from SL3ME import calculate_SLME
import pandas




def bg_hse_pbe_diff_fromdata(data_file):
    perovs_data = pandas.read_excel(data_file, engine='openpyxl')
    perovs_data_use = copy.deepcopy(perovs_data)
    bg_pbe_read = perovs_data_use.PBE_bandgap[:]
    bg_target_read = perovs_data_use.target_bandgap[:]
    bg_diff_read = bg_target_read - bg_pbe_read
    return perovs_data_use, bg_diff_read


def absorbance_load(absorb_file):
    absorb_file_name = absorb_file
    material_ev_for_absorbance_data, material_absorbance_data = np.loadtxt(absorb_file_name, usecols=[0, 1],
                                                                           unpack=True)
    # absorbance in 2nd column
    return material_ev_for_absorbance_data, material_absorbance_data


def build_absorbance_shift(bd_diff, material_ev, material_absorbance):
    material_ev_input = copy.deepcopy(material_ev)
    material_absorbance_input = copy.deepcopy(material_absorbance)
    shift_ev_array = [0, bd_diff - 0.01]
    shift_absorbance_array = [0, 0.001]
    material_ev_input_shifted = [i + bd_diff for i in material_ev_input]
    material_ev_shift = np.concatenate((shift_ev_array, material_ev_input_shifted))
    material_absorbance_shift = np.concatenate((shift_absorbance_array, material_absorbance_input))
    return material_ev_shift, material_absorbance_shift


if __name__ == "__main__":
    #SLME calculation function constants
    # defaultly defined by SL3ME.py module, L Yu et al.
    material_direct_allowed_gap = 2.304
    material_indirect_gap = 2.294
    perovs_hse_data, bg_diff = bg_hse_pbe_diff_fromdata('data.xlsx')
    perovs_hse_data.insert(perovs_hse_data.shape[1],"bandgap_diff",bg_diff)
    shifted_slme_data = []
    for x in range(len(perovs_hse_data.perovs_index[:])):
        absorb_input = "strut_loptics_550/" + str(perovs_hse_data.perovs_index[x]) + ".dat"
        material_eV, material_absorbance = absorbance_load(absorb_input)
        material_ev_shifted, material_absorbance_shifted = build_absorbance_shift(bg_diff[x], material_eV,
                                                                                  material_absorbance)
        SLME = calculate_SLME(material_eV_for_absorbance_data=material_ev_shifted,
                              material_absorbance_data=material_absorbance_shifted,
                              material_direct_allowed_gap=material_direct_allowed_gap,
                              material_indirect_gap=material_indirect_gap,
                              thickness=5E-6, T=293.15)
        print('File :', str(perovs_hse_data.perovs_index[x]))
        print('Shifted Standard SLME :', SLME)
        shifted_slme_data.append(SLME)

        ###################################################
        # SLME calculation done
        ###################################################
    shifted_slme_data = np.array(shifted_slme_data)
    perovs_hse_data.insert(perovs_hse_data.shape[1], "HSE_slme", shifted_slme_data)
    perovs_hse_data.to_excel('perovs_data_hse_shifted.xlsx')
