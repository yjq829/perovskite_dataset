###############################################################################
#
#             Python 3.9
#      Decomposition Energy Calculator
#      input: element fraction vector
#             element space:
#     A:
#     B:
#     X: Cl, Br, I
#
#      output: decomposition energy value, mixing entropy included
#
#
#       future development:
# (1) combine with fomula -> element fraction vector
# (2) ref state phase library, make it global when use
#
###############################################################################

# modules dependency
import numpy
import numpy as np
import pandas
import copy
import sympy
import matplotlib.pyplot as plt

# initialize the reference data
# TOTEN for reference compounds loading
AX_ref_HSE = pandas.read_excel('ref_data.xlsx', sheet_name='AX_HSE', engine='openpyxl')
AX_ref_PBE = pandas.read_excel('ref_data.xlsx', sheet_name='AX_PBE', engine='openpyxl')
BX2_ref_HSE = pandas.read_excel('ref_data.xlsx', sheet_name='BX2_HSE', engine='openpyxl')
BX2_ref_PBE = pandas.read_excel('ref_data.xlsx', sheet_name='BX2_PBE', engine='openpyxl')
AX_ref_HSE_dict = AX_ref_HSE.set_index('sys')['HSE_ref'].to_dict()
AX_ref_PBE_dict = AX_ref_PBE.set_index('sys')['PBE_ref'].to_dict()
BX2_ref_HSE_dict = BX2_ref_HSE.set_index('sys')['HSE_ref'].to_dict()
BX2_ref_PBE_dict = BX2_ref_PBE.set_index('sys')['PBE_ref'].to_dict()
ref_HSE_dict = AX_ref_HSE_dict
ref_HSE_dict.update(BX2_ref_HSE_dict)
ref_PBE_dict = AX_ref_PBE_dict
ref_PBE_dict.update(BX2_ref_PBE_dict)


# print(ref_HSE_dict)
# print(ref_PBE_dict)


# mixing analyze
# imput an np array or list, default is list
def element_exist_list(frac_vec_input):
    frac_mix = np.array(copy.deepcopy(frac_vec_input))
    # element_exits is the index of element, not fraction
    element_exist = np.nonzero(frac_mix)[0]
    element_space = A_element + B_element + X_element
    formula = [element_space[x] for x in element_exist]
    element_frac = [frac_mix[x] for x in element_exist]
    # collect elements for each site
    # [element number in A site,element number in B site,element number in X site]
    element_mixed_num = [0, 0, 0]
    for elem_num in range(len(element_exist)):
        if element_exist[elem_num] <= 4:
            element_mixed_num[0] += 1
        elif 5 <= element_exist[elem_num] <= 10:
            element_mixed_num[1] += 1
        elif element_exist[elem_num] >= 11:
            element_mixed_num[2] += 1
    return element_mixed_num, formula, element_frac


def mixing_ana(frac_vec_input):
    frac_mix = np.array(copy.deepcopy(frac_vec_input))
    element_mixed_num, formula, element_frac = element_exist_list(frac_mix)
    # find mixing, only consider pure, amix, bmix and xmix
    # only 1 site is mixed
    if element_mixed_num == [1, 1, 1]:
        mixing = 'Pure'
    elif element_mixed_num[1:] == [1, 1] and element_mixed_num[0] != 1:
        mixing = 'Amix'
    elif element_mixed_num[0:2] == [1, 1] and element_mixed_num[2] != 1:
        mixing = 'Xmix'
    else:
        mixing = 'Bmix'
    return element_mixed_num, mixing, formula, element_frac


# decomposed phase extract
def decomp_phase_ext(element_input, formula_input, mix, element_frac):
    element_tem = copy.deepcopy(element_input)
    formula = copy.deepcopy(formula_input)
    if mix == 'Pure':
        decomp_phase = [formula[0] + formula[2], formula[1] + formula[2] + '2']
        decomp_phase_frac = [1, 1]
    elif mix == 'Amix':
        A_num = element_tem[0]
        A_decomp = [formula[x] + formula[-1] for x in range(A_num)]
        decomp_phase = A_decomp + [formula[-2] + formula[-1] + '2']
        n_element = len(element_frac)
        decomp_phase_frac = element_frac[0:n_element - 1]
    elif mix == 'Bmix':
        B_num = element_tem[1]
        B_decomp = [formula[x + 1] + formula[-1] + '2' for x in range(B_num)]
        decomp_phase = [formula[0] + formula[-1]] + B_decomp
        n_element = len(element_frac)
        decomp_phase_frac = element_frac[0:n_element - 1]
    elif mix == 'Xmix':
        X_num = element_tem[2]
        A_decomp = [formula[0] + formula[-2 + x] for x in range(X_num)]
        B_decomp = [formula[1] + formula[-2 + x] + '2' for x in range(X_num)]
        decomp_phase = A_decomp + B_decomp
        decomp_phase_frac = [x / 3 for x in element_frac[2:]] + [y / 3 for y in element_frac[2:]]
    return decomp_phase, decomp_phase_frac


def entropy_calcs(decomp_frac):
    decomp_frac_test = copy.deepcopy(decomp_frac)
    # print(decomp_frac_test)
    # kB in eV/K, using 300K as room temperature
    k_b = 8.617e-5
    T_ref = 300
    # print(element, mix, fomula)
    frac_test_a=[]
    for a in decomp_frac_test:
        if a != 0:
            frac_test_a.append(a)
    # print(frac_test_a)
    mixing_entropy = k_b * T_ref * np.dot(np.array(frac_test_a), np.log(np.array(frac_test_a)))
    return mixing_entropy


# decomposition calculation function
def decomp_calc(frac_input, TOTEN_input, functional='HSE'):
    frac = copy.deepcopy(frac_input)
    TOTEN = copy.deepcopy(TOTEN_input)
    if functional == 'HSE':
        element, mix, fomula, element_frac = mixing_ana(frac)
        # print(element, mix, fomula)
        decomp_phase, decomp_phase_frac = decomp_phase_ext(element, fomula, mix, element_frac)
        print(decomp_phase)
        phase_ref_energy = []
        # use decom_phase to search ref energy in dictionary
        for i in decomp_phase:
            phase_ref_energy.append(ref_HSE_dict[i])
        print(decomp_phase_frac)
        print(phase_ref_energy)
        mixing_entropy = entropy_calcs(decomp_phase_frac)
        print(mixing_entropy)
        decomp_energy = TOTEN - np.dot(np.array(decomp_phase_frac), np.array(phase_ref_energy)) + mixing_entropy
    elif functional == 'PBE':
        element, mix, fomula, element_frac = mixing_ana(frac)
        # print(element, mix, fomula)
        decomp_phase, decomp_phase_frac = decomp_phase_ext(element, fomula, mix, element_frac)
        print(decomp_phase)
        phase_ref_energy = []
        # use decom_phase to search ref energy in dictionary
        for i in decomp_phase:
            phase_ref_energy.append(ref_PBE_dict[i])
        print(decomp_phase_frac)
        print(phase_ref_energy)
        mixing_entropy = entropy_calcs(decomp_phase_frac)
        print(mixing_entropy)
        decomp_energy = TOTEN - np.dot(np.array(decomp_phase_frac), np.array(phase_ref_energy)) + mixing_entropy
    return decomp_energy


def decomp_phase_solve(ele_frac):
    ele_frac_solve = copy.deepcopy(ele_frac)
    X1 = ele_frac_solve[-2]
    X2 = ele_frac_solve[-1]
    x, y = sympy.symbols('x y')
    ans = sympy.linsolve([x + 2 * y - X1, (1 - x) + 2 * (1 - y) - X2], (x, y))
    return ans


def decomp_calc_solve(decomp_phase, decomp_phase_frac, TOTEN, functional='PBE'):
    decomp_phase_tem = copy.deepcopy(decomp_phase)
    decomp_phase_frac_tem = copy.deepcopy(decomp_phase_frac)
    if functional == 'PBE':
        phase_ref_energy = []
        for i in decomp_phase_tem:
            phase_ref_energy.append(ref_PBE_dict[i])
        # print(decomp_phase_tem)
        # print(phase_ref_energy)
        mixing_entropy = entropy_calcs(decomp_phase_frac_tem)
        decomp_energy = TOTEN - np.dot(np.array(decomp_phase_frac_tem), np.array(phase_ref_energy))   + mixing_entropy
    elif functional == 'HSE':
        phase_ref_energy = []
        for i in decomp_phase_tem:
            phase_ref_energy.append(ref_HSE_dict[i])
        # print(decomp_phase_tem)
        # print(phase_ref_energy)
        mixing_entropy = entropy_calcs(decomp_phase_frac_tem)
        decomp_energy = TOTEN - np.dot(np.array(decomp_phase_frac_tem), np.array(phase_ref_energy))  + mixing_entropy
    return decomp_energy


# define list for element space
A_element = ['K', 'Rb', 'Cs', 'MA', 'FA']
B_element = ['Ca', 'Sr', 'Ba', 'Ge', 'Sn', 'Pb']
X_element = ['Cl', 'Br', 'I']

if __name__ == '__main__':
  
    decomp_result = []
    sample_input = pandas.read_excel('Decomp_HSErel_SOC_calcs.xlsx', engine='openpyxl')
    sample_list = sample_input.values.tolist()
    num_sample = len(sample_list)
    decomp_result_phase = []
    for sample in range(num_sample):
        print(sample)
        test_element, test_mix, test_fomula, test_decomp_ele_frac = mixing_ana(sample_list[sample][0:14])
        test_decomp_phase, test_decomp_fake_frac = decomp_phase_ext(test_element, test_fomula, test_mix,
                                                                    test_decomp_ele_frac)
        decomp_phase_ans = decomp_phase_solve(test_decomp_ele_frac)
        print(decomp_phase_ans)
        list_decomp = []
        list_phase = []
        for i in range(10001):
            y = 0.0001 * i
            decomp_phase_frac = [eval(str(decomp_phase_ans.args[0][0])),
                                 1 - eval(str(decomp_phase_ans.args[0][0])),
                                 eval(str(decomp_phase_ans.args[0][1])),
                                 1 - eval(str(decomp_phase_ans.args[0][1]))]
            decomp_phase_frac = [round(r, 4) for r in decomp_phase_frac]
            if (numpy.array(decomp_phase_frac) >= 0).all():
                test_decomp = decomp_calc_solve(test_decomp_phase, decomp_phase_frac, sample_list[sample][-1], 'HSE')
                list_decomp.append(test_decomp)
                list_phase.append(decomp_phase_frac)
        decomp_result.append(max(list_decomp))
        max_index=list_decomp.index(max(list_decomp))
        decomp_result_phase.append(list_phase[max_index])
        print(list_phase[max_index])
    result_out = sample_input
    result_out['decomp'] = decomp_result
    # result_out['p1']=decomp_result_phase[:,0]
    # result_out['p2']=decomp_result_phase[:,1]
    # result_out['p3']=decomp_result_phase[:,2]
    # result_out['p4']=decomp_result_phase[:,3]
    result_out.to_excel('Decomp_calcs_HSErel_SOC_results.xlsx', engine='openpyxl')
