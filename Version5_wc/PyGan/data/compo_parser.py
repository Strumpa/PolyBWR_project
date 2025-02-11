### COMPO parser to post treat DRAGON5 outputs
# Author : R. Guasch
# Date : 05/02/2025

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import serpentTools as st
import re
import sys
import lifo
import lcm
import cle2000 


def open_compo_in_PyGan(compo_file_name, path_to_compo):
    """
    compo_file_name : str : name of the compo file to be opened
    path_to_compo : str : path to the compo file : must be included in .access script with ln -s command
    """
    cwd = os.getcwd()
    print('Opening COMPO file :', compo_file_name)
    print(cwd)
    os.chdir(path_to_compo)
    cwd = os.getcwd()
    pyCOMPO = lcm.new('LCM_INP', compo_file_name, impx=0)
    os.chdir(cwd)
    return pyCOMPO

def parse_compo(py_compo_obj, targets, dirs, params={}):
    """
    compo_obj : obj : compo object from PyGan
    targets : list : list of physical quantities to be parsed
    dirs : list : list of directories to be parsed / must be included in compo
    params : dict/of dict : 1st level keys = dirs, 2nd level keys = params to be parsed, values = values to be assigned to params
    """

    for DIR in dirs:
        print('Parsing directory :', DIR)
        # Get the number of mixtures
        #py_compo_obj[DIR]['MIXTURES'].keys()
        keys0 = py_compo_obj[DIR].keys()
        print('Keys0 :', keys0)
        mixtures = True
        glob = False
        if mixtures:
            number_of_mixtures = py_compo_obj[DIR]['MIXTURES'].len()
            print('Number of mixtures :', number_of_mixtures)
            for mix_idx in range(number_of_mixtures):
                print('Mixture index :', mix_idx)
                keys1 = py_compo_obj[DIR]['MIXTURES'][mix_idx].keys()
                print('MIXTURES :', keys1)
                print(f"CLACULATIONS[0] are : {py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][0].keys()}")
                print(f"CLACULATIONS[1] are : {py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][1].keys()}")
                num_calculations = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'].len()
                print('Number of calculations :', num_calculations)
                for calc_idx in range(num_calculations):
                    print(f"calculation index : {calc_idx}")
                    keys2 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx].keys()
                    print('MIXTURES/CALCULATIONS :', keys2)
                    values1 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['ISOTOPESDENS']
                    print('MIXTURES/CALCULATIONS/ISOTOPEDENS :', values1)
                    values2 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['ISOTOPESMIX']
                    print('MIXTURES/CALCULATIONS/ISOTOPESMIX :', values2)
                    values3 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['ISOTOPESTYPE']
                    print('MIXTURES/CALCULATIONS/ISOTOPESTYPE :', values3)
                    values4 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['ISOTOPESVOL']
                    print('MIXTURES/CALCULATIONS/ISOTOPESVOL :', values4)
                    values4 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['MIXTURESVOL']
                    print('MIXTURES/CALCULATIONS/MIXTURESVOL :', values4)
                    values4 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['K-EFFECTIVE']
                    print('MIXTURES/CALCULATIONS/K-EFFECTIVE :', values4)
                    values4 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['ENERGY']
                    print('MIXTURES/CALCULATIONS/ENERGY :', values4)
                    values4 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['DELTAU']
                    print('MIXTURES/CALCULATIONS/DELTAU :', values4)
                    num_iso_in_CALCULATION = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['ISOTOPESLIST'].len()
                    print(f'number of isotopes in mixture {mix_idx} in calculation {calc_idx} :', num_iso_in_CALCULATION)
                    for iso_idx in range(num_iso_in_CALCULATION):
                        keys = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['ISOTOPESLIST'][iso_idx].keys()
                        print('MIXTURES/CALCULATIONS/ISOTOPESLIST :', keys)
                        values6 = py_compo_obj[DIR]['MIXTURES'][mix_idx]['CALCULATIONS'][calc_idx]['ISOTOPESLIST'][iso_idx]['ALIAS']
                        print('MIXTURES/CALCULATIONS/ISOTOPESLIST/ALIAS :', values6)
        if glob:
            keys1 = py_compo_obj[DIR]["GLOBAL"].keys()
            print('Keys1 :', keys1)
            keys2 = py_compo_obj[DIR]["GLOBAL"]["PARPAD"][0]  
            print('Keys2 :', keys2)
            keys2 = py_compo_obj[DIR]["GLOBAL"]["PARCAD"][0]  
            print('Keys2 :', keys2)

        
        # Get the number of isotopes per mixture
        num_iso = np.shape(py_compo_obj[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
        print('Number of isotopes per mixture :', num_iso)




if __name__ == '__main__':
    compo_file_name = 'COMPO_2x2_UOX'
    path_to_compo = 'Linux_aarch64'
    pyCOMPO = open_compo_in_PyGan(compo_file_name, path_to_compo)
    targets = ['NFTOT']
    dirs = ['HOM1g', 'HOM2g', 'C1g', 'C2g']
    parse_compo(pyCOMPO, targets, dirs)
    print('COMPO parsing completed -- normal comlpetion of PyGan/compo_parser.py')
