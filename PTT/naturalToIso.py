# Python3 script to convert from natural compostion to isotopic composition
# Author : R. Guasch
# Date : 2024-06-18
#
# This script assumes that mass fractions are given for natural elements. Converts them using the proportion of isotopes in the natural element.

import numpy as np
check_total = True
#natural_isotopes = [14000, 24000, 25000, 26000, 28000] # Si, Cr, Mn, Fe, Ni
mass_fractions = {"14000": 5.1000E-03, "24000": 1.7400E-01, "25000": 1.9900E-02, "26000": 6.8400E-01 , "28000": 1.1700E-01}
natural_proportions = {"14028":0.922545, "14029": 0.04672, "14030": 0.030735, "24050": 0.04345, "24052": 0.83789, "24053": 0.09501, "24054": 0.02365, "25055": 1.0, "26054": 0.05845, "26056": 0.91754, "26057": 0.02119, "26058": 0.00282, "28058": 0.680769, "28060": 0.262231, "28061": 0.011399, "28062": 0.036345, "28064": 0.009256}

mass_fractions_per_isotope = {}
for iso in natural_proportions.keys():
    element = str(int(iso[:2]))+"000"
    mass_fractions_per_isotope[iso] = mass_fractions[element]*natural_proportions[iso]

if check_total:
    print(mass_fractions_per_isotope)
    sum_14 = 0
    sum_24 = 0
    sum_25 = 0
    sum_26 = 0
    sum_28 = 0
    for iso in mass_fractions_per_isotope.keys():
        if iso[:2] == "14":
            sum_14 += mass_fractions_per_isotope[iso]
        elif iso[:2] == "24":
            sum_24 += mass_fractions_per_isotope[iso]
        elif iso[:2] == "25":
            sum_25 += mass_fractions_per_isotope[iso]
        elif iso[:2] == "26":
            sum_26 += mass_fractions_per_isotope[iso]
        elif iso[:2] == "28":
            sum_28 += mass_fractions_per_isotope[iso]
    #print(sum_14, sum_24, sum_25, sum_26, sum_28)
    sum_all = sum_14 + sum_24 + sum_25 + sum_26 + sum_28
    print(sum_all)

# Output:
for iso in mass_fractions_per_isotope.keys():
    print(f"{iso}.05c -{mass_fractions_per_isotope[iso]}") 