# Author : R. Guasch
# Date : 2024/09/04
# Purpose : plot D5-S2 bias for HOM_U5_U8 keff with RSE+CORR as a function of epsilon_svd
import numpy as np
import matplotlib.pyplot as plt

epsilons_svd_U5_U8 = [1.00E-01,5.00E-02,1.00E-02,5.00E-03,1.00E-03,5.00E-04,1.00E-04,5.00E-05,1.00E-05]
error_U5_U8 = [51.58,49.75,46.94,46.89,46.64,46.61,46.51,46.54,46.20]

dict_of_errors_U5_U8_CORR = {
"1.00E-01" : 45.433,
"9.00E-02" : 17.836,
"7.00E-02" : 14.867,
"5.00E-02" : -80.928,
"3.00E-02" : -79.653,
"1.00E-02" : 116.06,
"9.00E-03" : 116.004,
"7.00E-03" : 162.591,
"5.00E-03" : 161.459,
"3.00E-03" : 94.36847,
"1.00E-03": 263.59,
"9.00E-04" : 204.8,
"7.00E-04" : 206.12,
"5.00E-04" : 206.84,
"3.00E-04" : 256.82,
"1.00E-04" : 389.459,
"5.00E-05" : 417.413,
"9.00E-05" : 389.69,
"7.00E-05" : 385.76,
"3.00E-05" : 668.69,
"1.00E-05" : 770.75 }

epsilons_svd_U5_U8_CORR = [float(key) for key in dict_of_errors_U5_U8_CORR.keys()]
error_U5_U8_CORR = [dict_of_errors_U5_U8_CORR[key] for key in dict_of_errors_U5_U8_CORR.keys()]

# --- PLOTTING linear

plt.figure(figsize=(10,6))
plt.scatter(epsilons_svd_U5_U8_CORR,error_U5_U8_CORR,marker = "x", color='r',label='U5_U8_CORR')
plt.scatter(epsilons_svd_U5_U8,error_U5_U8, marker="x", color='b',label='U5_U8')
plt.xlabel('epsilon_svd')
plt.ylabel('D5-S2 bias on keff (pcm)')
plt.legend()
plt.title('D5-S2 bias for HOM_U5_U8 keff as a function of epsilon_svd')
plt.grid()
plt.savefig('D5-S2_bias_HOM_U5_U8_vs_epsilon_svd.png')
plt.show()
# --- END OF PLOTTING linear

# --- PLOTTING log

plt.figure(figsize=(10,6))
plt.scatter(-np.log(epsilons_svd_U5_U8_CORR)/np.log(10),error_U5_U8_CORR,marker = "x", color='r',label='U5_U8_CORR')
plt.scatter(-np.log(epsilons_svd_U5_U8)/np.log(10),error_U5_U8, marker="x", color='b',label='U5_U8')
plt.xlabel('-log(epsilon_svd)')
plt.ylabel('D5-S2 bias on keff (pcm)')
plt.legend()
plt.title('D5-S2 bias for HOM_U5_U8 keff vs logscale of epsilon_svd')
plt.grid()
plt.savefig('D5-S2_bias_HOM_U5_U8_vs_epsilon_svd_log.png')

epsilons = [1.00E-5, 3.00E-5, 5.00E-5, 7.00E-5, 9.00E-5, 1.00E-4, 3.00E-4, 5.00E-4, 7.00E-4, 9.00E-4, 1.00E-3, 3.00E-3, 5.00E-3, 7.00E-3, 9.00E-3, 1.00E-2, 3.00E-2, 5.00E-2, 7.00E-2, 9.00E-2, 1.00E-1]
