# -*- coding: utf-8 -*-
"""
Plot reference SN and nodal diffusion flux curves
"""

import lcm,numpy
import matplotlib.pyplot as plt

# access SN flux
print('access SN flux')
my_lcm=lcm.new('LCM_INP','FLX_CB1_plot.txt')
my_lcm.lib()
print('object signature=', my_lcm['SIGNATURE'])
state=my_lcm['STATE-VECTOR']
ngr=state[0]
nel=state[1]
print('ngr=', ngr, 'nel=', nel)
volume=my_lcm['VOLUME']
flux_case1 = numpy.empty([ngr, nel])
hfac_case1 = numpy.empty([ngr, nel])
power_sn=0.0
for ig in range(ngr):
    o2 = my_lcm['GROUP'][ig].keys()
    flux_case1[ig] = my_lcm['GROUP'][ig]['FLUX-INTG']
    if 'H-FACTOR' in o2:
        hfac_case1[ig] = my_lcm['GROUP'][ig]['H-FACTOR']
    else:
        hfac_case1[ig] = my_lcm['GROUP'][ig]['NUSIGF']
    for iel in range(nel):
        power_sn=power_sn+hfac_case1[ig][iel]*flux_case1[ig][iel]
        flux_case1[ig][iel]=flux_case1[ig][iel]/volume[iel]
xmidsn = numpy.empty([nel])
for iel in range(nel):
    if iel==0:
        xmidsn[iel]=0.0
    else:
        xmidsn[iel]=xmidsn[iel-1]+(volume[iel-1]+volume[iel])/2
flux_case1=flux_case1/power_sn # normalize power to one
print('power_sn=', power_sn)
del my_lcm

# access nodal flux
print('access nodal flux')
my_lcm=lcm.new('LCM_INP','VAL_CB1_plot.txt')
my_lcm.lib()
print('object signature=', my_lcm['SIGNATURE'])
state=my_lcm['STATE-VECTOR']
ngr=state[0]
nx=state[1]
print('ngr=', ngr, 'nx=', nx)
xmid=my_lcm['MXI']
flx1 = numpy.empty([ngr, nx])
for ig in range(ngr):
    flx1[ig][:nx] = my_lcm['FLUX'][ig][:nx]

# plot values
side=21.5
major_ticks = [0., side, 2.0*side, 3.0*side]
plt.rcParams["figure.figsize"] = (10,3)
fig, (ax1, ax2) = plt.subplots(1,2)
fig.suptitle("DF-RT: fast and thermal fluxes")
ax1.plot(xmidsn,flux_case1[0], "-b", label="SN reference flux", linewidth=0.7)
ax1.plot(xmid,flx1[0], "-r", label="Raviart-Thomas flux", linewidth=0.7)
ax1.set_xticks(major_ticks)
ax1.grid(True)
ax1.legend(loc="upper right")
ax1.set_xlabel("position (cm)")
ax2.plot(xmidsn,flux_case1[1], "-b", label="SN reference flux", linewidth=0.7)
ax2.plot(xmid,flx1[1], "-r", label="Raviart-Thomas flux", linewidth=0.7)
ax2.set_xticks(major_ticks)
ax2.grid(True)
ax2.legend(loc="upper right")
ax2.set_xlabel("position (cm)")
plt.savefig('SmallCore_BaffRefl.eps', bbox_inches='tight', format='eps', dpi=300)
#plt.show()
