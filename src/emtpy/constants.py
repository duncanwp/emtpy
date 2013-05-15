__author__ = 'pard'

eV = 1.6022E-19 # C
me = 9.109E-31 # Kg
hbar = 1.0546E-34 # J.s
h = 6.62606876E-34 # J.s
cvacuum = 299792458 # m/s
epsilon0 = 8.85E-12 #
units = 0.0761996 #=((hbar**2)*1E18)/(me*eV)

hbarev = hbar/eV # eV.s
kb = 1.3806503E-23 # J / kelvin
kbev = kb/eV # eV / kelvin

maxn = 6.0E6
maxnev = 100
maxncv = 200
ldv = maxn
iterations = 1E6
#Eigen problem solver parameters

maxpos = 1E6
specsize = 10000
NoTemps = 1
