"""

PROJET MEC8211
                            >FICHIER CLASSE<

AUTEUR: EROJ MOHAMMAD ISHOQ, COLIN BISSONNETTE-CAMPEAU, TRISTAN ANCEL-SÉGUIN
CRÉATION: 2 AVRIL 2024
MISE À JOUR: 6 AVRIL 2024

"""
    
#%%===================== IMPORTATION DES MODULES ==========================%%#


import numpy as np
import matplotlib.pyplot as plt
try:
    from fonction import *
    from classe import *
except:
    pass

#%%========================= Solution analytique ==========================%%#


x_anal, t = sol_analytique(prm)

plt.plot(t,x_anal,'k-',label='Analytique, zeta = ' + str(prm.zeta))



#%%========================= Méthode Velocity Verlet ==========================%%#

vec_x, vec_v, vec_a, t = verlet(prm)

plt.plot(t,vec_x,'r--',label='Verlet')
plt.legend()

L2 = f_L2(vec_x, x_anal)
print(L2)

