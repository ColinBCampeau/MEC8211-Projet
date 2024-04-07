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

L2_verlet = f_L2(vec_x, x_anal)
print('L2 Verlet :', L2_verlet, 'pour dt =', prm.dt)

#%%========================= Méthode Euler Explicite ==========================%%#

t_vect, x_vect, v_vect = euler(prm)

plt.plot(t_vect,x_vect,'b:',label='euler explicite')
plt.legend()

L2_euler = f_L2(x_vect, x_anal)
print('L2 Euler :', L2_euler, 'pour dt =', prm.dt)

