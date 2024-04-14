"""

PROJET MEC8211
                            >FICHIER MAIN<

AUTEUR: EROJ MOHAMMAD ISHOQ, COLIN BISSONNETTE-CAMPEAU, TRISTAN ANCEL-SÉGUIN
CRÉATION: 2 AVRIL 2024
MISE À JOUR: 14 AVRIL 2024

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

plt.figure
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
plt.show()

L2_euler = f_L2(x_vect, x_anal)
print('L2 Euler :', L2_euler, 'pour dt =', prm.dt)

#%%========================= Calcul de la fréquence ==========================%%#

freq = frequence(vec_x)
print('Fréquence :',freq,'Hz')

#%%========================= Calcul de u_num ==========================%%#

r = 2
u_num = incert_num(r,prm)
print('u_num :',u_num, 'pour dt_1 =', prm.dt, 'et r =', r)

#%%========================= Calcul de u_input et S ==========================%%#

n_lhs = 100
u_input,S = incert_input(n_lhs,prm)
print('u_input :',u_input, 'pour dt =', prm.dt, 'et n_lhs =', n_lhs)
print('S :',S, 'Hz')

#%%========================= Calcul de u_D ==========================%%#

u_D = incert_D(prm)
print('u_D :',u_D)
