"""

PROJET MEC8211
                    >FICHIER GRAPHIQUE RÉSULTATS<

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




#%%========================= Graphiques ==========================%%#

zeta = np.array([0, 0.1, 0.5])

prm.m = 1

prm.k = 100

prm.w_n = np.sqrt(prm.k/prm.m)

prm.x0 = 0

prm.v0 = 10

prm.t_fin = 5

for i in range(len(zeta)):
    
    prm.zeta = zeta[i]

    prm.c =  prm.zeta*(2*prm.w_n*prm.m)
    
    prm.w_d = prm.w_n*np.sqrt(1-prm.zeta**2) # Fréquence amortie [Hz]
    
    prm.a0 = -1/prm.m*(prm.c*prm.v0+prm.k*prm.x0) # Accélération initiale [m/s^2]
    
    x_anal, t = sol_analytique(prm)
    
    plt.plot(t,x_anal,'-',label='Analytique')
    
    vec_x, vec_v, vec_a, t = verlet(prm)

    plt.plot(t,vec_x,'--',label='Verlet')
    
    t_vect, x_vect, v_vect = euler(prm)

    plt.plot(t_vect,x_vect,':',label='Euler explicite')
    
    plt.legend()
    
    plt.grid()
    
    plt.title(r"Déplacement en fonction du temps pour $\zeta  = $"+str(zeta[i]))
    
    plt.xlabel(r"Temps [s]")
    
    plt.ylabel(r"Déplacement [m]")
    
    plt.savefig("graph_zeta_"+str(zeta[i])+".png", dpi=300,bbox_inches='tight')
    
    plt.show()



