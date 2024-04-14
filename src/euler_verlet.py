"""

OUTIL NUMÉRIQUE DE SIMULATION ...

            >FICHIER EULER EXPLICITE VS VERLET<

AUTEUR: EROJ MOHAMMAD ISHOQ, COLIN BISSONNETTE-CAMPEAU, TRISTAN ANCEL-SÉGUIN
CRÉATION: 2 MARS 2024
MISE À JOUR: 6 MARS 2024

"""

    
#%%===================== IMPORTATION DES MODULES ==========================%%#


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pandas as pd
import scipy as sp
import time
try:
    from fonction import *
    from classe import *
except:
    pass



#%%============= ANALYSE TEMPS DE CALCUL ===============%%#

delta_t = np.array([1e-1,1e-2,1e-3,1e-4, 1e-5])

time_verlet = np.array([])

time_euler = np.array([])

for i in range(0, len(delta_t)):
    
    prm.dt = delta_t[i]
    
    #%% verlet
    
    start_time = time.time()
    
    vec_x, vec_v, vec_a, t = verlet(prm)
    
    end_time = time.time()
    
    time_verlet = np.append(time_verlet, end_time - start_time)
    
    #%% euler
    
    start_time = time.time()
    
    t_vect, x_vect, v_vect = euler(prm)
    
    end_time = time.time()
    
    time_euler = np.append(time_euler, end_time - start_time)

#%% traitement des données

time_verlet = np.sum(time_verlet)

time_euler = np.sum(time_euler)

time = np.array([time_verlet, time_euler])

method = np.array(["Verlet", "Euler explicite"])
    
#%% generation de graph

plt.barh(method, time, color=['blue', 'green'])
plt.xlabel(r"Temps d'exécution [s]") 
plt.ylabel(r"Type de méthode")
plt.title(r"Somme des temps d'exécution pour $\Delta t \in [1 \times 10^{-1}, 1 \times 10^{-5}]$" +"\nen fonction des différentes méthodes utilisées")
plt.savefig("euler_verlet.png", dpi=300,bbox_inches='tight')
plt.show()
