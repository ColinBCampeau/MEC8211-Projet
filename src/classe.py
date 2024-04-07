"""

PROJET MEC8211
                            >FICHIER CLASSE<

AUTEUR: EROJ MOHAMMAD ISHOQ, COLIN BISSONNETTE-CAMPEAU, TRISTAN ANCEL-SÉGUIN
CRÉATION: 2 AVRIL 2024
MISE À JOUR: 6 AVRIL 2024

"""

#%%===================== IMPORTATION DES MODULES ==========================%%#

import numpy as np

#%%===================== Classe paramètre ==========================%%#
class prm:
    
    m = 0.1 # Masse [kg] 
    
    c = 1 # Amortissement [N*s/m]
    
    k = 40 # Constante de rappel [N/m]
    
    w_n = np.sqrt(k/m) # Fréquence naturelle [Hz]
    
    zeta = c/(2*w_n*m) # Coefficient d'amortissement
    
    w_d = w_n*np.sqrt(1-zeta**2) # Fréquence amortie [Hz]
    
    x0 = 0.1 # Position initiale [m]
    
    v0 = 10 # Vitesse initiale [m/s]
    
    a0 = -1/m*(c*v0+k*x0) # Accélération initiale [m/s^2]
    
    dt = 1e-5 # Pas de temps [s]
    
    t_fin = 1 # Temps simulé [s]

    