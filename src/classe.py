"""

PROJET MEC8211
                            >FICHIER CLASSE<

AUTEUR: EROJ MOHAMMAD ISHOQ, COLIN BISSONNETTE-CAMPEAU, TRISTAN ANCEL-SÉGUIN
CRÉATION: 2 AVRIL 2024
MISE À JOUR: 14 AVRIL 2024

"""

#%%===================== IMPORTATION DES MODULES ==========================%%#

import numpy as np

#%%===================== Classe paramètre ==========================%%#
class prm:
    
    m = 1.9 # Masse [kg] 
    
    c = 11.05 # Amortissement [N*s/m]
    
    k = 1666 # Constante de rappel [N/m]
    
    w_n = np.sqrt(k/m) # Fréquence naturelle [Hz]
    
    zeta = c/(2*w_n*m) # Coefficient d'amortissement
    
    # c = zeta*2*m*w_n
    
    w_d = w_n*np.sqrt(1-zeta**2) # Fréquence amortie [Hz]
    
    x0 = 0.1 # Position initiale [m]
    
    v0 = 10 # Vitesse initiale [m/s]
    
    a0 = -1/m*(c*v0+k*x0) # Accélération initiale [m/s^2]
    
    dt = 1e-5 # Pas de temps [s]
    
    t_fin = 1 # Temps simulé [s]
    
    delta_m = 0.01 # Incertitude expérimentale sur la masse [kg] #0.05
    
    delta_k = k*0.1 # Incertitude expérimentale sur la constante de rappel [N/m] #0.5
    
    err_t_reaction = 0.273 # Erreur commise par le temps de réaction humain sur le temps[s]
    
    err_t_lecture = 0.1 # Erreur commise par la lecture du chronomètre [s]
    
    # err_mesures_pos = 0.5 # Erreur commise sur la lecture de l'amplitude d'un cycle [mm]
    
    

    