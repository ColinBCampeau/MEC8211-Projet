"""

OUTIL NUMÉRIQUE DE SIMULATION DE LA DIFFUSION DU SEL DANS UN PILIER EN BÉTON POREUX

                        >FICHIER ANALYSE CONVERGENCE - COMSOL<

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
try:
    from fonction import *
    from classe import *
except:
    pass



#%%============= ANALYSE CONVERGENCE - TEMPORELLE ===============%%#

delta_t = np.array([1e-7,1e-6,1e-5,1e-4,1e-3,1e-2, 1e-1])

vec_l2_t = np.array([])

for i in range(0, len(delta_t)):
    
    prm.dt = delta_t[i]
    
    x_anal, t = sol_analytique(prm)

    vec_x, vec_v, vec_a, t = verlet(prm)

    L2_verlet = f_L2(vec_x, x_anal)

    vec_l2_t = np.append(vec_l2_t, L2_verlet)

#vec_l2_t = np.array([1.128497e-7,1.128515e-6,1.1287e-5,1.1305e-4,1.14912e-3,1.3517e-2])

# Ajuster une loi de puissance à toutes les valeurs (en utilisant np.polyfit avec logarithmes)
coefficients = np.polyfit(np.log(delta_t[0:-4]), np.log(vec_l2_t[0:-4]), 1)
exponent = coefficients[0]

# Fonction de régression en termes de logarithmes
fit_function_log = lambda x: exponent * x + coefficients[1]

# Fonction de régression en termes originaux
fit_function = lambda x: np.exp(fit_function_log(np.log(x)))

# Extrapoler la valeur prédite pour la dernière valeur de h_values
extrapolated_value = fit_function(delta_t[2])

# Tracer le graphique en échelle log-log avec des points et la courbe de régression extrapolée
plt.figure(figsize=(8, 6))
plt.scatter(delta_t, vec_l2_t, marker='o', color='b', label='Données numériques obtenues')
plt.plot(delta_t, fit_function(delta_t), linestyle='--', color='r', label='Régression en loi de puissance')

# Marquer la valeur extrapolée
#plt.scatter(h_values[-1], extrapolated_value, marker='x', color='g', label='Extrapolation')

# Ajouter des étiquettes et un titre au graphique
plt.title('Convergence d\'ordre 1\n de l\'erreur $L_2$ en fonction de $Δt$',
          fontsize=14, fontweight='bold', y=1.02)  # Le paramètre y règle la position verticale du titre

plt.xlabel('Taille d\'intervalle $Δt$ (s)', fontsize=12, fontweight='bold')  
plt.ylabel('Erreur $L_2$ (mol/m$^3$)', fontsize=12, fontweight='bold')

# Rendre les axes plus gras
plt.gca().spines['bottom'].set_linewidth(2)
plt.gca().spines['left'].set_linewidth(2)
plt.gca().spines['right'].set_linewidth(2)
plt.gca().spines['top'].set_linewidth(2)

# Placer les marques de coche à l'intérieur et les rendre un peu plus longues
plt.tick_params(width=2, which='both', direction='in', top=True, right=True, length=6)

# Afficher l'équation de la régression en loi de puissance
equation_text = f'$L_2 = {np.exp(coefficients[1]):.4f} \\times Δt^{{{exponent:.4f}}}$'
equation_text_obj = plt.text(0.05, 0.05, equation_text, fontsize=12, transform=plt.gca().transAxes, color='k')

# Déplacer la zone de texte
equation_text_obj.set_position((0.5, 0.4))


# Afficher le graphique
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.savefig("L2_verlet_dt.png", dpi=300,bbox_inches='tight')
plt.show()
