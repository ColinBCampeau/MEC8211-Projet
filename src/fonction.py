"""

PROJET MEC8211
                            >FICHIER CLASSE<

AUTEUR: EROJ MOHAMMAD ISHOQ, COLIN BISSONNETTE-CAMPEAU, TRISTAN ANCEL-SÉGUIN
CRÉATION: 2 AVRIL 2024
MISE À JOUR: 6 AVRIL 2024

"""
#%%===================== IMPORTATION DES MODULES ==========================%%#

import numpy as np

from classe import *

from matplotlib import pyplot as plt

#%%===================== FONCTION SOLUTION ANALYTIQUE ==========================%%#
def sol_analytique(prm):
    
    """
    Fonction qui permet de calculer la solution analytique du système
    
    Entrées:
               
        prm -> paramètres
        
    Sortie:
        
        x_anal -> vecteur des valeurs de position pour chaque pas de temps
        
        t -> vecteur des valeurs des pas de temps
    
    """
    
    dt = prm.dt

    t = np.arange(0,prm.t_fin+dt,dt)

    A = np.sqrt((prm.w_d**2*prm.x0**2+(prm.v0+prm.zeta*prm.w_n*prm.x0)**2)/prm.w_d**2)
    phi = np.arctan(prm.w_d*prm.x0/(prm.v0+prm.zeta*prm.w_n*prm.x0))

    x_anal = A*np.exp(-prm.zeta*prm.w_n*t)*np.sin(prm.w_d*t+phi)
    
    
    return x_anal, t

#%%================= FONCTION MÉTHODE VELOCITY VERLET ====================%%#
def verlet(prm):
    
    """
    Fonction qui permet de calculer la solution numérique du système avec la 
    méthode de verlet
    
    Entrées:
               
        prm -> paramètres
        
    Sortie:
        
        vec_x -> vecteur des valeurs de position pour chaque pas de temps
        
        vec_v -> vecteur des valeurs de vitesse pour chaque pas de temps
        
        vec_a -> vecteur des valeurs d'accélération pour chaque pas de temps
        
        t -> vecteur des valeurs des pas de temps
    
    """ 
    
    v_demidt = 0
    x_dt = 0
    v_dt = 0
    a_dt = 0
    dt = prm.dt
    x0 = prm.x0
    v0 = prm.v0
    a0 = prm.a0
    
    t = np.arange(0,prm.t_fin+dt,dt)

    vec_x = np.zeros(len(t))
    vec_v = np.zeros(len(t))
    vec_a = np.zeros(len(t))
    vec_x[0] = x0
    vec_v[0] = v0
    vec_a[0] = a0


    for i in range(len(t)-1):
        
        v_demidt = v0 + 0.5*dt*a0
        x_dt = x0 + dt*v_demidt
        a_dt = -(prm.c*v_demidt+prm.k*x_dt)/prm.m
        v_dt = v_demidt + 0.5*dt*a_dt
        
        vec_x[i+1] = x_dt
        vec_v[i+1] = v_dt
        vec_a[i+1] = a_dt
        
        x0 = x_dt
        v0 = v_dt
        a0 = a_dt
    
    return vec_x, vec_v, vec_a, t

#%%================= FONCTION MÉTHODE EULER EXPLICITE ====================%%#
def euler(prm):
    
    """
    Fonction qui permet de calculer la solution numérique du système avec la 
    méthode euler explicite
    
    Entrées:
               
        prm -> paramètres
        
    Sortie:
        
        x_vect -> vecteur des valeurs de position pour chaque pas de temps
        
        v_vect -> vecteur des valeurs de vitesse pour chaque pas de temps
        
        t_vect -> vecteur des valeurs des pas de temps
    
    """ 
       
    x = prm.x0 # Position initiale
    v = prm.v0 # Vitesse initiale


    t_vect = np.arange(0,prm.t_fin+prm.dt,prm.dt) # vecteur des valeurs des pas de temps
    x_vect = np.array([x]) # Initialisation vecteur des valeurs de position pour chaque pas de temps
    v_vect = np.array([v]) # Initialisation vecteur des valeurs de vitesse pour chaque pas de temps

    for i in range(len(t_vect)-1):
        
        # calcul de la position au temps t+dt
        
        x = x + prm.dt * v 
        
        # calcul de la vitesse au temps t+dt
        
        v = v - prm.dt * (prm.c * v + prm.k * x_vect[-1]) / prm.m
        
        # Stockage des résultats dans les vecteurs respectifs
        
        x_vect = np.append(x_vect, x)
        
        v_vect = np.append(v_vect, v)

    return t_vect, x_vect, v_vect

    

#%%=========================== FONCTION ERREUR L2 ==========================%%#

def f_L2(x_num, x_anal):
    
    """
    Fonction qui permet de calculer l'erreur L2
    
    Entrés:
        
        x_num -> vecteur de la position obtenue de manière numérique
        
        x_anal -> vecteur de la position obtenue de manière analytique
    
    Sortie:
        
        L2 -> erreur L2
    
    """
    
    L2 = (np.sum((x_num-x_anal)**2)/(len(x_num)))**0.5
    
    return L2
    



#%%=========================== FONCTION FRÉQUENCE ==========================%%#

def frequence(vec_x):
    
    """
    Fonction qui permet de calculer la fréquence de l'oscillation
    
    Entrés:
        
        vec_x -> vecteur de la position obtenue de manière numérique
        
    Sortie:
        
        freq -> fréquence d'oscillation
    
    """
    
    c = 0 
    
    for i in range(1,len(vec_x)-1):
        if vec_x[i] > vec_x[i-1] and vec_x[i] > vec_x[i+1]:
             
            if c == 1:
                freq = 1/(i*prm.dt-val_prec)
                #print('Fréquence :',freq, 'Hz')
            val_prec = i*prm.dt
            c += 1
    
    return freq