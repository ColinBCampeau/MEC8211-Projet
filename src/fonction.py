"""

PROJET MEC8211
                            >FICHIER FONCTION<

AUTEUR: EROJ MOHAMMAD ISHOQ, COLIN BISSONNETTE-CAMPEAU, TRISTAN ANCEL-SÉGUIN
CRÉATION: 2 AVRIL 2024
MISE À JOUR: 14 AVRIL 2024

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
    # print(len(t))

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



#%%=========================== FONCTION MHS ==========================%%#

def mhs(n, val_min, val_max):
    
    """
    Fonction qui permet de générer un vecteur de n valeurs aléatoire comprises entre deux valeurs
    
    Entrés:
        
        n -> taille de l'échantillon à générer
        
        val_min -> valeur minimale de l'échantillon
        
        val_max -> valeur maximale de l'échantillon
        
    Sortie:
        
        vec_rand -> vecteur de taille n contenant un échantillon aléatoire de données entre val_min et val_max
    
    """
    
    vec_rand = np.random.ran(n)*(val_max-val_min) + val_min
    
    return vec_rand



#%%=========================== FONCTION LHS ==========================%%#

def lhs(n):
    
    """
    Fonction qui permet de générer un vecteur de n valeurs aléatoire comprises entre deux valeurs
    
    Entrés:
        
        n -> taille de l'échantillon à générer
        
    Sortie:
        
        points -> vecteur de taille n contenant un échantillon aléatoire de données entre val_min et val_max
    
    """
    
    limit_min = np.arange (0, n) / n
    limit_max = np.arange (1, n+1) / n
    points = np.random.uniform(low=limit_min, high=limit_max, size=[2,n]).T
    np.random.shuffle(points[:, 1])
    return points




#%%=========================== FONCTION PROPAGATION ==========================%%#

def incert_input(n, prm):
    
    """
    Fonction qui permet de calculer l'incertitude des paramètres d'entré
    
    Entrés:
        
        n -> taille des échantillons épistémiques à générer
        
        prm -> paramètres
        
    Sortie:
        
        u_input -> incertitude des paramètres d'entré
    
    """
    vec_freq = np.zeros(n)
    pts_lhs = lhs(n)
    vec_m = pts_lhs[:,0]*prm.delta_m*2+(prm.m-prm.delta_m)
    vec_k = pts_lhs[:,1]*prm.delta_k*2+(prm.k-prm.delta_k)
    # print(vec_m)
    # print(vec_k)
    
    plt.figure
    plt.scatter(vec_m, vec_k, c='r')
    plt.xlabel('m [kg]')
    plt.ylabel('k [N/m]')
    plt.title('Paires de variables aléatoires')    
    plt.show()
    
    for i in range(len(vec_freq)):
        prm.m = vec_m[i]
        prm.k = vec_k[i]
        vec_x, vec_v, vec_a, t = verlet(prm)
        freq = frequence(vec_x)
        vec_freq[i] = freq
    
    S_bar = sum(vec_freq)/n
    
    sum_LHS = 0
    for i in range(len(vec_freq)):
        sum_LHS += (vec_freq[i]-S_bar)**2
    
    u_input = np.sqrt(sum_LHS/(n-1))
    S = np.mean(vec_freq)
    
    return u_input,S



#%%=========================== FONCTION U_NUM ==========================%%#

def incert_num(r, prm):
    
    """
    Fonction qui permet de calculer l'incertitude numérique
    
    Entrés:
        
        r -> taux de raffinement
        
        prm -> paramètres
        
    Sortie:
        
        u_num -> incertitude numérique
    
    """
    vec_freq = np.zeros(3)
    dti = prm.dt
    
    for i in range(3):
        dt = prm.dt
        # print(dt)
        vec_x, vec_v, vec_a, t = verlet(prm)
        freq = frequence(vec_x)
        # print(freq)
        vec_freq[i] = freq
        prm.dt = prm.dt*r
    
    prm.dt = dti    
    
    p_hat = np.log((abs(vec_freq[2]-vec_freq[1]))/(abs(vec_freq[1]-vec_freq[0])))/np.log(r)
    # print(p_hat)
    val_p = abs(p_hat-1)
    
    if val_p > 0.1:
        Fs = 3
        p = min(max(0.5,p_hat),1)
    else:
        Fs = 1.25
        p = 1
    
    GCI = Fs/(r**p-1)*abs(vec_freq[1]-vec_freq[0])
    # print(GCI)
    u_num = GCI/2
    
    return u_num


#%%=========================== FONCTION U_D ==========================%%#

def incert_D(prm):
    
    """
    Fonction qui permet de calculer l'incertitude expériementale
    
    Entrés:
        
        prm -> paramètres
        
    Sortie:
        
        u_D -> incertitude expérimentale
    
    """
    
    bi1 = prm.err_t_reaction
    bi2 = prm.err_t_lecture
    
    bi = np.sqrt(bi1*bi1+bi2*bi2)
    
    drdxi = (5/(1.2+bi)-5/(1.2-bi))/(2*bi)
    
    br = abs(drdxi*bi)
    
    u_D = br
    
    return u_D