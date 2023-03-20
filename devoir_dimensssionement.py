import numpy as np
import math as math
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from sympy.solvers import solve
from sympy import Symbol

######################
#   VALEURS FIXES    #
######################

tau     =   10         # valeur taux compression [-]
D       =   0.08       # valeur alesage   [m]
C       =   0.09       # valeur course [m]
L       =   0.18       # valeur longueur bielle@ [m]
mpiston =   0.8        # valeur masse piston@ [kg]
mbielle =   0.9        # valeur masse bielle  [kg]
Q_kg    =   2800000    # valeur chaleur emise par fuel par kg de melange admis@ [J/kg_inlet gas]
E       =   210e9      # Module de Young de l'acier [pa]. Source : cours de Mécanique des Structures.

global zzz
zzz = 0
rad_degre = math.pi/180

#Valeurs liées au volume

R = C / 2              # Longueur de manivelle
V_critique = (( math.pi * D**2) / 4 ) * (2 * R) # OK
masse_air = 0.029 * V_critique *1e5/(8.31451*303.15) # En [pa]
beta = L / R           # rapport entre la longueur de la bielle et la longueur de la manivelle
V_min = (1 / ( tau - 1)) * V_critique
V_max = (tau / (tau -1)) * V_critique


def vol (theta):
    
    """""""""""""""""
    Calcul du V_output en fonction de thêta
    
    """""""""""""""""
    V_output = np.zeros(theta.size);

    for i in range(len(theta)):
        V_output[i] = V_critique/2 * (1 - math.cos(theta[i]*rad_degre) + beta - np.sqrt(beta**2 - math.sin(theta[i]*rad_degre)**2)) + V_min
    
    if(theta.size != 1): plot_vol(V_output)
    return V_output

#graphe volume
def plot_vol(V_output):
    plt.plot(theta, V_output/V_max)
    
    # Axes
    fs_text = 16 # Taille du texte
    plt.xlabel("$Thêta$ [degré]", fontsize=fs_text)
    plt.ylabel("$V_{output}$ [m^3]", fontsize=fs_text) #FAUT PAS MARQUER DIVISER PAR V_max?
    
    # Titre
    plt.title("Volume/V_max en fonction de l'angle de vilebrequin", fontsize=fs_text)
    plt.show();
    return;

def q(theta):
    """""""""""""""""
    Calcul du Q_output en fonction de thêta 
    """""""""""""""""
    Q_output = np.zeros(theta.size)
    Q_max = float('-inf')
    Q_min = float('inf')
    for i in range (len(theta)):
        if (theta[i] < thetaC or theta [i] > thetaC + deltaThetaC):
            Q_output[i] = 0
        
        else:
            fraction = (theta[i] + thetaC)/deltaThetaC # Pas besoin de mettre en radiant car fraction.
            Q_output[i] = (Q_kg*s* masse_air)/ 2 * (1 - math.cos(math.pi * fraction))
        
        if (Q_output[i] > Q_max) :
            Q_max = Q_output[i] 
            print("Q est parti trop haut en " + str(i))
        if (Q_output[i] < Q_min) :
            Q_min = Q_output[i]
    
    return Q_output

#Graphe de Q_output en fct de theta
def plot_Q_output(theta, Q_output):
    plt.plot(theta, Q_output)
    
    # Axes
    fs_text = 16 # Taille du texte
    plt.xlabel("$Thèta$ [degré]", fontsize=fs_text)
    plt.ylabel("$Q_{output}$ [J]", fontsize=fs_text)
    
    # Titre
    plt.title("Q_max = " + str(Q_max) + "    Q_min = " + str(Q_min))
    plt.suptitle("Chaleur en fonction de l'angle de vilebrequin", fontsize=fs_text)
    plt.show()
    return

def dVol(theta, theta_rad, racine) :
    return V_critique/2*(math.sin(theta_rad) + (math.sin(theta_rad)*math.cos(theta_rad))/(racine))

def dQ(theta, thetaC, deltaThetaC):
    fraction = (theta-thetaC)/deltaThetaC
    return (math.pi*Q_kg*masse_air*s)/(2*deltaThetaC) * (math.sin(math.pi*fraction))

def data(theta, V, dp, dV, dQ):
    if(dp > 0): strin = ' '
    else: strin = '';
    
    if(int(theta*1000)%13 == 0):
        print('pour theta = {:7.2f}, taux du volume : {:3.2f}, dp/dtheta : {}{:E}, dQ/dtheta : {:7.2f}, dV/dtheta : {:5.3f}'.format(theta, V/V_max, strin, dp[0], dQ, dV/V_max))

def pression(theta, s, thetaC, deltaThetaC):
    """""""""""""""""
    Calcul du p_output en fonction de thêta 
    """""""""""""""""
    
    #dérivée de p en fonction de l'angle du coup aussi
    def model(p,theta):
        """
        ici p et theta sont des floats et plus des array!!!
        """
        theta_rad = theta*rad_degre #math.radians(theta)
        racine = np.sqrt(beta*beta - (math.sin(theta_rad)*math.sin(theta_rad)))
        
        theta_array = np.ones(1)*theta
        V_formule = vol(theta_array)[0];
        #V_formule = V_critique/2 * (1 - math.cos(theta_rad) + beta - racine) + V_min # Je peux pas prendre vol(theta) car renvoie un array et j'ai besoin d'un float.
        
        #print('pour theta = ' + str(theta) + ", le taux du volume est de : " + str(V_formule/V_max))
        
        dVdtheta = dVol(theta, theta_rad, racine)
        dQdtheta = dQ(theta, thetaC, deltaThetaC)
        
        dpdtheta = (- 1.3*p * dVdtheta + 0.3*dQdtheta)/V_formule
        
        """
        if (thetaC <= theta and theta < thetaC + deltaThetaC) : # Si dans combustion
            dpdtheta = - 1.3 *p/V_formule * dVdtheta + (1.3 - 1)/V_formule * dQdtheta
        else :
            if(int(theta*1000)%13 == 0): print('True')
            dpdtheta = - 1.3 *p/V_formule * dVdtheta
        """
        
        data(theta, V_formule, dpdtheta, dVdtheta, dQdtheta)
        #if(zzz != int(theta)):
        #print('pour theta = ' + str(theta) + ", taux du volume : " + str(V_formule/V_max), "dp/dtheta : " + str(dpdtheta))
        #    zzz = int(theta)
        return dpdtheta[0]
 
    p_output = odeint(model,s*1e5,theta)
    plt.plot(theta,p_output)
    plt.title("Pression")
    plt.show()
    return p_output
     

def F_tete(theta, p_output, w) : 
    """""""""""""""""
    Calcul du F_tete_output en fonction de thêta 
    """""""""""""""""
    F_tete_output = np.zeros(theta.size)
    for k in range(len(theta)) : 
        F_tete_output[k] = -(math.pi*D**2)/4 *p_output[k] + (mpiston + mbielle)*R*w**2*math.cos(math.radians(theta[k]))
     
            
    plt.plot(theta, F_tete_output)
    plt.title("F_tete_output")
    plt.show()

    return F_tete_output


def F_pied(theta, p_output, w) : 
    """""""""""""""""
    Calcul du F_pied_output en fonction de thêta 
    """""""""""""""""
    F_pied_output = np.zeros(theta.size)
    for j in range(len(theta)) : 
        F_pied_output[j] = (math.pi*D**2)/4 *p_output[j] - mpiston*R*w**2*math.cos(math.radians(theta[j]))
        
    plt.plot(theta, F_pied_output)
    plt.title("F_pied_output")
    plt.show()

    return F_pied_output


def calcul_t() : 
    """""""""""""""""
    Calcul du t 
    """""""""""""""""
    t = Symbol('t')
    Ixx = 419/12 * t # Calculé à la main. Ok avec vérif du cours de Mécanique des Structures.
    Iyy = 131/12 * t

    Aire = 11*t**2

    Imax = np.max(Ixx, Iyy)
    if (Imax == Ixx) : # Pas sure !!! Voir notes ! 
        x = 1
    else : 
        x = 0.5

    F_critique = math.pi*math.pi*E*Imax/(L**2 * x**2)
    Fmax = 0 # A trouver de foorce dans tete et pied de bielle !! Prendre la force max !

    # Solve : fonction permettant de résoudre une équation. 1er argument : équation à résoudre égale à 0.
    # 2eme argument : variable.

    tab = solve(Fmax - F_critique, t)

    t = min(tab)
    print("t = " + str(t))
    return t
        

def myfunc(rpm, s, theta, thetaC, deltaThetaC):
    """ 
    dimBielle dimensionnement d'une bielle
    dimBielle(rpm, s, theta, thetaC, deltaThetaC) calcules les données thermodynamiques
    et les forces d'un système Piston-Bielle-Vilebrequin afin de dimensionner la section
    d'une bielle.
    INPUTS :
    rpm : vitesse angulaire du moteur [rotation per minute]
    s : surcharge du moteur [-]
    theta : angle auxquels renvoyer les données [°]
    thetaC : angle d'allumage [°]
    deltaThetaC : durée de la combustion (en angle) [°] 
    theta : angle auxquels renvoyer les données [°]
    (Les angles sont donnés entre 0 et 720°)
    OUTPUTS :
    t : section de la bielle [m]
    V(theta) : Volume de la chambre de combustion en fonction de theta [m3]
    Q(theta) : chaleur dégagée par la combustion en fonction de theta [J]
    F_pied(theta) : [N]
    F_tete(theta) : [N]
    F_inertie(theta) : [N]
    (une force est positive si dirigée vers le haut).
    """

    w = rpm/60*2*math.pi # Transformation de revolution/minute en rad/s.    

    V_output = vol (theta)
    Q_output = q(theta)
    p_output = pression(theta, s, thetaC, deltaThetaC)
    F_tete_output = F_tete(theta, p_output, w)
    F_pied_output = F_pied(theta, p_output, w)
    t = calcul_t()
    print(t)
    
    #V_output = np.ones_like(theta)
    #Q_output      = 2*np.ones_like(theta);
    #F_pied_output = 3*np.ones_like(theta);
    #F_tete_output = 4*np.ones_like(theta);
    #p_output      = 5*np.ones_like(theta);
    #t             = 6;

    return ( F_pied_output, F_tete_output, p_output, t  )



if __name__ == '__main__':
    print("""À faire:
    # Faire fonction pour dV et dQ pour avoir plus facile.
    # Faire docstrings.
    # Questions pour t.
    # Problème sur p --> réglé avec décomposition des fonctions ?

Remarques:
    # Faudrait pas mettre dans l'axe du graphe V_output/V_max, pour éviter de porter à confusion que ce serait 1 m^3?
    """)
    
    rpm = 2000
    s = 1.5
    Pin = s*1E5
    thetaC = 25
    deltaThetaC = 50
    theta = np.linspace(-180., 180., 100)

    #print(myfunc(2555,1.9,np.linspace(-1*180,180,10000), 26,43))
    print(myfunc(rpm, s, theta, thetaC, deltaThetaC))
    quit()
    
    

