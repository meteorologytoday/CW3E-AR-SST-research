import numpy as np


rho0 = 1026.0
cp = 3996.0
alpha = 2e-4
g0 = 9.81

NK_m = 0.5
NK_n = 0.2


def calDelta(h, U10, wb_prime, zeta, F_sol):
    
    u_star = calu_star(U10)
    B = calB(h, wb_prime, zeta, F_sol)

    Delta = 2 * NK_m * u_star**3 - 0.5 * ( (1 - NK_n) * abs(B) + (1 + NK_n) * B )

    return Delta

def calu_star(U10):
 
    # Wu, J. (1982). Wind‐stress coefficients over sea surface from breeze to 
    # hurricane. Journal of Geophysical Research: Oceans, 87(C12), 9704-9706.
    u_star = U10 * np.sqrt( (0.8 + 0.065 * U10) * 1e-3 ) 

    return u_star
   

def calB(h, wb_prime, zeta, F_sol):

    B = - wb_prime + alpha * g0 / (rho0 * cp) * F_sol * (1 + calI(- h, zeta) - 2 / h * calIntI(- h, zeta))
    
    return B


def calI(z, zeta):
    return np.exp(z/zeta)

def calIntI(z, zeta):
    
    return zeta * (1.0 - np.exp(z / zeta))
