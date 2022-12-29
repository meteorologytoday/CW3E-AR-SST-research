import numpy as np

g0 = 9.80616     # m / s**2      copied from models/csm_share/shr/shr_const_mod.F90

# Linear TS and buoyancy parameterization
alpha_T = 2e-4
alpha_S = 8e-4
T_ref = 0.0
S_ref = 35.0

def TS2b(T, S):

    dT = T - T_ref
    dS = S - S_ref

    return g0 * (alpha_T * dT - alpha_S * dS) 

