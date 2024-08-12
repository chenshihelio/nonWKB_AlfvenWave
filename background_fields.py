import numpy as np 
from numpy import exp,sqrt
from parameters import *
from scipy.optimize import root_scalar

def calc_rho(x):
    # Calculate the rho/rho0 according to Equi (6) of Vallietal&Velli, 1991
    rho = exp(-alpha*(1. - 1./x))/(1. + beta*(x - 1.))**2.
    return rho


def calc_Va(x,rho = None):
    # Calculate the Va according to Equi (6) of Vallietal&Velli, 1991

    if rho == None:
        rho = calc_rho(x)

    Va = Va0 * (1./x)**2. * sqrt(1./rho)
    return Va

def calc_U(x,rho=None,Va=None):
    #Calculate the U according to Equi (6) of Vallietal&Velli, 1991
    #Notice that, in the equation for U, it should be exp(-alpha) rather than exp(-alpha/2)

    if rho==None:
        rho = calc_rho(x)

    if Va==None:
        Va = calc_Va(x,rho=rho)

    # U = Uinf/beta**2. * EXP(-alpha/2.) * (x*Va_temp/Va0)**2.  
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    U = Uinf/beta**2. * exp(-alpha) * (x*Va/Va0)**2.  
    return U



def calc_drho(x,rho=None):
    # Calculate the derivative of rho
    if rho==None:
        rho = calc_rho(x)
    
    drho = -rho*(alpha/x**2. + 2* beta/(1.+beta*(x-1.) ))
    return drho

def calc_dVa(x,rho=None,drho=None,Va=None):
    # Calculate the derivative of va
   
    if rho==None:
        rho = calc_rho(x)
    
    if drho==None:
        drho = calc_drho(x,rho=rho)

    if Va==None:
        Va = calc_Va(x,rho=rho)

    dVa = -Va*(2./x + drho/(2.*rho))
    return dVa

def calc_dU(x,rho = None, drho=None, Va=None, U=None,dVa=None):
    # Calculate the derivative of U
    if rho==None:
        rho = calc_rho(x)

    if drho==None:
        drho = calc_drho(x,rho=rho)

    if Va==None:
        Va = calc_Va(x,rho=rho)
        
    if U==None:
        U = calc_U(x,rho=rho,Va=Va)

    if dVa == None:
        dVa = calc_dVa(x,rho=rho,drho=drho,Va=Va)

    dU = 2* U * (1./x + dVa/Va)
    return dU


######### 2nd order derivatives of the background fields #######
# Used in calculating dz-/dr at the Alfven point (L'Hopital's rule)
def calc_ddrho(x,rho=None,drho=None):
    # Calculate the second-order derivative of rho
    if rho==None:
        rho = calc_rho(x)

    if drho==None:
        drho = calc_drho(x,rho=rho)

    ddrho = -(drho*(alpha/x**2. + 2.*beta/(1.+beta*(x-1.))) + 
            rho* (-2.*alpha/x**3. - 2. * beta**2. / (1.+beta*(x-1.))**2. ))

    return ddrho

def calc_ddVa(x,rho=None,drho=None,ddrho=None,Va=None,dVa=None):
    # Calculate the second-order derivative of rho

    if rho==None: 
        rho = calc_rho(x)

    if drho == None:
        drho = calc_drho(x,rho=rho)
    
    if ddrho == None:
        ddrho = calc_ddrho(x,rho=rho,drho=drho)

    if Va == None:
        Va = calc_Va(x,rho=rho)

    if dVa == None:
        dVa = calc_dVa(x,rho=rho,drho=drho,Va=Va)

    ddVa = -(dVa*(2./x + drho/(2.*rho)) + Va*(-2./x**2. +
            ddrho/(2.*rho) - (drho/rho)**2. / 2.))
    return ddVa


def calc_ddU(x,rho=None,Va=None,U=None,drho=None,dVa=None,
    dU=None,ddrho=None,ddVa=None):
    # Calculate the second-order derivative of rho
    
    if rho==None:
        rho = calc_rho(x)

    if Va==None:
        Va = calc_Va(x,rho=rho)

    if U==None:
        U = calc_U(x,rho=rho,Va=Va)


    if drho==None:
        drho = calc_drho(x,rho=rho)

    if dVa == None:
        dVa = calc_dVa(x,rho=rho,drho=drho,Va=Va)
    
    if dU==None:
        dU = calc_dU(x,rho=rho,drho=drho,Va=Va,U=U,dVa=dVa)

    if ddrho==None:
        ddrho = calc_ddrho(x,rho=rho,drho=drho)

    if ddVa==None:
        ddVa = calc_ddVa(x,rho=rho,drho=drho,ddrho=ddrho,Va=Va,dVa=dVa)


    ddU =2. * (dU*(1./x + dVa/Va) + U*(-1./x**2. +ddVa/Va - (dVa/Va)**2.))

    return ddU


def diff_between_U_Va(x):
    rho = calc_rho(x)

    Va = calc_Va(x,rho=rho)
    U = calc_U(x,rho=rho,Va=Va)

    return U-Va


def find_alfven_point(x1=x1,x2=x2):
    sol = root_scalar(diff_between_U_Va,bracket=[x1,x2],method='bisect')

    if sol.converged:
        return sol.root
    else:
        print('Fail to find Alfven point!')
        return None


if __name__ == "__main__":
    print('Yuo are running background_fields.py: Module file, nothing will happen!')