from parameters import *
from background_fields import *
from scipy.integrate import solve_ivp

def calc_zm_at_Alfven_point(x_alf,zp_alf,omega=None):
    x = x_alf
    rho = calc_rho(x)
    Va = calc_Va(x,rho=rho)
    U = calc_U(x,rho=rho,Va=Va)
    
    drho = calc_drho(x, rho=rho)
    dVa = calc_dVa(x, rho=rho, drho=drho,Va=Va)
    dU = calc_dU(x,rho=rho,drho=drho,Va=Va,dVa=dVa,U=U)


    r = x * Rsun # Rs to km
    dU_dr = dU/Rsun # so that dU/dr will have the unit s^{-1}
    dVa_dr = dVa/Rsun 

    #  1/(2r^2) * d[r^2 * (Va+U/2)]/dr
    term1 = (Va+U/2)/r + (dVa_dr + dU_dr/2)/2

    top = (term1 - (U + Va)/r) * zp_alf
    bottom = term1 -1j*omega 

    return top/bottom


def calc_dzp_dr(x,zp,zm,omega=None,x_alf=None):
    rho = calc_rho(x)
    Va = calc_Va(x,rho=rho)
    U = calc_U(x,rho=rho,Va=Va)
    
    drho = calc_drho(x, rho=rho)
    dVa = calc_dVa(x, rho=rho, drho=drho,Va=Va)
    dU = calc_dU(x,rho=rho,drho=drho,Va=Va,dVa=dVa,U=U)


    r = x * Rsun # Rs to km
    dU_dr = dU/Rsun # so that dU/dr will have the unit s^{-1}
    dVa_dr = dVa/Rsun 

    # 1/(2r^2) * d[r^2 * (Va-U/2)]/dr
    term1 = (Va-U/2)/r + (dVa_dr - dU_dr/2)/2

    dzp_dr = (1j * omega  * zp - (U - Va)/r * zm 
        - term1 * (zm-zp))/(U+Va)

    return dzp_dr


def calc_dzm_dr(x,zp,zm,dzp_dr,omega=None,x_alf=None):
    if x_alf==None:
        print('x_alf is needed for calc_dzm_dr!!!')
        return None

    rho = calc_rho(x)
    Va = calc_Va(x,rho=rho)
    U = calc_U(x,rho=rho,Va=Va)
    
    drho = calc_drho(x, rho=rho)
    dVa = calc_dVa(x, rho=rho, drho=drho,Va=Va)
    dU = calc_dU(x,rho=rho,drho=drho,Va=Va,dVa=dVa,U=U)

    r = x * Rsun # Rs to km
    dU_dr = dU/Rsun # so that dU/dr will have the unit s^{-1}
    dVa_dr = dVa/Rsun 

    # 1/(2r^2) * d[r^2 * (Va+U/2)]/dr
    term1 = (Va+U/2)/r + (dVa_dr + dU_dr/2)/2

    if np.abs(x-x_alf)<=0.05:
        # very close to the Alfven point
        ddrho = calc_ddrho(x,rho=rho,drho=drho)

        ddVa = calc_ddVa(x,rho=rho,drho=drho,ddrho=ddrho,Va=Va,dVa=dVa)
        ddU = calc_ddU(x,rho=rho,Va=Va,U=U,drho=drho,dVa=dVa,dU=dU,ddrho=ddrho,ddVa=ddVa)

        ddVa_dr2 = ddVa/Rsun/Rsun # so the unit is s^{-2}
        ddU_dr2 = ddU/Rsun/Rsun

        top = -((dU_dr+dVa_dr)/r - (U+Va)/r/r)*zp + \
            (dVa_dr/2 + dU_dr/4 - U/r/2) * dzp_dr - \
            (ddVa_dr2/2 + ddU_dr2/4 + (dVa_dr + dU_dr/2)/r 
            - (Va+U/2)/r/r) * (zm-zp)
        bottom = (dU_dr - dVa_dr) - 1j*omega + term1
        dzm_dr = top/bottom
    else:
        dzm_dr = (1j * omega * zm - (U + Va)/r * zp 
            - term1 * (zm-zp))/(U - Va)
    
    return dzm_dr


def calc_deriv(x,zpzm,omega=None,x_alf=None):
    if x_alf==None:
        print('x_alf is needed for calc_deriv!!!')
        return None

    zp = zpzm[0]
    zm = zpzm[1]

    dzp_dr = calc_dzp_dr(x,zp,zm,omega=omega,x_alf=x_alf)
    dzm_dr = calc_dzm_dr(x,zp,zm,dzp_dr,omega=omega,x_alf=x_alf)

    dzp_dx = dzp_dr * Rsun
    dzm_dx = dzm_dr * Rsun

    return [dzp_dx,dzm_dx]


def integrate_zp_zm(zp0,zm0,x_start,x_end,x_output,omega=None,x_alf=None):
    if x_alf==None:
        print('x_alf is needed for integrate_zp_zm!!!')
        return None

    result = solve_ivp(calc_deriv,[x_start,x_end],[zp0,zm0],
        args=(omega,x_alf),t_eval=x_output)

    if result.status != 0:
        print('Integration Failed!!!')
        return None


    x_result = result.t 
    zpzm_result = result.y

    zp_arr = zpzm_result[0,:]
    zm_arr = zpzm_result[1,:]

    return {'x':x_result,'zp':zp_arr,'zm':zm_arr}
