from parameters import *
from background_fields import *
from functions_for_integral import *
import matplotlib.pyplot as plt
from sys import exit

def horizontal_line():
    print('-------------------------')
    return

def line_breaks(n=1):
    for i in range(n):
        print('')
    return

def write_parameters():
    horizontal_line()
    print('Background field parameters:')
    print('Uinf = {:.2f} km/s'.format(Uinf))
    print('Va0 = {:.2f} km/s'.format(Va0))
    print('alpha = {:.3f}, beta = {:.3f}'.format(alpha,beta))
    horizontal_line()

    line_breaks()

    horizontal_line()
    print('Grid parameters:')
    print('x1 = {:.3f}, x2 = {:.3f}'.format(x1,x2))
    print('dx = {:.4f}'.format(dx))
    horizontal_line()

    line_breaks()

    return



if __name__ == "__main__":
    write_parameters()

    # generate grid
    x_arr = np.arange(x1,x2 + dx/10, dx)
    
    # Find alfven point ----
    x_alf = find_alfven_point()
    print('Alfven point is x_alf = {:.3f} Rs'.format(x_alf))


    # find closest point in x_arr to the Alfven point
    for i in range(len(x_arr)-1):
        if (x_arr[i]-x_alf)*(x_arr[i+1]-x_alf)<=0:
            ind_alf = i 
            break 
    # print(x_arr[ind_alf],x_arr[ind_alf+1])


    # loop for different wave periods --- 
    for iome in range(len(T_wave_arr_hours)):
        horizontal_line()
        print('Calculating Twave = {:.2f} hrs......'.format(T_wave_arr_hours[iome]))

        T_wave = T_wave_arr_hours[iome] * 3600 # seconds
        omega = 2 * np.pi / T_wave

        # start integral from the Alfven point -------------
        zp_alf = 100 + 0j # km/s, this is arbitrary (both amplitude and phase)

        # calculate zm at the Alfven point.
        zm_alf = calc_zm_at_Alfven_point(x_alf,zp_alf,omega=omega)
        # print(zm_alf)

        # integrate from the Alfven point inward
        result = integrate_zp_zm(zp_alf,zm_alf,x_alf,x_arr[0],
            np.flip(x_arr[0:(ind_alf+1)]),omega=omega,x_alf=x_alf)
        if result==None:
            print('Integration inward failed!!!')
            exit()
        
        x_inward = result['x']
        zp_inward = result['zp']
        zm_inward = result['zm']

        x_inward = np.flip(x_inward[1:])
        zp_inward = np.flip(zp_inward[1:])
        zm_inward = np.flip(zm_inward[1:])

        # integrate from the Alfven point outward
        result = integrate_zp_zm(zp_alf,zm_alf,x_alf,x_arr[-1],
            x_arr[(ind_alf+1):],omega=omega,x_alf=x_alf)
        if result==None:
            print('Integration outward failed!!!')
            exit()
        
        x_outward = result['x']
        zp_outward = result['zp']
        zm_outward = result['zm']


        # concatenate inward & outward
        x_sol = np.concatenate([x_inward,x_outward])
        zp_sol = np.concatenate([zp_inward,zp_outward])
        zm_sol = np.concatenate([zm_inward,zm_outward])        

        # calculate background field
        rho_arr = np.zeros(x_sol.shape)
        U_arr = np.zeros(x_sol.shape)
        Va_arr = np.zeros(x_sol.shape)

        for i in range(len(x_sol)):
            rho = calc_rho(x_sol[i])
            Va = calc_Va(x_sol[i],rho=rho)
            U = calc_U(x_sol[i],rho=rho,Va=Va)

            rho_arr[i] = rho
            Va_arr[i] = Va 
            U_arr[i] = U 

        # save the solution and background field
        np.save('./solution_T_{:.2f}_hr.npy'.format(T_wave_arr_hours[iome]),
            np.array([x_sol,rho_arr,Va_arr,U_arr,zp_sol,zm_sol],dtype=complex))



        horizontal_line()
        line_breaks()
