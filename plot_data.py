from parameters import *
from background_fields import *
import matplotlib.pyplot as plt 
from matplotlib.gridspec import GridSpec

if __name__ == "__main__":

    # fig = plt.figure()
    # sub = fig.add_subplot(111)

    x_alf = find_alfven_point()

    for iome in range(len(T_wave_arr_hours)):
        T_wave_hours = T_wave_arr_hours[iome]
        data = np.load('./solution_T_{:.2f}_hr.npy'.format(T_wave_hours))

        x_arr = np.real(data[0])
        rho_arr = np.real(data[1])
        Va_arr = np.real(data[2])
        U_arr = np.real(data[3])

        zp_arr = data[4]
        zm_arr = data[5]


        # calculate energies 
        Ep = (np.abs(zp_arr))**2
        Em = (np.abs(zm_arr))**2

        sigma_c = (Ep - Em)/(Ep + Em)


        va_fluc = (zm_arr - zp_arr)/2
        u_fluc = (zp_arr + zm_arr)/2

        Ek = (np.abs(u_fluc))**2
        Eb = (np.abs(va_fluc))**2

        sigma_r = (Ek-Eb)/(Ek+Eb)

        
        # calculate Mach number
        Mach_a = U_arr/Va_arr # Alfven mach number
        Mach_a_scaling = Mach_a/(Mach_a+1)**2



        # sub.plot(x_arr, sigma_c, label=r'$T={:.2f}$ hr'.format(T_wave_hours),
        #     color=plt.cm.plasma(T_wave_hours / 105))

        # plot the data
        fig = plt.figure(figsize=[8,6])
        gs = GridSpec(4,1)

        xlim = [0,100]

        sub = fig.add_subplot(gs[0,0])
        sub.plot(x_arr,U_arr,color='C0',label=r'$U$')
        sub.plot(x_arr,Va_arr,color='C1',label=r'$V_a$')
        sub.legend()
        sub.axvline(x=x_alf,color='red',ls='--')
        sub.set_xlim(xlim)
        sub.set_xticklabels([])
        sub.tick_params(axis='x',which='both',direction='in')

        sub = fig.add_subplot(gs[1,0])
        sub.plot(x_arr,sigma_c,label=r'$\sigma_c$')
        sub.plot(x_arr,sigma_r,label=r'$\sigma_r$')
        sub.legend()
        sub.axvline(x_alf,color='red',ls='--')
        sub.axhline(y=0,color='grey',ls='--')
        sub.set_xlim(xlim)
        sub.set_xticklabels([])
        sub.tick_params(axis='x',which='both',direction='in')
    
        sub = fig.add_subplot(gs[2,0])
        sub.plot(x_arr,Ep,label=r'${z_+}^2$')
        sub.plot(x_arr,Em,label=r'${z_-}^2$')
        sub.legend(loc='upper right')
        sub.axvline(x_alf,color='red',ls='--')
        subt = sub.twinx()
        subt.plot(x_arr,Mach_a_scaling,ls='--',color='k',label=r'$M_a/(M_a+1)^2$')
        subt.legend(loc='lower right')
        sub.set_xlim(xlim)
        sub.set_xticklabels([])
        sub.tick_params(axis='x',which='both',direction='in')

        sub = fig.add_subplot(gs[3,0])
        sub.plot(x_arr,Ek,label=r'$u^2$')
        sub.plot(x_arr,Eb,label=r'$v_a^2$')
        sub.legend(loc='upper right')
        sub.axvline(x_alf,color='red',ls='--')
        sub.set_xlim(xlim)
        subt = sub.twinx()
        subt.plot(x_arr,Mach_a_scaling,ls='--',color='k',label=r'$M_a/(M_a+1)^2$')
        subt.legend(loc='lower right')
        sub.set_xlabel(r'$r/R_s$',fontsize=12)

        fig.suptitle(r'$T = {:.1f}$ hrs'.format(T_wave_hours))
        # fig.tight_layout(rect=[0,0,1,0.95])
        fig.subplots_adjust(hspace=0,bottom=0.08,top=0.95,left=0.09,right=0.94)
        # plt.show()
        # exit()
        fig.savefig('./fig_T_{:.2f}_hr.png'.format(T_wave_hours))
        plt.close(fig)

    # sub.set_xlabel(r'$x$')
    # sub.set_ylabel(r'$\sigma_c$')
    # sub.legend()
    # plt.show()