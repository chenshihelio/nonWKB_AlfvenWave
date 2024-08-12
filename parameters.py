# wave periods to be solved
T_wave_arr_hours = [5,10,30,50,100]

# background fields---------
Uinf = 600 # km/s
Va0 = 500 # km/s
rho0 = 1e10 # We only need U and Va. Density is not used.

alpha = 2.0
beta = 8.5

# domain--------------
x1 = 1
x2 = 215
dx = 0.1

# constants------------
Rsun = 6.957E5 # km, solar radius

if __name__ == "__main__":
    print('You are running parameters.py. It is an input file, nothing will happen!')