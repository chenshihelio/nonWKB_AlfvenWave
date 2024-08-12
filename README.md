# nonWKB_AlfvenWave solver
Solve the two Elsasser variables given background solar wind profiles. 1D spherical coordinate system is used. Please refer to **nonWKB_Alfven_wave_solver.pdf** for the equation set that is solved by this solver.

**parameters.py**: Input parameters.

**background_fields.py**: Functions that calculate background fields, including solar wind speed $U(r)$, Alfven speed $V_A(r)$, and density $\rho(r)$, and their 1st-order and 2nd-order derivatives. Note that density is not directly used in integrating the equation but is used in calculating $V_A(r)$ and $U(r)$. The 2nd-order derivatives are used only at the Aflven point.
By default, the models are taken from *Waves from the Sun?* by Velli et al. (1991). The users can define their own background fields either analytically or from interpolated satellite data (hint: scipy.interpolate package can return a callable function from interpolated data). 
Note: It is assumed that $U(r)$ and $V_A(r)$ only have one crossing, i.e. there is only one Alfven point. 

**functions_for_integral.py**: Core functions that are used in integrating the equation set. *The users do not need to modify this file.*

**main.py**: Main file to be run. By default, a set of output files will be generated in ***.npy*** format.

**plot_data.py**: Generate plots from the ***.npy*** data.


---
Copyright (C) 2024  Chen Shi (cshi1993@ucla.edu)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
