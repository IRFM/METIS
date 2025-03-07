# -*- coding: iso-8859-1 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
# If you need to restart vim tabs then :retab
'''
    In this file you have a function that gets some
    equilibrium IDS fields and transform them in
    matlab objects stored in a .mat file
'''
# Standard python modules
from __future__ import (unicode_literals, absolute_import,  \
                        print_function, division)
import argparse
import numpy as np
import scipy
import scipy.interpolate as interp
import scipy.io
import os

# Local python modules
import imas

def equilibrium_python_to_matlab(shot, run, occurrence, user, machine):
    ''' Function reading equilibrium IMAS data and writing it to
        matlab format, it creates a .mat file
    '''
    grid_points_interp = 60 # Number of points in R and Z for cartesian grid interpolation
    print('shot       =', shot)
    print('run        =', run)
    print('occurrence =', occurrence)
    print('user       =', user)
    print('machine    =', machine)
    print('Reading data...')

    # Open shot and specific run of machine
    idd = imas.ids(shot, run)
    idd.open_env(user, machine, '3')

    # Get IDS
    if (occurrence == 0):
        idd.equilibrium.get()
        idd.wall.get()
    else:
        idd.equilibrium.get(occurrence)
        idd.wall.get(occurrence)

    print(' ')
    print('Time in equilibrium IDS =', idd.equilibrium.time)
    timeEquiIDS = idd.equilibrium.time

    # Array with all times requested
    lenArrTimes = len(idd.equilibrium.time)

    equi_tSlice = idd.equilibrium.time_slice[0]
    equi_space  = idd.equilibrium.time_slice[0].ggd[0]

    # For EQUINOX run
    NbrPoints = len(equi_space.grid.space[0].objects_per_dimension[0].object)
    print('NbrPoints (number of grid points) =', NbrPoints)

    # Declaration of arrays
    ip       = np.zeros(lenArrTimes)
    q_95     = np.zeros(lenArrTimes)
    q_axis   = np.zeros(lenArrTimes)
    li_3     = np.zeros(lenArrTimes)
    beta_pol = np.zeros(lenArrTimes)
    w_mhd    = np.zeros(lenArrTimes)
    mag_ax_R = np.zeros(lenArrTimes)
    mag_ax_Z = np.zeros(lenArrTimes)
    RNodes   = np.zeros(NbrPoints)
    ZNodes   = np.zeros(NbrPoints)

    Psi_val    = np.zeros((lenArrTimes, NbrPoints))
    Psi_interp = np.zeros((lenArrTimes, grid_points_interp, grid_points_interp))

    boundPlasma = np.zeros((lenArrTimes, 2, 201))

    xPoint   = np.zeros((lenArrTimes, 2))
    limPoint = np.zeros((lenArrTimes, 2))

    wall = np.zeros((2, \
           len(idd.wall.description_2d[0].limiter.unit[0].outline.r)))

    b0 = np.zeros(lenArrTimes)

    lenProf1d = len(equi_tSlice.profiles_1d.rho_tor)
    print('length of rho_tor = ', lenProf1d)
    rho_tor     = np.zeros((lenArrTimes, lenProf1d))
    q_safFac    = np.zeros((lenArrTimes, lenProf1d))
    elong       = np.zeros((lenArrTimes, lenProf1d))
    triang_up   = np.zeros((lenArrTimes, lenProf1d))
    triang_low  = np.zeros((lenArrTimes, lenProf1d))
    j_tor       = np.zeros((lenArrTimes, lenProf1d))
    pressure    = np.zeros((lenArrTimes, lenProf1d))
    f_df_dpsi   = np.zeros((lenArrTimes, lenProf1d))
    dpress_dpsi = np.zeros((lenArrTimes, lenProf1d))
    psi_prof    = np.zeros((lenArrTimes, lenProf1d))

    # Put data into declared arrays
    # Wall
    wall[0, :] = idd.wall.description_2d[0].limiter.unit[0].outline.r
    wall[1, :] = idd.wall.description_2d[0].limiter.unit[0].outline.z

    # b0 vacuum toroidal field and r0
    b0 = idd.equilibrium.vacuum_toroidal_field.b0
    r0 = idd.equilibrium.vacuum_toroidal_field.r0

    # Read [R,Z] coordinates of nodes
    for i in range(NbrPoints):
        RNodes[i]  = equi_space.grid.space[0].objects_per_dimension[0]. \
                     object[i].geometry[0]
        ZNodes[i]  = equi_space.grid.space[0].objects_per_dimension[0]. \
                     object[i].geometry[1]

    # Min and max values used later for cartesian interpolation
    min_RNodes  = np.min(RNodes)
    max_RNodes  = np.max(RNodes)
    min_ZNodes  = np.min(ZNodes)
    max_ZNodes  = np.max(ZNodes)

    # Create cartesian grid for interpolation of the data
    RInterp, ZInterp = np.mgrid[ \
                       min_RNodes:max_RNodes:complex(0,grid_points_interp), \
                       min_ZNodes:max_ZNodes:complex(0,grid_points_interp)]

    print('RInterp.shape =', RInterp.shape)
    print('ZInterp.shape =', ZInterp.shape)

    ptsNodes = np.vstack((RNodes, ZNodes)).T
    print('ptsNodes.shape =', ptsNodes.shape)
    print(' ')

    # Loop over time
    for timeit in range(lenArrTimes):

        print('-----')
        print('It:', timeit, ', time in equilibrium IDS =', \
                                idd.equilibrium.time[timeit])

        equi_tSlice = idd.equilibrium.time_slice[timeit]
        equi_space  = idd.equilibrium.time_slice[timeit].ggd[0]

        ip[timeit]       = equi_tSlice.global_quantities.ip
        q_95[timeit]     = equi_tSlice.global_quantities.q_95
        q_axis[timeit]   = equi_tSlice.global_quantities.q_axis
        li_3[timeit]     = equi_tSlice.global_quantities.li_3
        beta_pol[timeit] = equi_tSlice.global_quantities.beta_pol
        w_mhd[timeit]    = equi_tSlice.global_quantities.w_mhd
        mag_ax_R[timeit] = equi_tSlice.global_quantities.magnetic_axis.r
        mag_ax_Z[timeit] = equi_tSlice.global_quantities.magnetic_axis.z

        # Psi and plasma boundary
        Psi_val[timeit, :] = equi_space.psi[0].values

        print('len equi_tSlice.boundary.outline.r =', \
                   len(equi_tSlice.boundary.outline.r))
        boundPlasma[timeit, 0, :] = np.interp(np.linspace(0, 1, 201), \
                    np.linspace(0, 1, len(equi_tSlice.boundary.outline.r)), \
                    equi_tSlice.boundary.outline.r)
        boundPlasma[timeit, 1, :] = np.interp(np.linspace(0, 1, 201), \
                    np.linspace(0, 1, len(equi_tSlice.boundary.outline.z)), \
                    equi_tSlice.boundary.outline.z)

        if (equi_tSlice.boundary.x_point[0].r != 0):
            xPoint[timeit, 0] = equi_tSlice.boundary.x_point[0].r
            xPoint[timeit, 1] = equi_tSlice.boundary.x_point[0].z
        else:
            xPoint[timeit, 0] = None
            xPoint[timeit, 1] = None

        if (equi_tSlice.boundary.active_limiter_point.r != 0):
            limPoint[timeit, 0] = equi_tSlice.boundary.active_limiter_point.r
            limPoint[timeit, 1] = equi_tSlice.boundary.active_limiter_point.z
        else:
            limPoint[timeit, 0] = None
            limPoint[timeit, 1] = None

        # Interpolation
        Psi_interp[timeit, :, :] = interp.griddata(ptsNodes, \
                                                   Psi_val[timeit, :], \
                                                   (RInterp, ZInterp), \
                                                   method='linear')

        # Compute profiles 1d
        rho_tor[timeit, :]     = equi_tSlice.profiles_1d.rho_tor
        q_safFac[timeit, :]    = equi_tSlice.profiles_1d.q
        elong[timeit, :]       = equi_tSlice.profiles_1d.elongation
        triang_up[timeit, :]   = equi_tSlice.profiles_1d.triangularity_upper
        triang_low[timeit, :]  = equi_tSlice.profiles_1d.triangularity_lower
        j_tor[timeit, :]       = equi_tSlice.profiles_1d.j_tor
        pressure[timeit, :]    = equi_tSlice.profiles_1d.pressure
        f_df_dpsi[timeit, :]   = equi_tSlice.profiles_1d.f_df_dpsi
        dpress_dpsi[timeit, :] = equi_tSlice.profiles_1d.dpressure_dpsi
        psi_prof[timeit, :]    = equi_tSlice.profiles_1d.psi

    # Write .mat file
    filepath = os.path.dirname(os.path.realpath(__file__))
    print(' ')
    print('filepath =', filepath)
    print(' ')
    scipy.io.savemat('data_vactheqx_Shot' + \
        str(shot) + '_Run' +  str(run) + '_Occ' + str(occurrence) + '_' + \
        user + '_' + machine + '.mat', \
        {'timeEquiIDS': timeEquiIDS, \
         'wall': wall, \
         'ip': ip, \
         'q_95': q_95, \
         'q_axis': q_axis, \
         'li_3': li_3, \
         'beta_pol': beta_pol, \
         'w_mhd': w_mhd, \
         'mag_ax_R': mag_ax_R, \
         'mag_ax_Z': mag_ax_Z, \
         'RNodes': RNodes, \
         'ZNodes': ZNodes, \
         'Psi_val': Psi_val, \
         'boundPlasma': boundPlasma, \
         'xPoint': xPoint, \
         'RInterp': RInterp, \
         'ZInterp': ZInterp, \
         'Psi_interp': Psi_interp, \
         'rho_tor': rho_tor, \
         'q_safFac': q_safFac, \
         'elong': elong, \
         'triang_up': triang_up, \
         'triang_low': triang_low, \
         'j_tor': j_tor, \
         'pressure': pressure, \
         'f_df_dpsi': f_df_dpsi, \
         'dpress_dpsi': dpress_dpsi, \
         'psi_prof': psi_prof, \
         'limPoint': limPoint, \
        })

if (__name__ == '__main__'):

    # Parse input arguments
    parser = argparse.ArgumentParser(description= \
            'Function reading equilibrium IDS and writing matlab .mat file')
    parser.add_argument('shot', type=int, nargs='?', default=52187, \
                        help='shot, default=52187')
    parser.add_argument('run', type=int, nargs='?', default=0, \
                        help='run, default=0')
    parser.add_argument('occurrence', type=int, nargs='?', default=0, \
                        help='occurrence, default=0')
    parser.add_argument('user', type=str, nargs='?', default='imas_public', \
                        help='user, default=imas_public')
    parser.add_argument('machine', type=str, nargs='?', default='west', \
                        help='machine, default=west')
    #parser.add_argument('--fast', action='store_true', \
    #                    help='fast calculation')

    args = parser.parse_args()

    equilibrium_python_to_matlab(args.shot, args.run, args.occurrence, \
                                 args.user, args.machine)

