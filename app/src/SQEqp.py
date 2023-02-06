from os import environ
environ["MKL_NUM_THREADS"] = "1"
environ["NUMEXPR_NUM_THREADS"] = "1"
environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from numba import jit
from math import erf
from numba.core import types
from numba.typed import Dict
import json
import scipy
import sys


def calculate_charges(molecule,impl='orig'):
    if impl is None or impl == 'orig':
        charges = orig_sqeqp_calculate(molecule.bonds,
                              molecule.precalc_bond_hardnesses,
                              molecule.coordinates,
                              molecule.total_chg,
                              molecule.precalc_params)
    elif impl == 'new':
        charges = new_sqeqp_calculate(molecule.bonds,
                              molecule.precalc_bond_hardnesses,
                              molecule.coordinates,
                              molecule.total_chg,
                              molecule.precalc_params)

    elif impl == 'triang':
        charges = triang_sqeqp_calculate(molecule.bonds,
                              molecule.precalc_bond_hardnesses,
                              molecule.coordinates,
                              molecule.total_chg,
                              molecule.precalc_params)

    else:
        raise RuntimeError(f"{impl}_calculate_charges() not available")

    return charges[:molecule.calculated_atoms]


def mynan(name,x):
	print(name,x.shape,np.any(np.isnan(x)),file=sys.stderr)

# @jit(nopython=False, cache=True)
def new_sqeqp_calculate(bonds,
                    precalc_bond_hardnesses,
                    coordinates,
                    total_chg,
                    precalc_params):
    # this method is the same for both SQEqp and SQEqps
    electronegativities = precalc_params[:, 0]
    hardnesses = precalc_params[:, 1]
    radiuses = precalc_params[:, 2]
    initial_charges = precalc_params[:, 3]

    num_of_ats = len(coordinates)
    num_of_bonds = len(bonds)
    T = np.zeros((num_of_bonds, num_of_ats), dtype=np.float64)

#    for i,(a1, a2, _) in enumerate(bonds):
#        T[i, a1] += 1
#        T[i, a2] -= 1
    bidx = np.arange(num_of_bonds)
    T[(bidx,bonds[:,0])] = 1
    T[(bidx,bonds[:,1])] = -1
    

#    matrix[np.diag_indices(num_of_ats)] = hardnesses[i]
    baux = np.broadcast_to(radiuses,(num_of_ats,num_of_ats))
    d0 = np.sqrt(baux + baux.T)


    d2 = np.resize(coordinates,(num_of_ats,num_of_ats,3))
#    d2 = np.empty((num_of_ats,num_of_ats,3))
#    for i in range(num_of_ats):
#        d2[i,:,:] = coordinates

    d2 = d2 - np.swapaxes(d2,0,1)
    distances = np.linalg.norm(d2,axis=2)
    matrix = scipy.special.erf(distances/d0)/distances

#    for i in range(num_of_ats):
#        matrix[i, i] = hardnesses[i]
    matrix[np.diag_indices(num_of_ats)] = hardnesses

#        i_radius = radiuses[i]
#        ix, iy, iz = coordinates[i]
#        for j, (j_radius, (jx,jy,jz)) in enumerate(zip(radiuses[i+1:],
#                                                       coordinates[i+1:]),
#                                                       i + 1):
#            d0 = np.sqrt(i_radius + j_radius)
#            distance = np.sqrt((ix - jx) ** 2 + (iy - jy) ** 2 + (iz - jz) ** 2)
#            matrix[i, j] = matrix[j, i] = erf(distance / d0) / distance
#            matrix[i, j] = matrix[j, i] = erf(distance / d0[i,j]) / distance
    initial_charges -= (np.sum(initial_charges) - total_chg) / len(initial_charges)
    A_sqe = np.dot(T, np.dot(matrix, T.T))
    for i, hardness in enumerate(precalc_bond_hardnesses):
        A_sqe[i, i] += hardness
    electronegativities -= np.dot(matrix, initial_charges)
    electronegativities += hardnesses * initial_charges
    B_sqe = np.dot(T, electronegativities)
    r = np.dot(np.linalg.solve(A_sqe, B_sqe), T) + initial_charges
    return r

# @jit(nopython=False, cache=True)
def triang_sqeqp_calculate(bonds,
                    precalc_bond_hardnesses,
                    coordinates,
                    total_chg,
                    precalc_params):
    # this method is the same for both SQEqp and SQEqps
    electronegativities = precalc_params[:, 0].astype(np.float32)
    hardnesses = precalc_params[:, 1].astype(np.float32)
    radiuses = precalc_params[:, 2].astype(np.float32)
    initial_charges = precalc_params[:, 3].astype(np.float32)

    num_of_ats = len(coordinates)
    num_of_bonds = len(bonds)
    T = np.zeros((num_of_bonds, num_of_ats), dtype=np.float32)
    matrix = np.zeros((num_of_ats, num_of_ats), dtype=np.float32)

    bidx = np.arange(num_of_bonds)
    T[(bidx,bonds[:,0])] = 1.0
    T[(bidx,bonds[:,1])] = -1.0
    
    idx = np.triu_indices(num_of_ats,k=1)

    baux = np.broadcast_to(radiuses,(num_of_ats,num_of_ats))
#    d0 = np.sqrt(baux + baux.T)
    d0 = np.sqrt((baux + baux.T)[idx])

    #mynan('d0',d0)

    d2 = np.resize(coordinates,(num_of_ats,num_of_ats,3))
#    d2 = d2 - np.swapaxes(d2,0,1)
    d2 = (d2 - np.swapaxes(d2,0,1))[idx[0],idx[1],:]

#    distances = np.linalg.norm(d2,axis=2)
    distances = np.linalg.norm(d2,axis=1)
    #mynan('distances',distances)

#    matrix = scipy.special.erf(distances/d0)/distances
   
    uf = scipy.special.erf(distances/d0)/distances
    #print('uf:',len(uf),list(uf),file=sys.stderr)
    #print('atoms',num_of_ats,file=sys.stderr)
    #mynan('uf',uf)

    #mynan('matrix0',matrix)
    #print('idx',idx[0].shape,file=sys.stderr)
    np.put(matrix,idx[0]*num_of_ats+idx[1],uf)
    matrix += matrix.T
    #mynan('matrix',matrix)

    matrix[np.diag_indices(num_of_ats)] = hardnesses
    #mynan('matrix2',matrix)

    initial_charges -= (np.sum(initial_charges) - total_chg) / len(initial_charges)
    A_sqe = np.dot(T, np.dot(matrix, T.T))
    for i, hardness in enumerate(precalc_bond_hardnesses):
        A_sqe[i, i] += hardness
    electronegativities -= np.dot(matrix, initial_charges)
    electronegativities += hardnesses * initial_charges
    B_sqe = np.dot(T, electronegativities)
    r = np.dot(np.linalg.solve(A_sqe, B_sqe), T) + initial_charges
    return r

@jit(nopython=True, cache=True)
def orig_sqeqp_calculate(bonds,
                    precalc_bond_hardnesses,
                    coordinates,
                    total_chg,
                    precalc_params):
    # this method is the same for both SQEqp and SQEqps
    electronegativities = precalc_params[:, 0]
    hardnesses = precalc_params[:, 1]
    radiuses = precalc_params[:, 2]
    initial_charges = precalc_params[:, 3]
    num_of_ats = len(coordinates)
    num_of_bonds = len(bonds)
    T = np.zeros((num_of_bonds, num_of_ats), dtype=np.float64)
    matrix = np.empty((num_of_ats, num_of_ats), dtype=np.float64)
    for i,(a1, a2, _) in enumerate(bonds):
        T[i, a1] += 1
        T[i, a2] -= 1
    for i in range(num_of_ats):
        matrix[i, i] = hardnesses[i]
        i_radius = radiuses[i]
        ix, iy, iz = coordinates[i]
        for j, (j_radius, (jx,jy,jz)) in enumerate(zip(radiuses[i+1:],
                                                       coordinates[i+1:]),
                                                       i + 1):
            d0 = np.sqrt(i_radius + j_radius)
            distance = np.sqrt((ix - jx) ** 2 + (iy - jy) ** 2 + (iz - jz) ** 2)
            matrix[i, j] = matrix[j, i] = erf(distance / d0) / distance
    initial_charges -= (np.sum(initial_charges) - total_chg) / len(initial_charges)
    A_sqe = np.dot(T, np.dot(matrix, T.T))
    for i, hardness in enumerate(precalc_bond_hardnesses):
        A_sqe[i, i] += hardness
    electronegativities -= np.dot(matrix, initial_charges)
    electronegativities += hardnesses * initial_charges
    B_sqe = np.dot(T, electronegativities)
    r = np.dot(np.linalg.solve(A_sqe, B_sqe), T) + initial_charges
    return r


@jit(nopython=True, cache=True)
def precalculate_parameters_SQEqps(atomic_types, bonds_types, surfaces, parameters, bond_hardnesses):
    n_atoms = len(atomic_types)
    n_bonds = len(bonds_types)
    precalc_params = np.empty((n_atoms, 4), dtype=np.float32)
    precalc_bond_hardnesses = np.empty(n_bonds, dtype=np.float32)
    for i in range(n_atoms):
        symbol_i = atomic_types[i]
        surface = surfaces[i]
        electronegativity, hardness, width, q0, q0_cor, hardness_cor, electronegativity_cor, width_cor = parameters[symbol_i]
        precalc_params[i] = (-electronegativity + electronegativity_cor * surface,
                         hardness + hardness_cor * surface,
                         2 * (width + width_cor * surfaces[i]) ** 2,
                         q0 + q0_cor * surface)
    for i in range(n_bonds):
        precalc_bond_hardnesses[i] = bond_hardnesses[bonds_types[i]]
    return precalc_params, precalc_bond_hardnesses


@jit(nopython=True, cache=True)
def precalculate_parameters_SQEqp(atomic_types, bonds_types, parameters, bond_hardnesses):
    n_atoms = len(atomic_types)
    n_bonds = len(bonds_types)
    precalc_params = np.empty((n_atoms, 4), dtype=np.float64)
    precalc_bond_hardnesses = np.empty(n_bonds, dtype=np.float64)
    for i in range(n_atoms):
        symbol_i = atomic_types[i]
        electronegativity, hardness, width, q0, = parameters[symbol_i]
        precalc_params[i] = np.array((-electronegativity,
                             hardness,
                             2 * width ** 2,
                             q0), dtype=np.float64)
    for i in range(n_bonds):
        if bonds_types[i] not in bond_hardnesses:
            print(bonds_types[i])
        precalc_bond_hardnesses[i] = bond_hardnesses[bonds_types[i]]
    return precalc_params, precalc_bond_hardnesses

def load_parameters(root_dir):
    params_SQEqps = json.load(open(f"{root_dir}/parameters/parameters_SQEqps.json"))
    parameters_SQEqps = Dict.empty(key_type=types.unicode_type,
                                   value_type=types.float32[:])
    for key, value in params_SQEqps["atom"]["data"].items():
        parameters_SQEqps[key] = np.array(value, dtype=np.float32)
    bond_hardnesses_SQEqps = Dict.empty(key_type=types.unicode_type,
                                        value_type=types.float64)
    for key, value in params_SQEqps["bond"]["data"].items():
        bond_hardnesses_SQEqps[key] = value
    params_SQEqp = json.load(open(f"{root_dir}/parameters/parameters_SQEqp.json"))
    parameters_SQEqp = Dict.empty(key_type=types.unicode_type,
                                   value_type=types.float32[:])
    for key, value in params_SQEqp["atom"]["data"].items():
        parameters_SQEqp[key] = np.array(value, dtype=np.float32)
    bond_hardnesses_SQEqp = Dict.empty(key_type=types.unicode_type,
                                        value_type=types.float64)
    for key, value in params_SQEqp["bond"]["data"].items():
        bond_hardnesses_SQEqp[key] = value
    return parameters_SQEqp, bond_hardnesses_SQEqp, parameters_SQEqps, bond_hardnesses_SQEqps
