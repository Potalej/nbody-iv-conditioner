"""
This is an API to facilitate the comunication with the Fortran subroutines
and functions.
"""
import numpy as np
from ._core import *

def fix_vector (vector:np.array):
    if not vector.flags.f_contiguous:
        return np.asfortranarray(vector)
    else:
        return vector

class Utils:
    def center_of_mass (ms:np.array, qs:np.array):
        """
        Evaluates the center of mass of a system of N bodies.
        """
        qs = fix_vector(qs)
        return utils_mod.center_of_mass(ms, qs)
    
    def total_linear_momentum (ps:np.array):
        """
        Total linear momentum
        """
        ps = fix_vector(ps)
        return utils_mod.total_linear_momentum(ps)
    
    def kinect_energy (ms:np.array, ps:np.array):
        """
        Kinect energy
        """
        ps = fix_vector(ps)
        return utils_mod.kinect_energy(ms, ps)
    
    def potential_energy (ms:np.array, qs:np.array, G:float=1.0, soft:float=0.0):
        """
        Potential energy
        """
        qs = fix_vector(qs)
        return utils_mod.potential_energy(ms, qs, G, soft)
    
    def total_energy (ms:np.array, qs:np.array, ps:np.array, G:float=1.0, soft:float=0.0):
        """
        Total energy
        """
        qs = fix_vector(qs)
        ps = fix_vector(ps)
        return utils_mod.total_energy(ms, qs, ps, G, soft)

    def total_angular_momentum (qs:np.array, ps:np.array):
        """
        Total angular momentum
        """
        qs = fix_vector(qs)
        ps = fix_vector(ps)
        return utils_mod.total_angular_momentum(qs, ps)
    
    def individual_inertia_tensor (m:float, qs:np.array):
        """
        Individual inertia tensor
        """
        return utils_mod.individual_inertia_tensor(m, qs)
    
    def general_inertia_tensor (ms:np.array, qs:np.array):
        """
        General inertia tensor
        """
        qs = fix_vector(qs)
        return utils_mod.total_angular_momentum(ms, qs)


 
class Conditioners:
    def comori (ms:np.array, qs:np.array):
        """
        Move the center of mass to the origin
        """
        qs = fix_vector(qs)
        conditioners_mod.comori(ms, qs)
        return ms, qs
    
    def contotangmom (J:np.array, ms:np.array, qs:np.array, ps:np.array):
        """
        Conditionates the total angular momentum
        """
        qs = fix_vector(qs)
        ps = fix_vector(ps)
        ier = 0
        conditioners_mod.contotangmom(J, ms, qs, ps, ier)
        return ms, qs, ps
    
    def contotene (E:float, ms:np.array, qs:np.array, ps:np.array, G:float=1.0, soft:float=0.0):
        """
        Conditionates the total energy
        """
        qs = fix_vector(qs)
        ps = fix_vector(ps)
        conditioners_mod.contotene(E, ms, qs, ps, G, soft)
        return ms, qs, ps
    
    def conpotenesof (E:float, ms:np.array, qs:np.array, ps:np.array, initial_guess:float, G:float=1.0, soft:float=0.0):
        """
        Conditionates the potential energy softened
        """
        qs = fix_vector(qs)
        ps = fix_vector(ps)
        conditioners_mod.conpotenesof(E, ms, qs, ps, G, soft, initial_guess)
        return ms, qs, ps
    
    def conpotenesofequ (E:float, ms:np.array, qs:np.array, initial_guess:float, G:float=1.0, soft:float=0.0):
        """
        Conditionates the potential energy softened and in equilibrium
        """
        qs = fix_vector(qs)
        conditioners_mod.conpotenesofequ(E, ms, qs, G, soft, initial_guess)
        return ms, qs
    
    def contotlinmom (P:np.array, ps:np.array, ws:np.array):
        """
        Conditionates the total linear momentum
        """
        ps = fix_vector(ps)
        conditioners_mod.contotlinmom(P, ps, ws)
        return ps
    
    ##########
    # methods
    def coninivalite (ms:np.array, qs:np.array, ps:np.array, G:float=1.0, soft:float=0.0, E:float=0.0, J:np.array=np.zeros(3), P:np.array=np.zeros(3), num_iter_max:int=10):
        """
        Conditionates the first integrals iteratively
        """
        qs = fix_vector(qs)
        ps = fix_vector(ps)
        ier = 0
        conditioners_mod.coninivalite(ms, qs, ps, G, soft, ier, E, J, P, num_iter_max)
        return ms, qs, ps
    
    def coninivaldir (ms:np.array, qs:np.array, ps:np.array, G:float=1.0, E:float=0.0, J:np.array=np.zeros(3), P:np.array=np.zeros(3)):
        """
        Conditionates the first integrals directly
        """
        qs = fix_vector(qs)
        ps = fix_vector(ps)
        ier = 0
        conditioners_mod.coninivaldir(ms, qs, ps, G, E, J, P)
        return ms, qs, ps
    
    def coninivalaar (ms:np.array, qs:np.array, ps:np.array, G:float=1.0):
        """
        Conditionates the system to Henon units using the Aarseth method
        """
        qs = fix_vector(qs)
        ps = fix_vector(ps)
        ier = 0
        conditioners_mod.coninivalaar(ms, qs, ps, G, ier)
        return ms, qs, ps
    
    def coninivalaarmod (ms:np.array, qs:np.array, ps:np.array, G:float=1.0, soft:float=0.0):
        """
        Conditionates the system to Henon units using the modified Aarseth method
        """
        qs = fix_vector(qs)
        ps = fix_vector(ps)
        ier = 0
        conditioners_mod.coninivalaarmod(ms, qs, ps, G, soft, ier)
        return ms, qs, ps