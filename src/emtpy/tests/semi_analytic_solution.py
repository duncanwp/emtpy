__author__ = 'duncan'
import numpy as np
from math import sqrt
from emtpy.constants import *


def infinite_square_box_energy(n, L):
    """
        Return the energies od a particle in an infinite box
    @param n:
    @param L:
    @return:
    """
    from emtpy.constants import hbar, me, eV
    from math import pi
    return ((hbar**2 * pi**2 * n**2)/(2.0*me*L**2))/eV


def initial_guess_at_roots(func, delta, a, b):
    """
        Static method designed to have a first go at finding the zeros of a function
    @param func:
    @param delta:
    @param a:
    @param b:
    @return:
    """
    n = int((b-a)/delta)
    increments = np.linspace(a, b, n)
    zeros = []
    prev = a
    ar = np.array([func(inc) for inc in increments])
    zero_crossings = np.where(np.diff(np.sign(ar)))[0] + 1
    return zero_crossings*delta
    # return [ (zc*delta, (zc-1)*delta) for zc in zero_crossings ]


def finite_square_box_energies(L, V_0, tol=1.0E-3):
    """
        Helper method
    @param L:
    @param V_0:
    @param tol:
    @return:
    """
    AS = AnalyticSolution(L, V_0, tol)
    return AS.finite_square_box_energies()


class AnalyticSolution(object):

    def __init__(self, L, V_0, tol=1.0E-3):
        self.L = L
        self.V_0 = V_0
        self.tol = tol

        self.u_0_sq = (me * L**2 * V_0)/(2 * hbar**2)

        #u_0_sq = 20

    def a(self, v):
        return np.sqrt(self.u_0_sq - v**2)

    def b_symm(self, v):
        return v * np.tan(v)

    def b_anti_symm(self, v):
        return - v / np.tan(v)

    def symmetric_solution_vec(self, v):
        a = np.sqrt(self.u_0_sq - v**2)
        where_are_nans = np.isnan(a)
        a[where_are_nans] = 0.0
        b = (v * np.tan(v))
        val = a - b
        return val

    def symmetric_solution(self, v):
        if v**2 < self.u_0_sq:
            a = sqrt(self.u_0_sq - v**2)
        else:
            a = sqrt(v**2 - self.u_0_sq)
        b = (v * np.tan(v))
        val = a - b
        return val

    def anti_symmetric_solution_vec(self, v):
        from numpy import isnan
        a = np.sqrt(self.u_0_sq - v**2)
        where_are_nans = isnan(a)
        a[where_are_nans] = 0.0
        return a + (v / np.tan(v))


    def anti_symmetric_solution(self, v):
        if v**2 < self.u_0_sq:
            a = sqrt(self.u_0_sq - v**2)
        else:
            a = sqrt(v**2 - self.u_0_sq)
        b = v / np.tan(v)
        return a + b

    def energy(self, v):
        return ((2.0 * hbar**2 * v**2)/(me*self.L**2))/eV

    def finite_square_box_energies(self):
        from emtpy.constants import me, eV, hbarev, hbar
        from numpy import sqrt, tan, arange, linspace, fromfunction, array
        from scipy.optimize import fsolve, bisect, brentq, newton
        import numpy as np
        import math


        # initial_guess = sqrt(linspace(0.0 +  sqrt(u_0_sq)/5.0, sqrt(u_0_sq),5))
        # initial_guess = linspace(0.0 + (u_0_sq)/5.0, (u_0_sq),5)
        # symm_sol = brentq(symmetric_solution, 0.0 + (u_0_sq)/50.0, (u_0_sq))
        # sym_en = energy(symm_sol)
        #anti_symm_sol = energy(bisect(anti_symmetric_solution, 0.0, u_0_sq))
        # grid = np.linspace(0.0, math.sqrt(self.u_0_sq), 100)
        symm_guess = initial_guess_at_roots(self.symmetric_solution, sqrt(self.u_0_sq)/100.0, 0.0, sqrt(self.u_0_sq))

        # symm_sol = brentq(self.symmetric_solution, initial_guess[0][0], initial_guess[0][1])

        solutions, info, ierr, mesg = fsolve(self.symmetric_solution_vec, symm_guess, full_output=True)
        if ierr != 1:
            print mesg
            assert False

        anti_symm_guess = initial_guess_at_roots(self.anti_symmetric_solution, sqrt(self.u_0_sq)/100.0, 0.0, sqrt(self.u_0_sq))

        anti_solutions, info, ierr, mesg = fsolve(self.anti_symmetric_solution_vec, anti_symm_guess, full_output=True)
        if ierr != 1:
            print mesg
            assert False

        # solutions, info, ierr, mesg = fsolve(symmetric_solution, (linspace(0.0 +  (u_0_sq**2)/5.0, (u_0_sq**2),5)), full_output=True)
        # anti_solutions, info, ierr, mesg = fsolve(anti_symmetric_solution, (linspace(0.0 +  (u_0_sq**2)/5.0, (u_0_sq**2),5)), full_output=True)

        sols = set(list(abs(solutions.round(5)))+list(abs(anti_solutions.round(5))))
        energies = [self.energy(e) for e in sols if e < sqrt(self.u_0_sq)]
        energies.sort()
        return energies