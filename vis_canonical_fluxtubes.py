<<<<<<< HEAD
#! /home/jensv/anaconda/bin/python
=======
#! /Users/vonderlinden2/anaconda/bin/python
>>>>>>> canonical_flux_tubes/master
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 14:38:10 2016

@author: Jens von der Linden
"""

from datetime import date
from datetime import datetime
import numpy as np
import os
from scipy.constants import elementary_charge, proton_mass
from glob import glob
import sys
<<<<<<< HEAD
visit_path1 = "/home/jensv/visit/visit2_10_3.linux-x86_64/2.10.3/linux-x86_64/lib/site-packages"
visit_path2 = "/home/jensv/visit/visit2_10_3.linux-x86_64/bin/"
sys.path.append("/home/jensv/visit/visit2_10_3.linux-x86_64/2.10.3/linux-x86_64/lib/site-packages")
sys.path.append("/home/jensv/visit/visit2_10_3.linux-x86_64/bin/")
os.environ["PATH"] += os.pathsep + visit_path1
os.environ["PATH"] += os.pathsep + visit_path2
=======
>>>>>>> canonical_flux_tubes/master
import visit
import argparse

tan = (209, 178, 111, 255)
olive = (110, 117, 14, 255)
<<<<<<< HEAD
dim_grey =(105, 105, 105, 255)
black = (0, 0, 0, 255)
dark_grey = (169, 169, 169, 255)
red = (255, 0, 0, 255)
green = (0, 154, 0, 255)
navy = (0, 0, 128, 255)
aqua = (0, 255, 255, 255)
=======
red = (255, 0, 0, 255)
green = (0, 154, 0, 255)
>>>>>>> canonical_flux_tubes/master

def define_expressions(visit, alpha=8.1e5):
    r"""
    Define Visit expressions with repective alpha.
    """
    visit.DefineVectorExpression("B", "{B_x, B_y, B_z}")
    visit.DefineVectorExpression("B_norm", "{B_norm_x, B_norm_y, "
                                 "B_norm_z}")
    visit.DefineVectorExpression("B_perp", "{B_x, B_y, 0}")
    visit.DefineVectorExpression("B_para", "{0, 0, B_z}")
    visit.DefineScalarExpression("B_para_scalar", "B_z")

<<<<<<< HEAD
    visit.DefineVectorExpression("A_vacuum", "{A_vacuum_x, A_vacuum_y, 0}")
    visit.DefineVectorExpression("A", "{A_x, A_y, A_z}")

=======
>>>>>>> canonical_flux_tubes/master
    visit.DefineVectorExpression("J_smooth", "{j_x, j_y, j_z}")
    visit.DefineScalarExpression("J_smooth_mag", "sqrt(j_x^2 + j_y^2 + j_z^2)")
    visit.DefineVectorExpression("J_smooth_perp", "{j_x, j_y, 0}")
    visit.DefineVectorExpression("J_smooth_para", "{0, 0, j_z}")
<<<<<<< HEAD
=======
    visit.DefineScalarExpression("J_smooth_para_scalar", "j_z")
>>>>>>> canonical_flux_tubes/master

    visit.DefineVectorExpression("J_raw", "{j_raw_x, j_raw_y, j_raw_z}")
    visit.DefineScalarExpression("J_raw_mag", "sqrt(j_raw_x^2 +" +
                                                   "j_raw_y^2 +" +
                                                   "j_raw_z^2)")
    visit.DefineVectorExpression("J_raw_perp", "{j_raw_x, j_raw_y, 0}")
<<<<<<< HEAD
    visit.DefineVectorExpression("J_raw_para", "{0, 0, j_raw_z}")

    visit.DefineScalarExpression("divergence_B", "divergence(B)")
    visit.DefineScalarExpression("divergence_Omega_i_raw_plus",
                                 "divergence(Omega_i_raw_plus)")
    visit.DefineScalarExpression("divergence_Omega_i_plus",
                                 "divergence(Omega_i_plus)")
=======


    normalize_scalar(visit, "Te_raw", "Te_raw_normalized")
    normalize_scalar(visit, "Te", "Te_smooth_normalized")
>>>>>>> canonical_flux_tubes/master

    visit.DefineVectorExpression("J_raw_filtered_by_Te",
                                  "J_raw * Te_raw_normalized")
    visit.DefineVectorExpression("J_raw_filtered_by_Te^2",
                                  "J_raw * Te_raw_normalized^2")
    visit.DefineVectorExpression("J_raw_filtered_by_Te^3",
                                  "J_raw * Te_raw_normalized^3")
    visit.DefineVectorExpression("J_raw_filtered_by_Te^4",
                                  "J_raw * Te_raw_normalized^4")
    visit.DefineVectorExpression("J_raw_filtered_by_Te_smooth",
                                  "J_raw * Te_smooth_normalized")
    visit.DefineVectorExpression("J_raw_filtered_by_Te_smooth^2",
                                  "J_raw * Te_smooth_normalized^2")
    visit.DefineVectorExpression("J_raw_filtered_by_Te_smooth^3",
                                  "J_raw * Te_smooth_normalized^3")
    visit.DefineVectorExpression("J_raw_filtered_by_Te_smooth^4",
                                  "J_raw * Te_smooth_normalized^4")

<<<<<<< HEAD
    visit.DefineVectorExpression("u_i_plus", "{u_i_x_plus, u_i_y, u_i_z}")
    visit.DefineVectorExpression("u_i_plus_perp", "{dot(u_i_plus, {1, 0, 0}), dot(u_i_plus, "
                                 "{0, 1, 0}), 0}")
    visit.DefineVectorExpression("u_i_plus_para", "{0, 0, dot(u_i_plus, {0, 0, 1})}")

    visit.DefineVectorExpression("omega_i_plus", "{w_i_x_plus, w_i_y_plus, w_i_z_plus}")
    visit.DefineVectorExpression("omega_i_plus_perp", "{dot(omega_i_plus, {1, 0, 0}),"
                                                 "dot(omega_i_plus, {0, 1, 0}), 0})")
    visit.DefineVectorExpression("omega_i_plus_para", "{0, 0, dot(omega_i_plus, "
                                 "{0, 0, 1})}")

    visit.DefineVectorExpression("omega_i_raw_plus",
                                 "{w_i_raw_x_plus, w_i_raw_y_plus, w_i_raw_z_plus}")
    visit.DefineVectorExpression("omega_i_raw_plus_perp", "{dot(omega_i_raw_plus, {1, 0, 0}),"
                                 "dot(omega_i_raw_plus, {0, 1, 0}), 0})")
    visit.DefineVectorExpression("omega_i_raw_plus_para", "{0, 0, dot(omega_i_raw_plus, "
                                 "{0, 0, 1})}")

    visit.DefineVectorExpression("Omega_i_plus",
                                 str(elementary_charge) + "*B +" + str(proton_mass) +
                                 "*omega_i_plus")
    visit.DefineVectorExpression("Omega_i_plus_perp", "{dot(Omega_i_plus, {1, 0, 0}), "
                                 "dot(Omega_i_plus, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_plus_para", "{0, 0, dot(Omega_i_plus, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_plus_para_scalar",
                                 "dot(Omega_i_plus, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_raw_plus",
                                 str(elementary_charge) + "*B +" + str(proton_mass) +
                                 "*omega_i_raw_plus")
    visit.DefineVectorExpression("Omega_i_raw_plus_perp", "{dot(Omega_i_raw_plus, {1, 0, 0}), "
                                 "dot(Omega_i_raw_plus, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_raw_plus_para", "{0, 0, dot(Omega_i_raw_plus, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_raw_plus_para_scalar",
                                 "dot(Omega_i_raw_plus, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_plus_density_dependence",
                                 "n*(%e *B + %e *omega_i_plus) +"
                                 "cross(gradient(n), %e *A + %e * u_i_plus)"
                                 % (elementary_charge, proton_mass,
                                    elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_plus_density_dependence_perp",
                                 "{dot(Omega_i_plus_density_dependence, {1, 0, 0}), "
                                 "dot(Omega_i_plus_density_dependence, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_plus_density_dependence_para",
                                 "{0, 0, dot(Omega_i_plus_density_dependence, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_plus_density_dependence_para_scalar",
                                 "dot(Omega_i_plus_density_dependence, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_raw_plus_density_dependence",
                                 "n*(%e *B + %e *omega_i_raw_plus) +"
                                 " cross(gradient(n), %e *A + %e * u_i_plus)"
                                 % (elementary_charge, proton_mass,
                                    elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_raw_plus_density_dependence_perp",
                                 "{dot(Omega_i_raw_plus_density_dependence, {1, 0, 0}), "
                                 "dot(Omega_i_raw_plus_density_dependence, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_raw_plus_density_dependence_para",
                                 "{0, 0, dot(Omega_i_raw_plus_density_dependence, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_raw_plus_density_dependence_para_scalar",
                                 "dot(Omega_i_raw_plus_density_dependence, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_plus_times_density",
                                 "n*(%e *B + %e *omega_i_plus)" %
                                 (elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_plus_times_density_perp",
                                 "{dot(Omega_i_plus_times_density, {1, 0, 0}), "
                                 "dot(Omega_i_plus_times_density, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_plus_times_density_para",
                                 "{0, 0, dot(Omega_i_plus_density_dependence, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_plus_times_density_para_scalar",
                                 "dot(Omega_i_plus_times_density, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_raw_plus_times_density",
                                 "n*(%e *B + %e *omega_i_raw_plus)" %
                                 (elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_raw_plus_times_density_perp",
                                 "{dot(Omega_i_raw_plus_density_dependence, {1, 0, 0}), "
                                 "dot(Omega_i_raw_plus_times_density, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_raw_plus_times_density_para",
                                 "{0, 0, dot(Omega_i_raw_plus_times_density, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_raw_plus_times_density_para_scalar",
                                 "dot(Omega_i_raw_plus_times_density, {0, 0, 1})")


    ## Omega_e density dependence
    ##
    visit.DefineVectorExpression("Omega_e_density_dependence",
                                 "n*B + cross(gradient(n), A)")
    visit.DefineVectorExpression("Omega_e_density_dependence_perp",
                                 "{dot(Omega_e_density_dependence, {1, 0, 0}), "
                                 "dot(Omega_e_density_dependence, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_e_density_dependence_para",
                                 "{0, 0, dot(Omega_e_plus_density_dependence, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_e_density_dependence_para_scalar",
                                 "dot(Omega_e_density_dependence, {0, 0, 1})")



    ## u_i_x(t) = u_i_y(t MINUS tau*0.25)
    ##
    visit.DefineVectorExpression("u_i_minus",
                                 "{u_i_x_minus, u_i_y, u_i_z}")
    visit.DefineVectorExpression("u_i_minus_perp",
                                 "{dot(u_i_minus, {1, 0, 0}), dot(u_i_minus, "
                                 "{0, 1, 0}), 0}")
    visit.DefineVectorExpression("u_i_minus_para", "{0, 0, dot(u_i_minus, {0, 0, 1})}")
    visit.DefineScalarExpression("u_i_minus_para_scalar", "dot(u_i_minus, {0, 0, 1})")

    visit.DefineVectorExpression("omega_i_minus", "{w_i_minus_x, w_i_minus_y, w_i_minus_z}")
    visit.DefineVectorExpression("omega_i_minus_perp", "{dot(omega_i_minus, {1, 0, 0}),"
                                                 "dot(omega_i_minus, {0, 1, 0}), 0})")
    visit.DefineVectorExpression("omega_i_minus_para", "{0, 0, dot(omega_i_minus, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("omega_i_minus_para_scalar", "dot(omega_i_minus, "
                                 "{0, 0, 1})")

    visit.DefineVectorExpression("omega_i_raw_minus",
                                 "{w_i_raw_x_minus, w_i_raw_y_minus, w_i_raw_z_minus}")
    visit.DefineVectorExpression("omega_i_minus_raw_perp",
                                 "{dot(omega_i_raw_minus, {1, 0, 0}),"
                                 "dot(omega_i_raw_minus, {0, 1, 0}), 0})")
    visit.DefineVectorExpression("omega_i_minus_raw_para",
                                 "{0, 0, dot(omega_i_raw_minus, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("omega_i_raw_minus_para_scalar",
                                 "dot(omega_i_raw_minus, "
                                 "{0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_minus", str(elementary_charge) +
                                 "*B +" + str(proton_mass) +
                                 "*omega_i_minus")
    visit.DefineVectorExpression("Omega_i_minus_perp", "{dot(Omega_i_minus, {1, 0, 0}), "
                                 "dot(Omega_i_minus, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_minus_para", "{0, 0, dot(Omega_i_minus, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_minus_para_scalar",
                                 "dot(Omega_i_minus, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_raw_minus",
                                 str(elementary_charge) + "*B +" + str(proton_mass) +
                                 "*omega_i_raw_minus")
    visit.DefineVectorExpression("Omega_i_raw_minus_perp", "{dot(Omega_i_raw_minus, {1, 0, 0}), "
                                 "dot(Omega_i_raw_minus, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_raw_minus_para", "{0, 0, dot(Omega_i_raw_minus, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_raw_minus_para_scalar",
                                 "dot(Omega_i_raw_minus, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_minus_density_dependence",
                                 "n*(%e *B + %e *omega_i_minus) +"
                                 " cross(gradient(n), %e *A + %e * u_i_minus)" %
                                 (elementary_charge, proton_mass,
                                  elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_minus_density_dependence_perp",
                                 "{dot(Omega_i_minus_density_dependence, {1, 0, 0}), "
                                 "dot(Omega_i_minus_density_dependence, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_minus_density_dependence_para",
                                 "{0, 0, dot(Omega_i_minus_density_dependence, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_minus_density_dependence_para_scalar",
                                 "dot(Omega_i_minus_density_dependence, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_raw_minus_density_dependence",
                                 "n*(%e *B + %e *omega_i_raw_minus) +"
                                 " cross(gradient(n), %e *A + %e * u_i_minus)" %
                                 (elementary_charge, proton_mass,
                                  elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_raw_minus_density_dependence_perp",
                                 "{dot(Omega_i_raw_minus_density_dependence, {1, 0, 0}), "
                                 "dot(Omega_i_raw_minus_density_dependence, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_raw_minus_density_dependence_para",
                                 "{0, 0, dot(Omega_i_raw_minus_density_dependence, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_raw_minus_density_dependence_para_scalar",
                                 "dot(Omega_i_raw_minus_density_dependence, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_minus_times_density",
                                 "n*(%e *B + %e *omega_i_minus)" %
                                 (elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_minus_times_density_perp",
                                 "{dot(Omega_i_plus_times_density, {1, 0, 0}), "
                                 "dot(Omega_i_minus_times_density, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_minus_times_density_para",
                                 "{0, 0, dot(Omega_i_plus_density_dependence, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_minus_times_density_para_scalar",
                                 "dot(Omega_i_minus_times_density, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_raw_minus_times_density",
                                 "n*(%e *B + %e *omega_i_raw_minus)" %
                                 (elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_raw_minus_times_density_perp",
                                 "{dot(Omega_i_raw_minus_density_dependence, {1, 0, 0}), "
                                 "dot(Omega_i_raw_minus_times_density, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_raw_minus_times_density_para",
                                 "{0, 0, dot(Omega_i_raw_minus_times_density, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_raw_minus_times_density_para_scalar",
                                 "dot(Omega_i_raw_minus_times_density, {0, 0, 1})")

    ## Canonical momentum fields
    visit.DefineVectorExpression("P_i", "%e*A + %e*u_i_plus" %
                                 (elementary_charge, proton_mass))
    visit.DefineVectorExpression("P_i_times_density", "n*P_i")

    ## Reference fields
    visit.DefineVectorExpression("B_ref", "{B_ref_x, B_ref_y, B_ref_z}")
    visit.DefineVectorExpression("A_ref", "{A_ref_x, A_ref_y, A_ref_z}")
    visit.DefineVectorExpression("u_i_ref",
                                 "{u_i_ref_x, u_i_ref_y, u_i_ref_z}")
    visit.DefineVectorExpression("omega_i_ref",
                                 "{omega_i_ref_x, omega_i_ref_y, omega_i_ref_z}")
    visit.DefineVectorExpression("u_i_ref_raw_vort",
                                 "{u_i_raw_ref_x, u_i_raw_ref_y, u_i_raw_ref_z}")
    visit.DefineVectorExpression("omega_i_ref_raw",
                                 "{w_i_raw_ref_x, w_i_raw_ref_y, w_i_raw_ref_z}")
    visit.DefineVectorExpression("P_i_ref", "%e*A_ref + %e*u_i_ref" %
                                 (elementary_charge, proton_mass))
    visit.DefineVectorExpression("P_i_ref_times_density", "n*P_i_ref")
    visit.DefineVectorExpression("P_i_ref_raw_vort", "%e*A_ref + %e*u_i_ref_raw_vort" %
                                 (elementary_charge, proton_mass))
    visit.DefineVectorExpression("P_i_ref_raw_vort_times_density", "n*P_i_ref_raw_vort")
    visit.DefineVectorExpression("Omega_i_ref", "%e*B_ref + %e*omega_i_ref_raw" %
                                 (elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_ref_times_density", "n*Omega_i_ref")
    visit.DefineVectorExpression("Omega_i_ref_density_dependence",
                                 "n*(%e *B_ref + %e *omega_i_ref) +"
                                 "cross(gradient(n), %e *A_ref + %e * u_i_ref)"
                                 % (elementary_charge,proton_mass,
                                    elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_ref_raw_vort",
                                 "%e*B_ref + %e*omega_i_ref_raw" %
                                 (elementary_charge, proton_mass))
    visit.DefineVectorExpression("Omega_i_ref_raw_vort_time_density",
                                 "n*Omega_i_ref_raw_vort")
    visit.DefineVectorExpression("Omega_i_ref_raw_vort_density_dependence",
                                 "n*(%e *B_ref + %e *omega_i_ref_raw) +"
                                 "cross(gradient(n), %e *A_ref + %e * u_i_ref_raw_vort)"
                                 % (elementary_charge,proton_mass,
                                    elementary_charge, proton_mass))


    ## Relative fields
    visit.DefineVectorExpression("B_rel", "B + B_ref")
    visit.DefineVectorExpression("A_rel", "A - A_ref")
    visit.DefineVectorExpression("u_i_rel", "u_i_plus - u_i_ref")
    visit.DefineVectorExpression("u_i_rel_raw_vort", "u_i_plus - u_i_ref_raw_vort")
    visit.DefineVectorExpression("omega_i_rel", "omega_i_plus + omega_i_ref")
    visit.DefineVectorExpression("omega_i_rel_raw", "omega_i_raw_plus + omega_i_ref_raw")
    visit.DefineVectorExpression("P_i_rel", "P_i - P_i_ref")
    visit.DefineVectorExpression("P_i_rel_times_density",
                                 "P_i_times_density - P_i_ref_times_density")
    visit.DefineVectorExpression("P_i_rel_raw_vort", "P_i - P_i_ref_raw_vort")
    visit.DefineVectorExpression("P_i_rel_raw_vort_times_density",
                                 "P_i_rel_raw_vort_times_density - P_i_ref_raw_vort_times_density")
    visit.DefineVectorExpression("Omega_i_rel", "Omega_i_plus + Omega_i_ref")
    visit.DefineVectorExpression("Omega_i_rel_times_density",
                                 "Omega_i_plus_times_density + Omega_i_ref_times_density")
    visit.DefineVectorExpression("Omega_i_rel_density_dependence",
                                 "Omega_i_plus_density_dependence + Omega_i_ref_density_dependence")
    visit.DefineVectorExpression("Omega_i_raw_rel",
                                 "Omega_i_raw_plus + Omega_i_ref_raw_vort")
    visit.DefineVectorExpression("Omega_i_raw_rel_times_density",
                                 "Omega_i_raw_plus_times_density +"
                                 "Omega_i_ref_raw_vort_times_density")
    visit.DefineVectorExpression("Omega_i_ref_raw_density_dependence",
                                 "Omega_i_raw_plus_density_dependence +"
                                 "Omega_i_ref_raw_vort_density_dependence")
=======
    visit.DefineVectorExpression("u_e_norm", "{u_e_norm_x, u_e_norm_y, "
                                 "u_e_norm_z}")
    visit.DefineVectorExpression("u_e", "%1.2e * u_e_norm" % alpha)
    visit.DefineVectorExpression("u_e_perp", "{dot(u_e, {1, 0, 0}), "
                                 "dot(u_e, {0, 1, 0}), 0}")
    visit.DefineScalarExpression("u_e_para_scalar", "dot(u_e, {0, 0, 1})")

    visit.DefineVectorExpression("u_e_fitted_alpha_z03", "{u_e_x_fitted_alpha_z03, " +
                                                          "u_e_y_fitted_alpha_z03, " +
                                                          "u_e_z_fitted_alpha_z03}")
    visit.DefineVectorExpression("u_e_fitted_alpha_z04", "{u_e_x_fitted_alpha_z04, "
                                                          "u_e_y_fitted_alpha_z04, "
                                                          "u_e_z_fitted_alpha_z04}")
    visit.DefineVectorExpression("u_e_fitted_alpha_both_planes",
                                 "{u_e_x_fitted_alpha_both_planes, "
                                 "u_e_y_fitted_alpha_both_planes, "
                                 "u_e_z_fitted_alpha_both_planes}")

    visit.DefineVectorExpression("u_i_term1", "{u_i_term1_x, u_i_term1_y, "
                                 "u_i_term1_z}")
    visit.DefineVectorExpression("u_i_term1_perp", "{u_i_term1_x, "
                                 "u_i_term1_y, 0}")
    visit.DefineVectorExpression("u_i_term1_para", "{0, 0, u_i_term1_z}")
    visit.DefineScalarExpression("u_i_term1_para_scalar", "u_i_term1_z")

    visit.DefineVectorExpression("omega_i_term1", "{w_i_term1_x, w_i_term1_y, "
                                 "w_i_term1_z}")
    visit.DefineVectorExpression("omega_i_term1_perp", "{w_i_term1_x, "
                                 "w_i_term1_y, 0}")
    visit.DefineVectorExpression("omega_i_term1_para", "{0, 0, w_i_term1_z}")
    visit.DefineScalarExpression("omega_i_term1_para_scalar", "w_i_term1_z")

    visit.DefineVectorExpression("omega_i_term2", "{w_i_term2_x, w_i_term2_y, "
                                 "w_i_term2_z} * %1.2e" % alpha)
    visit.DefineVectorExpression("omega_i_term2_perp", "{w_i_term2_x, "
                                 "w_i_term2_y, 0} * %1.2e" % alpha)
    visit.DefineVectorExpression("omega_i_term2_para", "{0, 0, "
                                 "w_i_term2_z} * %1.2e" % alpha)
    visit.DefineScalarExpression("omega_i_term2_para_scalar",
                                 "w_i_term2_z * %1.2e" % alpha)

    visit.DefineVectorExpression("omega_i_term2_fitted_alpha_z03", "{w_i_term2_x_fitted_alpha_z03,"
                                 "w_i_term2_y_fitted_alpha_z03, w_i_term2_z_fitted_alpha_z03}")
    visit.DefineVectorExpression("omega_i_term2_fitted_alpha_z04", "{w_i_term2_x_fitted_alpha_z04,"
                                 "w_i_term2_y_fitted_alpha_z04, w_i_term2_z_fitted_alpha_z04}")
    visit.DefineVectorExpression("omega_i_term2_fitted_alpha_both_planes",
                                 "{w_i_term2_x_fitted_alpha_both_planes, "
                                 "w_i_term2_y_fitted_alpha_both_planes, "
                                 "w_i_term2_z_fitted_alpha_both_planes}")


    visit.DefineVectorExpression("u_i", "u_i_term1 + u_e")
    visit.DefineVectorExpression("u_i_perp", "{dot(u_i, {1, 0, 0}), dot(u_i, "
                                 "{0, 1, 0}), 0}")
    visit.DefineVectorExpression("u_i_para", "{0, 0, dot(u_i, {0, 0, 1})}")
    visit.DefineScalarExpression("u_i_para_scalar", "dot(u_i, {0, 0, 1})")


    visit.DefineVectorExpression("u_i_fitted_alpha_z03", "u_i_term1 + u_e_fitted_alpha_z03")
    visit.DefineVectorExpression("u_i_fitted_alpha_z04", "u_i_term1 + u_e_fitted_alpha_z04")
    visit.DefineVectorExpression("u_i_fitted_alpha_both_planes", "u_i_term1 + "
                                 "u_e_fitted_alpha_both_planes")

    visit.DefineVectorExpression("omega_i", "omega_i_term1 + "
                                 "omega_i_term2")
    visit.DefineVectorExpression("omega_i_perp", "{dot(omega_i, {1, 0, 0}), "
                                 "dot(omega_i, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("omega_i_para", "{0, 0, dot(omega_i, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("omega_i_para_scalar", "dot(omega_i, "
                                 "{0, 0, 1})")

    visit.DefineVectorExpression("omega_i_fitted_alpha_z03", "omega_i_term1 + "
                                 "omega_i_term2_fitted_alpha_z03")
    visit.DefineVectorExpression("omega_i_fitted_alpha_z04", "omega_i_term1 + "
                                 "omega_i_term2_fitted_alpha_z04")
    visit.DefineVectorExpression("omega_i_fitted_alpha_both_planes", "omega_i_term1 + "
                                 "omega_i_term2_fitted_alpha_both_planes")

    visit.DefineVectorExpression("Omega_i", "B +" + str(proton_mass) + "/" +
                                 str(elementary_charge) + "*omega_i")
    visit.DefineVectorExpression("Omega_i_perp", "{dot(Omega_i, {1, 0, 0}), "
                                 "dot(Omega_i, {0, 1, 0}), 0}")
    visit.DefineVectorExpression("Omega_i_para", "{0, 0, dot(Omega_i, "
                                 "{0, 0, 1})}")
    visit.DefineScalarExpression("Omega_i_para_scalar",
                                 "dot(Omega_i, {0, 0, 1})")

    visit.DefineVectorExpression("Omega_i_fitted_alpha_z03", "B +" + str(proton_mass) + "/" +
                                 str(elementary_charge) + "*omega_i_fitted_alpha_z03")
    visit.DefineVectorExpression("Omega_i_fitted_alpha_z04", "B +" + str(proton_mass) + "/" +
                                 str(elementary_charge) + "*omega_i_fitted_alpha_z03")
    visit.DefineVectorExpression("Omega_i_fitted_alpha_both_planes", "B +" + str(proton_mass) + "/" +
                                 str(elementary_charge) + "*omega_i_fitted_alpha_both_planes")
>>>>>>> canonical_flux_tubes/master


def normalize_scalar(visit, scalar_name,
                     normalized_scalar_name):
    r"""
    Determine max of scalar.
    """
    visit.AddPlot("Pseudocolor", scalar_name)
    visit.DrawPlots()
    visit.Query("Max")
    max = visit.GetQueryOutputValue()
    visit.DefineScalarExpression(normalized_scalar_name,
                                 "%s - %g" % (scalar_name, max))
    visit.DeleteActivePlots()


def launch_points_inner_outer(center, plane=0.249,
<<<<<<< HEAD
                              radius_inner=0.001, radius_outer=0.005,
=======
                              radius_inner=0.005, radius_outer=0.01,
>>>>>>> canonical_flux_tubes/master
                              num_inner=80, num_outer=60):
    r"""
    Calculate points on a circle outline for a given center point.
    """
    thetas = circle_with_cut_thetas(num_outer)
    points_outer = launch_points(center, thetas, radius=radius_outer,
                                 plane=plane)

    thetas = full_circle_thetas(num_inner)
    points_inner = launch_points(center, thetas, radius=radius_inner,
                                 plane=plane)

    return points_outer, points_inner


def full_circle_thetas(num_points):
    r"""
    """
    thetas = np.linspace(0, 2.*np.pi, num_points)
    return thetas


def circle_with_cut_thetas(num_points):
    r"""
    """
    thetas = np.linspace(0, 3./4.*np.pi, num_points)
    thetas = np.concatenate((thetas, np.linspace(5./4.*np.pi, 2.*np.pi,
                                                 num_points)))
    return thetas


<<<<<<< HEAD
def launch_points(center, thetas, radius=0.003,
=======
def launch_points(center, thetas, radius=0.03,
>>>>>>> canonical_flux_tubes/master
                  plane=0.249):
    r"""
    """
    x_points = radius * np.cos(thetas) + center[0]
    y_points = radius * np.sin(thetas) + center[1]
    z_points = plane * np.ones(x_points.size)
    points = np.empty((x_points.size + y_points.size + z_points.size))
    points[0::3] = x_points
    points[1::3] = y_points
    points[2::3] = z_points
    points = tuple(points)
    return points


<<<<<<< HEAD
def setup_scalar_isosurface(visit, quantity,
                            colortable="PuRd", max_val=1., min_val=0.9):
    r"""
    Setup iso_surface. Works best if quantity is plane normalized.
    """
    visit.AddPlot("Pseudocolor", quantity, 1, 0)
    PseudocolorAtts = visit.PseudocolorAttributes()
    PseudocolorAtts.colorTableName = colortable
    PseudocolorAtts.opacityType = PseudocolorAtts.Constant
    PseudocolorAtts.opacity = 0.25
    PseudocolorAtts.smoothingLevel = 0
    PseudocolorAtts.legendFlag = 0
    visit.SetPlotOptions(PseudocolorAtts)

    visit.AddOperator("Isosurface", 0)
    IsosurfaceAtts = visit.IsosurfaceAttributes()
    IsosurfaceAtts.contourNLevels = 5
    IsosurfaceAtts.contourValue = ()
    IsosurfaceAtts.contourPercent = ()
    IsosurfaceAtts.contourMethod = IsosurfaceAtts.Level
    IsosurfaceAtts.minFlag = 1
    IsosurfaceAtts.min = min_val
    IsosurfaceAtts.maxFlag = 1
    IsosurfaceAtts.max = max_val
    visit.SetOperatorOptions(IsosurfaceAtts, 0)
    return PseudocolorAtts, IsosurfaceAtts


=======
>>>>>>> canonical_flux_tubes/master
def setup_current_pseudocolor(visit, colortable="Greens", max_val=None,
                              min_val=None):
    r"""
    Setup pseudocolor current plot.
    """
<<<<<<< HEAD
    visit.AddPlot("Pseudocolor", "J_smooth_mag", 1, 0)
=======
    visit.AddPlot("Pseudocolor", "J_raw_mag", 1, 0)
>>>>>>> canonical_flux_tubes/master
    PseudocolorAtts = visit.PseudocolorAttributes()
    PseudocolorAtts.scaling = PseudocolorAtts.Linear
    PseudocolorAtts.limitsMode = PseudocolorAtts.OriginalData
    PseudocolorAtts.colorTableName = colortable

    if max_val:
        PseudocolorAtts.maxFlag = 1
        PseudocolorAtts.max = max_val
    if min_val:
        PseudocolorAtts.minFlag = 1
        PseudocolorAtts.min = min_val

    visit.SetPlotOptions(PseudocolorAtts)
    visit.AddOperator("Slice", 0)
    SliceAtts = visit.SliceAttributes()
    SliceAtts.originType = SliceAtts.Intercept
    SliceAtts.originIntercept = 0.249
    SliceAtts.axisType = SliceAtts.ZAxis
    SliceAtts.project2d = 0
    visit.SetOperatorOptions(SliceAtts, 0)

    return PseudocolorAtts, SliceAtts


def setup_massless_electron_canonical_flux_tubes(visit, points_outer,
                                                 points_inner):
    r"""
    Setup two massless electron canonical flux tubes i.e. magnetic flux tubes.
    Intended to be inner and outer flux tubes.
    """
    visit.AddPlot("Streamline", "B", 1, 0)
    StreamlineAtts_outer = visit.StreamlineAttributes()
    StreamlineAtts_outer.sourceType = StreamlineAtts_outer.SpecifiedPointList
    StreamlineAtts_outer.SetPointList(points_outer)
    StreamlineAtts_outer.coloringMethod = StreamlineAtts_outer.Solid
    StreamlineAtts_outer.colorTableName = "Default"
    StreamlineAtts_outer.singleColor = (255, 0, 0, 255)
    StreamlineAtts_outer.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_outer)

    visit.AddPlot("Streamline", "B", 1, 0)
    StreamlineAtts_inner = visit.StreamlineAttributes()
    StreamlineAtts_inner.sourceType = StreamlineAtts_inner.SpecifiedPointList
    StreamlineAtts_inner.SetPointList(points_inner)
    StreamlineAtts_inner.coloringMethod = StreamlineAtts_inner.Solid
    StreamlineAtts_inner.colorTableName = "Default"
    StreamlineAtts_inner.singleColor = (190, 64, 0, 255)
    StreamlineAtts_inner.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_inner)

    return StreamlineAtts_outer, StreamlineAtts_inner


<<<<<<< HEAD
def setup_outer_inner_ion_canonical_flux_tubes(visit, quantity, points_outer,
                                               points_inner,
                                               outer_color=dark_grey,
                                               inner_color=black):
=======
def setup_inner_outer_ion_canonical_flux_tubes(visit, points_outer,
                                               points_inner,
                                               outer_color=tan,
                                               inner_color=olive):
>>>>>>> canonical_flux_tubes/master
    r"""
    Setup two ion canonical flux tubes.
    Inteteded to be inner and outer flux tubes.
    """
<<<<<<< HEAD
    visit.AddPlot("Streamline", quantity, 1, 0)
=======
    visit.AddPlot("Streamline", "Omega_i", 1, 0)
>>>>>>> canonical_flux_tubes/master
    StreamlineAtts_outer = visit.StreamlineAttributes()
    StreamlineAtts_outer.sourceType = StreamlineAtts_outer.SpecifiedPointList
    StreamlineAtts_outer.SetPointList(points_outer)
    StreamlineAtts_outer.coloringMethod = StreamlineAtts_outer.Solid
    StreamlineAtts_outer.colorTableName = "Default"
    StreamlineAtts_outer.singleColor = outer_color
    StreamlineAtts_outer.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_outer)

<<<<<<< HEAD
    visit.AddPlot("Streamline", quantity, 1, 0)
=======
    visit.AddPlot("Streamline", "Omega_i", 1, 0)
>>>>>>> canonical_flux_tubes/master
    StreamlineAtts_inner = visit.StreamlineAttributes()
    StreamlineAtts_inner.sourceType = StreamlineAtts_inner.SpecifiedPointList
    StreamlineAtts_inner.SetPointList(points_inner)
    StreamlineAtts_inner.coloringMethod = StreamlineAtts_inner.Solid
    StreamlineAtts_inner.colorTableName = "Default"
    StreamlineAtts_inner.singleColor = inner_color
    StreamlineAtts_inner.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_inner)

    return StreamlineAtts_outer, StreamlineAtts_inner


def setup_forward_backward_ion_canonical_flux_tubes(visit, points_foward,
                                                    points_backward,
                                                    forward_color=tan,
                                                    backward_color=olive):
    r"""
    Setup two ion canonical flux tubes, one integrating in the forward
    direction, one integrating in the backward direction.
    """
    visit.AddPlot("Streamline", "Omega_i", 1, 0)

    StreamlineAtts_forward = visit.StreamlineAttributes()
    StreamlineAtts_forward.sourceType = StreamlineAtts_forward.SpecifiedPointList
    StreamlineAtts_forward.SetPointList(points_foward)
    StreamlineAtts_forward.coloringMethod = StreamlineAtts_forward.Solid
    StreamlineAtts_forward.colorTableName = "Default"
    StreamlineAtts_forward.singleColor = forward_color
    StreamlineAtts_forward.integrationDirection = StreamlineAtts_forward.Forward
    StreamlineAtts_forward.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_forward)

    visit.AddPlot("Streamline", "Omega_i", 1, 0)
    StreamlineAtts_backward = visit.StreamlineAttributes()
    StreamlineAtts_backward.sourceType = StreamlineAtts_backward.SpecifiedPointList
    StreamlineAtts_backward.SetPointList(points_backward)
    StreamlineAtts_backward.coloringMethod = StreamlineAtts_backward.Solid
    StreamlineAtts_backward.colorTableName = "Default"
    StreamlineAtts_backward.singleColor = backward_color
    StreamlineAtts_backward.integrationDirection = StreamlineAtts_backward.Backward
    StreamlineAtts_backward.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_backward)

    return StreamlineAtts_forward, StreamlineAtts_backward


def setup_backward_and_B_stream(visit, name, launch_points,
                                B_launch_points, color=green, B_color=red):
    r"""
    """
    visit.AddPlot("Streamline", 'B', 1, 0)
    StreamlineAtts_B = visit.StreamlineAttributes()
    StreamlineAtts_B.sourceType = StreamlineAtts_B.SpecifiedPointList
    StreamlineAtts_B.SetPointList(B_launch_points)
    StreamlineAtts_B.coloringMethod = StreamlineAtts_B.Solid
    StreamlineAtts_B.colorTableName = "Default"
    StreamlineAtts_B.singleColor = B_color
    StreamlineAtts_B.integrationDirection = StreamlineAtts_B.Forward
    StreamlineAtts_B.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_B)


    visit.AddPlot("Streamline", name, 1, 0)
    StreamlineAtts_backward = visit.StreamlineAttributes()
    StreamlineAtts_backward.sourceType = StreamlineAtts_backward.SpecifiedPointList
    StreamlineAtts_backward.SetPointList(launch_points)
    StreamlineAtts_backward.coloringMethod = StreamlineAtts_backward.Solid
    StreamlineAtts_backward.colorTableName = "Default"
    StreamlineAtts_backward.singleColor = color
    StreamlineAtts_backward.integrationDirection = StreamlineAtts_backward.Backward
    StreamlineAtts_backward.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_backward)

    return StreamlineAtts_B, StreamlineAtts_backward


def setup_forward_backward_alpha_fitted_ion_canonical_flux_tubes(visit, points_foward,
<<<<<<< HEAD
                                                                 points_backward,
                                                                 forward_color=tan,
                                                                 backward_color=olive):
=======
                                                    points_backward,
                                                    forward_color=tan,
                                                    backward_color=olive):
>>>>>>> canonical_flux_tubes/master
    r"""
    Setup two ion canonical flux tubes, one integrating in the forward
    direction, one integrating in the backward direction.
    """
    visit.AddPlot("Streamline", "Omega_i_alpha_fitted", 1, 0)
    StreamlineAtts_forward = visit.StreamlineAttributes()
    StreamlineAtts_forward.sourceType = StreamlineAtts_forward.SpecifiedPointList
    StreamlineAtts_forward.SetPointList(points_foward)
    StreamlineAtts_forward.coloringMethod = StreamlineAtts_forward.Solid
    StreamlineAtts_forward.colorTableName = "Default"
    StreamlineAtts_forward.singleColor = forward_color
    StreamlineAtts_forward.integrationDirection = StreamlineAtts_forward.Forward
    StreamlineAtts_forward.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_forward)

    visit.AddPlot("Streamline", "Omega_i_alpha_fitted", 1, 0)
    StreamlineAtts_backward = visit.StreamlineAttributes()
    StreamlineAtts_backward.sourceType = StreamlineAtts_backward.SpecifiedPointList
    StreamlineAtts_backward.SetPointList(points_backward)
    StreamlineAtts_backward.coloringMethod = StreamlineAtts_backward.Solid
    StreamlineAtts_backward.colorTableName = "Default"
    StreamlineAtts_backward.singleColor = backward_color
    StreamlineAtts_backward.integrationDirection = StreamlineAtts_backward.Backward
    StreamlineAtts_backward.legendFlag = 0
    visit.SetPlotOptions(StreamlineAtts_backward)

    return StreamlineAtts_forward, StreamlineAtts_backward


<<<<<<< HEAD
def setup_field_line(visit, quantity,
                     launch_point=(0.01, 0.01), launch_z=0.249,
                     color=black):
    r"""
    Setup single field line plot to better see twistedness.
    """
    visit.AddPlot("Streamline", quantity, 1, 0)
    StreamlineAtts_line = visit.StreamlineAttributes()
    StreamlineAtts_line.sourceType = StreamlineAtts_line.SpecifiedPoint
    StreamlineAtts_line.pointSource = (launch_point[0], launch_point[1], launch_z)
    StreamlineAtts_line.coloringMethod = StreamlineAtts_line.Solid
    StreamlineAtts_line.singleColor = color
=======
def setup_field_line(visit, center=(0.01, 0.01, 0.249),
                     outer_radius=0.01):
    r"""
    Setup single field line plot to better see twistedness.
    """
    visit.AddPlot("Streamline", "B", 1, 0)
    StreamlineAtts_line = visit.StreamlineAttributes()
    StreamlineAtts_line.sourceType = StreamlineAtts_line.SpecifiedPoint
    StreamlineAtts_line.pointSource = (center[0], center[1] + outer_radius, center[2])
    StreamlineAtts_line.coloringMethod = StreamlineAtts_line.Solid
    StreamlineAtts_line.singleColor = (255, 255, 153, 255)
>>>>>>> canonical_flux_tubes/master
    StreamlineAtts_line.legendFlag = 0
    StreamlineAtts_line.showSeeds = 0
    visit.SetPlotOptions(StreamlineAtts_line)
    return StreamlineAtts_line


def setup_annotations(visit, time_scale=1):
    r"""
    Setup Annotations: scale tick font size, label font size,
    hide unecessary text.
    """
    AnnotationAtts = visit.AnnotationAttributes()
    AnnotationAtts.axes3D.autoSetScaling = 0
    AnnotationAtts.axes3D.xAxis.title.visible = 0
    AnnotationAtts.axes3D.yAxis.title.visible = 0
    AnnotationAtts.axes3D.zAxis.title.visible = 0
    AnnotationAtts.axes3D.xAxis.label.font.scale = 3
    AnnotationAtts.axes3D.xAxis.label.scaling = -2
    AnnotationAtts.axes3D.yAxis.label.font.scale = 3
    AnnotationAtts.axes3D.yAxis.label.scaling = -2
    AnnotationAtts.axes3D.zAxis.label.font.scale = 3
    AnnotationAtts.axes3D.zAxis.label.scaling = -2
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.databaseInfoTimeScale = time_scale
    AnnotationAtts.databaseInfoTimeOffset = 0
    visit.SetAnnotationAttributes(AnnotationAtts)
    return AnnotationAtts


def set_default_view(visit):
    r"""
    Set default view for viewing fluxtubes.
    If view needs to be modified it is best to align visit with gui and save
    parameters from visit.GetView3D().
    """
    view = visit.GetView3D()
    view.SetViewNormal((-0.731293, 0.40847, 0.546227))
    view.SetFocus((0.00202222, 0.000976744, 0.331997))
    view.SetViewUp((0.322268, 0.91274, -0.251095))
    view.SetViewAngle(30)
    view.SetParallelScale(0.088383)
    view.SetNearPlane(-0.176766)
    view.SetImagePan((0, 0))
    view.SetImageZoom(1)
    view.SetPerspective(1)
    view.SetEyeAngle(2)
    view.SetCenterOfRotationSet(0)
    view.SetCenterOfRotation((0.00202222, 0.000976744, 0.331997))
    view.SetAxis3DScaleFlag(0)
    view.SetAxis3DScales((1, 1, 1))
    view.SetShear((0, 0, 1))
    view.SetWindowValid(0)
    visit.SetView3D(view)


<<<<<<< HEAD
def set_default_view_lower_angle(visit):
    r"""
    """
    view = visit.GetView3D()
    view.setViewNormal((-0.776189, 0.193398, 0.600106))
    view.setFocus((0.00202222, 0.000976744, 0.331997))
    view.setViewUp((0.138771, 0.980856, -0.136615))
    view.setViewAngle(30)
    view.setParallelScale(0.088383)
    view.setNearPlane(-0.176766)
    view.setFarPlane(0.175437)
    view.setImagePan((0, 0))
    view.setImageZoom(1)
    view.setPerspective(1)
    view.setEyeAngle(2)
    view.setCenterOfRotationSet(0)
    view.setCenterOfRotation((0.00202222, 0.000976744, 0.331997))
    view.setAxis3DScaleFlag(0)
    view.setAxis3DScales((1, 1, 1))
    view.setShear((0, 0, 1))
    view.setWindowValid(0)
    SetView3D(view)


def set_positive_x_view(visit):
    r"""
    """
    view = visit.GetView3D()
    view.SetViewNormal((-0.00997631, 0.0600385, 0.0335938))
    view.SetFocus((0.00202222, -0.00202703, 0.331997))
    view.SetViewUp((0.0598395, 0.998184, 0.00689852))
    view.SetViewAngle(30)
    view.SetParallelScale(0.0877186)
    view.SetNearPlane(-0.175437)
    view.SetImagePan((0, 0))
    view.SetImageZoom(1)
    view.SetPerspective(1)
    view.SetEyeAngle(2)
    view.SetCenterOfRotationSet(0)
    view.SetCenterOfRotation((0.00202222, -0.00202703, 0.331997))
    view.SetAxis3DScaleFlag(0)
    view.SetAxis3DScales((1, 1, 1))
    view.SetShear((0, 0, 1))
    view.SetWindowValid(1)
    visit.SetView3D(view)


def set_postive_z_view(visit):
    r"""
    """
    view = visit.GetView3D()
    view.SetViewNormal((0.00944856, 0.0379894, 0.999233))
    view.SetFocus((0.00202222, -0.00202703, 0.331997))
    view.SetViewUp((-0.00367716, 0.9999274, 0.0037961))
    view.SetViewAngle(30)
    view.SetParallelScale(0.0877186)
    view.SetNearPlane(-0.175437)
    view.SetImagePan((0, 0))
    view.SetImageZoom(2.14359)
    view.SetPerspective(1)
    view.SetEyeAngle(2)
    view.SetCenterOfRotationSet(0)
    view.SetCenterOfRotation((0.00202222, -0.00202703, 0.331997))
    view.SetAxis3DScaleFlag(0)
    view.SetAxis3DScales((1, 1, 1))
    view.SetShear((0, 0, 1))
    view.SetWindowValid(1)
    visit.SetView3D(view)


def set_negative_z_view(visit):
    r"""
    """
    view = visit.GetView3D()
    view.SetViewNormal((-0.00894299, -0.00985814, 0.999911))
    view.SetFocus((0.00202222, 0.000976744, 0.331997))
    view.SetViewUp((0.00367716, 0.999944, 0.00989136))
    view.SetViewAngle(30)
    view.SetParallelScale(0.0877186)
    view.SetNearPlane(-0.175437)
    view.SetImagePan((0, 0))
    view.SetImageZoom(2.14359)
    view.SetPerspective(1)
    view.SetEyeAngle(2)
    view.SetCenterOfRotationSet(0)
    view.SetCenterOfRotation((0.00202222, -0.00202703, 0.331997))
    view.SetAxis3DScaleFlag(0)
    view.SetAxis3DScales((1, 1, 1))
    view.SetShear((0, 0, 1))
    view.SetWindowValid(1)
    visit.SetView3D(view)


=======
>>>>>>> canonical_flux_tubes/master
def determine_j_mag_extrema(database_path, plane_num=0):
    r"""
    Determine extrema over time of current across all shots.
    Can be sued to set min and max values for colormaps.
    """
    numpy_archives =  glob(database_path + '*.npz')
    data = np.load(numpy_archives[0])
    j_x = data['j_x'][:, :, plane_num]
    j_y = data['j_y'][:, :, plane_num]
    j_z = data['j_z'][:, :, plane_num]
    j_mag = np.sqrt(j_x**2. + j_y**2. + j_z**2.)
    j_mag_max = np.nanmax(j_mag)
    j_mag_min = np.nanmin(j_mag)
    for archive in numpy_archives[1:]:
        data = np.load(archive)
        j_x = data['j_x'][:, :, plane_num]
        j_y = data['j_y'][:, :, plane_num]
        j_z = data['j_z'][:, :, plane_num]
        j_mag = np.sqrt(j_x**2. + j_y**2. + j_z**2.)
        j_mag_max = np.nanmax(j_mag) if np.nanmax(j_mag) > j_mag_max else j_mag_max
        j_mag_min = np.nanmin(j_mag) if np.nanmin(j_mag) < j_mag_min else j_mag_min
    return j_mag_max, j_mag_min


def set_save_settings(visit):
    r"""
    Set and return save_atts.
    """
    save_atts = visit.SaveWindowAttributes()
    save_atts.format = save_atts.PNG
    save_atts.height = 1080
    save_atts.width = 1920
    save_atts.family = 0
    visit.SetSaveWindowAttributes(save_atts)
    return save_atts

def main():
    r"""
    """
    args = parse_args()
    database_prefix = args.database_prefix + args.database_date
    visit.Launch()
<<<<<<< HEAD
    today = datetime.now().strftime('%Y-%m-%d-%H-%M')
=======
    today = str(date.today())
>>>>>>> canonical_flux_tubes/master
    out_dir = '../output/' + today
    try:
       os.makedirs(out_dir)
    except:
        pass


    if args.interactive_session:
       visit.OpenDatabase(database_prefix + args.database_postfix)
       define_expressions(visit, args.alpha_constant)
       visit.OpenGUI()
       return

    output_path = out_dir + '/' + args.output_prefix
    visit.OpenDatabase(database_prefix + args.database_postfix)
    define_expressions(visit, args.alpha_constant)
    field_nulls = np.loadtxt(args.field_nulls)
<<<<<<< HEAD
    AnnotationAtts = setup_annotations(visit, time_scale=args.time_scale)

    plot_count = 0

    if args.current_plane:
        PseudocolorAtts, SliceAtts = setup_current_pseudocolor(visit, max_val=args.current_max,
                                                               min_val=args.current_min)
        plot_count += 1
    if args.temperature_tubes:
        setup_scalar_isosurface(visit, "Te_plane_normalized", colortable="PuRd")
        plot_count += 1
    if args.density_tubes:
        setup_scalar_isosurface(visit, "n_plane_normalized", colortable="Greys")
        plot_count += 1

    points_outer, points_inner = launch_points_inner_outer(field_nulls[0])
    if args.ion:
        (StreamlineAtts_ion_outer,
         StreamlineAtts_ion_inner) = setup_outer_inner_ion_canonical_flux_tubes(visit,
                                                                                args.omega_to_use,
                                                                                points_outer,
                                                                                points_inner)
        plot_count += 2
    if args.electron:
        stream_line_func = setup_massless_electron_canonical_flux_tubes
        (StreamlineAtts_electron_outer,
         StreamlineAtts_electron_inner) = setup_massless_electron_canonical_flux_tubes(visit,
                                                                                       points_outer,
                                                                                       points_inner)
        plot_count += 2

    if args.velocity:
        (velocity_stream_1,
         velocity_stream_2) = setup_outer_inner_ion_canonical_flux_tubes(visit,
                                                                         'u_i_plus',
                                                                         points_inner,
                                                                         points_outer,
                                                                         outer_color=aqua,
                                                                         inner_color=navy)
        plot_count += 2
    #elif args.current:
    #    stream_line_func = setup_backward_and_B_stream
    #    current_thetas = full_circle_thetas(20)
    #    B_thetas = full_circle_thetas(10)
    #    current_launch_points = launch_points(field_nulls[0],
    #                                          current_thetas,
    #                                          radius=0.001)
    #    B_launch_points = launch_points(field_nulls[0],
    #                                    B_thetas,
    #                                    radius=0.005)
    #    params = {'visit': visit, 'name': args.current_to_use,
    #              'launch_points': current_launch_points,
    #              'B_launch_points': B_launch_points}

    if args.view == 'default':
        set_default_view(visit)
    elif args.view == 'default_lower_angle':
        set_default_view_lower_angle(visit)
    elif args.view == 'positive_z':
        set_postive_z_view(visit)
    elif args.view == 'negative_z':
        set_negative_z_view(visit)
    elif args.view == 'positive_x':
        set_positive_x_view(visit)


    if args.double_stream:
        stream_launch_point = (field_nulls[args.start_time_point][0] + args.x_offset,
                               field_nulls[args.start_time_point][1] + args.y_offset)
        setup_field_line(visit, args.double_stream,
                         launch_point=(stream_launch_point[0] + 0.001,
                                       stream_launch_point[1] + 0.))
        setup_field_line(visit, args.double_stream,
                         launch_point=(stream_launch_point[0] + 0.,
                                       stream_launch_point[1] + 0.001))
        setup_field_line(visit, args.double_stream,
                         launch_point=(stream_launch_point[0] - 0.001,
                                       stream_launch_point[1] + 0.))
        setup_field_line(visit, args.double_stream,
                         launch_point=(stream_launch_point[0] + 0.,
                                       stream_launch_point[1] - 0.001))

        setup_field_line(visit, args.double_stream,
                         launch_point=(stream_launch_point[0] + 0.005,
                                       stream_launch_point[1] + 0.))
        setup_field_line(visit, args.double_stream,
                         launch_point=(stream_launch_point[0] + 0.,
                                       stream_launch_point[1] + 0.005))
        setup_field_line(visit, args.double_stream,
                         launch_point=(stream_launch_point[0] - 0.005,
                                       stream_launch_point[1] + 0.))
        setup_field_line(visit, args.double_stream,
                         launch_point=(stream_launch_point[0] + 0.,
                                       stream_launch_point[1] - 0.005))

        setup_field_line(visit, args.double_stream,
                         launch_point=stream_launch_point)

        setup_field_line(visit, 'B',
                         launch_point=stream_launch_point, color=red)


    visit.DrawPlots()
    save_atts = set_save_settings(visit)
    ending = '.png'
    visit.SetTimeSliderState(args.start_time_point)
    if args.wait_for_manual_settings:
        visit.OpenGUI()
        comment = raw_input()

    for time_point in xrange(args.start_time_point, args.end_time_point):
        print time_point
        plot_number = 0
        save_atts.fileName = output_path + str(time_point + args.file_number_offset).zfill(4) + ending
        visit.SetSaveWindowAttributes(save_atts)

        if args.current_plane:
            plot_number += 1
        if args.temperature_tubes:
            plot_number += 1
        if args.density_tubes:
            plot_number += 1

        points_outer, points_inner = launch_points_inner_outer(field_nulls[time_point])

        if args.ion:
            visit.SetActivePlots(plot_number)
            StreamlineAtts_ion_outer.SetPointList(points_outer)
            visit.SetPlotOptions(StreamlineAtts_ion_outer)
            plot_number += 1

            visit.SetActivePlots(plot_number)
            StreamlineAtts_ion_inner.SetPointList(points_inner)
            visit.SetPlotOptions(StreamlineAtts_ion_inner)
            plot_number += 1

        if args.electron:
            visit.SetActivePlots(plot_number)
            StreamlineAtts_electron_outer.SetPointList(points_outer)
            visit.SetPlotOptions(StreamlineAtts_electron_outer)
            plot_number += 1

            visit.SetActivePlots(plot_number)
            StreamlineAtts_electron_inner.SetPointList(points_inner)
            visit.SetPlotOptions(StreamlineAtts_electron_inner)
            plot_number +=1

        if args.velocity:
            visit.SetActivePlots(plot_number)
            velocity_stream_1.SetPointList(points_outer)
            visit.SetPlotOptions(velocity_stream_1)
            plot_number +=1

            visit.SetActivePlots(plot_number)
            velocity_stream_2.SetPointList(points_inner)
            visit.SetPlotOptions(velocity_stream_2)
            plot_number +=1
=======
    PseudocolorAtts, SliceAtts = setup_current_pseudocolor(visit, max_val=args.current_max, min_val=args.current_min)
    points_outer, points_inner = launch_points_inner_outer(field_nulls[0])
    AnnotationAtts = setup_annotations(visit, time_scale=args.time_scale)

    if args.ion_forward_backward:
        points_outer, points_inner = launch_points_inner_outer(field_nulls[0])
        stream_line_func = setup_forward_backward_ion_canonical_flux_tubes
    elif args.ion:
        points_outer, points_inner = launch_points_inner_outer(field_nulls[0])
        stream_line_func = setup_inner_outer_ion_canonical_flux_tubes
    elif args.electron:
        points_outer, points_inner = launch_points_inner_outer(field_nulls[0])
        stream_line_func = setup_massless_electron_canonical_flux_tubes
        points_outer = points_inner
    elif args.current:
        stream_line_func = setup_backward_and_B_stream
        current_thetas = full_circle_thetas(20)
        B_thetas = full_circle_thetas(10)
        current_launch_points = launch_points(field_nulls[0],
                                              current_thetas,
                                              radius=0.001)
        B_launch_points = launch_points(field_nulls[0],
                                        B_thetas,
                                        radius=0.005)
        params = {'visit': visit, 'name': args.current_to_use,
                  'launch_points': current_launch_points,
                  'B_launch_points': B_launch_points}

    (StreamlineAtts_flux_1,
     StreamlineAtts_flux_2) = stream_line_func(**params)

    set_default_view(visit)
    visit.DrawPlots()
    save_atts = set_save_settings(visit)
    ending = '.png'

    for time_point in xrange(args.start_time_point, args.end_time_point):
        print time_point
        save_atts.fileName = output_path + str(time_point).zfill(4) + ending
        visit.SetSaveWindowAttributes(save_atts)

        points_outer, points_inner = launch_points_inner_outer(field_nulls[time_point])

        if args.ion_forward_backward:
            points_outer = points_inner

        if args.current:
            points_inner = launch_points(field_nulls[time_point],
                                         current_thetas,
                                         radius=0.001)
            points_outer = launch_points(field_nulls[time_point],
                                         B_thetas,
                                         radius=0.005)

        visit.SetActivePlots(1)
        StreamlineAtts_flux_1.SetPointList(points_outer)
        visit.SetPlotOptions(StreamlineAtts_flux_1)

        visit.SetActivePlots(2)
        StreamlineAtts_flux_2.SetPointList(points_inner)
        visit.SetPlotOptions(StreamlineAtts_flux_2)
>>>>>>> canonical_flux_tubes/master

        visit.SetTimeSliderState(time_point)

        name = visit.SaveWindow()


def parse_args():
    parser = argparse.ArgumentParser(description="Generate time step plots of canonical flux tubes.")
    parser.add_argument('--database_prefix', help='path to visit database i.e. vtk files',
                        default='/home/jensv/rsx/jens_analysis/writing_to_vtk/output/')
    parser.add_argument('--database_postfix', help='path to visit database i.e. vtk files',
                        default='/Bdot_triple_probe_quantities*.vtk database')
    parser.add_argument('database_date', help='date of data run YYYY-MM-DD-mm-ss')
    parser.add_argument('--output_prefix', help='output_file_prefix',
                        default='electron_canonical_flux_tubes_')
    parser.add_argument('--current_min', help='minimum for current color map', default=0.0)
    parser.add_argument('--current_max', help='maximum for current color map', default=5.1e5)
<<<<<<< HEAD
    parser.add_argument('--start_time_point', help='time point of first output frame', type=int, default=0)
    parser.add_argument('--end_time_point', help='time point of last output frame', type=int, default=250)
    parser.add_argument('--field_nulls', help='path to file listing field_nulls (launching centers)',
                        default='/home/jensv/rsx/jens_analysis/centroid_fitting/output/2016-08-12/field_nulls.txt')
    parser.add_argument('--time_scale', help='time scale of time steps', default=0.068)
    parser.add_argument('--alpha_constant', help='value of spatially constant alpha',
                        type=int, default=8.1e5)
    parser.add_argument('--current_plane', help='plot temperature contours',
                        action='store_true', default=False)
    parser.add_argument('--temperature_tubes', help='plot temperature isosurfaces',
                        action='store_true', default=False)
    parser.add_argument('--density_tubes', help='plot density isosurfaces',
                        action='store_true', default=False)
    parser.add_argument('--electron', help='plot canonical electron flux tubes',
                        action='store_true', default=False)
    parser.add_argument('--ion', help='plot canonical ion flux tubes', action='store_true', default=False)
    parser.add_argument('--current',
                        help='plot thin current flux tube surrounded by electron / magnetic flux tube',
                        action='store_true', default=False)
    parser.add_argument('--interactive_session', action='store_true', default=False)
    parser.add_argument('--current_to_use', default='J_raw')
    parser.add_argument('--omega_to_use', default='Omega_i_raw_plus')
    parser.add_argument('--view', help='pre-configured_views: default, default_lower_angle, positive_z, negative_z, positive_x', default='default')
    parser.add_argument('--wait_for_manual_settings',
                        help='flag makes program wait for input before rendering time series.',
                        default=False, action='store_true')
    parser.add_argument('--double_stream', help='plot canonical streamline and magnetic of given variable', default=None)
    parser.add_argument('--x_offset', help='x offset of single streamline', default=0, type=int)
    parser.add_argument('--y_offset', help='y offset of single streamline', default=0, type=int)
    parser.add_argument('--file_number_offset', help='offset in file numbering', default=0, type=int)
    parser.add_argument('--turn_off_density_start', help='time step at which to start turning off density cloud.', type=int, default=None)
    parser.add_argument('--turn_off_density_end', help='time step at which to end turning off density cloud', type=int, default=None)
    parser.add_argument('--velocity', action='store_true', default=False)
=======
    parser.add_argument('--start_time_point', help='time point of first output frame', default=0)
    parser.add_argument('--end_time_point', help='time point of last output frame', default=250)
    parser.add_argument('--field_nulls', help='path to file listing field_nulls (launching centers)',
                        default='/home/jensv/rsx/jens_analysis/centroid_fitting/output/2016-08-12/field_nulls.txt')
    parser.add_argument('--time_scale', help='time scale of time steps', default=0.068)
    parser.add_argument('--alpha_constant', help='value of spatially constant alpha', type=int, default=8.1e5)
    task = parser.add_mutually_exclusive_group()
    task.add_argument('--electron', help='plot canonical electron flux tubes', action='store_true', default=False)
    task.add_argument('--ion', help='plot canonical ion flux tubes', action='store_true', default=False)
    task.add_argument('--ion_forward_backward',
                      help='plot canonical ion flux tubes',
                      action='store_true', default=False)
    task.add_argument('--current',
                      help='plot thin current flux tube surrounded by electron / magnetic flux tube',
                      action='store_true', default=False)
    parser.add_argument('--interactive_session', action='store_true', default=False)
    parser.add_argument('--current_to_use', default='J_raw')
>>>>>>> canonical_flux_tubes/master
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    main()
