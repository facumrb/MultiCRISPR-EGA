#########################################################
# code to extract model parameters from SA fit
#########################################################
import numpy as np
import pandas as pd


def change_concentration(epsilon_1nM, forward_rates_1nM, new_concentration, ref_concentration=1.0):
    '''
    translate parameters to new concentration
    Here it is assumed you started with paramters at 1nM
    (can change using the 'ref_concentration')
    :param epsilon_1nM:
    :param forward_rates_1nM:
    :param new_concentration:
    :param ref_concentration:
    :return:
    '''

    # -- adjust epsilon PAM ---
    epsilon_new = np.array(epsilon_1nM).copy()
    epsilon_new[0] -= np.log(new_concentration/ref_concentration)

    # --- adjust binding rate from solution ---
    forward_rates_new = np.array(forward_rates_1nM).copy()
    forward_rates_new[0] *= new_concentration/ref_concentration
    return epsilon_new, forward_rates_new




