import numpy as np
from kinetic_model import active_Cas
from kinetic_model import dead_Cas
from kinetic_model.change_concentration import  change_concentration
from scipy import linalg
from scipy.optimize import curve_fit

################################################################################
# mimic the NucleaSeq from Finkelsteinlab
#
#
################################################################################


def calc_cleavage_rate(guide,target,epsilon,forward_rates, Cas):
    '''
    Effective cleavage rate as measured in NucleaSeq experiments (Jones et al.)
    -----
    1. Using saturating concentrations of Cas9-sgRNA,
    2. determine the fraction of cut DNA as a function of time
    3. Fit exponential decay (straight line in logarithmic space) ---> decay constant (slope) defines the effective cleavage rate

    :param guide: Guide sequence
    :param target: Target sequence
    :param epsilon: Energetic parameters
    :param forward_rates: Forward rate parameters
    :param Cas: instance of CRISPRclass.CRISPR()
    :return:
    '''
    # --- NucleaSeq is performed at saturating Cas-sgRNA concentrations ---
    epsilon_saturate, forward_rates_saturate = change_concentration(epsilon, forward_rates,
                                                                           new_concentration=10**9.,
                                                                           ref_concentration=1.0)

    # --- prepare Master Equation ---
    guide_length = Cas.guide_length
    mismatch_positions = active_Cas.get_mismatch_positions(guide, target, Cas)
    rate_matrix = active_Cas.get_master_equation(epsilon_saturate,forward_rates_saturate, mismatch_positions, guide_length)


    # --- Determine fraction of uncut DNA at specified timepoints (dictated by experiment) ---
    times = np.array([0.0, 12.0, 60.0, 180.0, 600.0, 1800.0, 6000.0, 18000.0, 60000.0])
    prob_uncleaved = np.zeros(len(times))
    everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))
    for i,time in enumerate(times):
        probs = dead_Cas.get_Probability(rate_matrix, everything_unbound,time)
        prob_uncleaved[i] = np.sum(probs)


    # --- Fit an exponential decay (straight line in log-space) ---
    # --- take the logarithm of the probabilities (zero values should not be considered)
    end = len(times)
    for i in range(len(times)):
        if prob_uncleaved[i] == 0:
            end = i
            break
    times = times[0:end]
    prob_uncleaved = prob_uncleaved[0:end]
    prob_log = np.log(prob_uncleaved)



    # fit the log of the probability to a linear function,
    # yielding the cleavage rate



    k, error = curve_fit(line, times, prob_log)
    return k[0]


def line(x, k):
    return -k * x