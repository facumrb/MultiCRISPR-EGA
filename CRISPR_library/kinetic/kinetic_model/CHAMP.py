import numpy as np
from scipy.optimize import curve_fit
from kinetic_model import dead_Cas
from kinetic_model import active_Cas
################################################################################
# mimic the CHAMP from Finkelsteinlab
#
################################################################################


def calc_ABA(guide, target, epsilon, forward_rates, Cas):
    '''
    Apparant Binding Affinity (ABA) as measured in CHAMP experiments (Jones et al.)
    ----
    1. Determine fraction of DNA bound by dCas9-sgRNA after 10 minutes of exposure as a function of concentration
    2. Fit Hill-curve to determine the dissociation constant
    3. ABA is defined as the natural logarithm of the dissociation constant

    :param guide: Guide sequence
    :param target: Target sequence
    :param epsilon: Energetic parameters
    :param forward_rates: Forward rate parameters
    :param Cas: instance of CRISPRclass.CRISPR()
    :return:
    '''

    concentrations = np.array([0.1,0.3,1.,3.,10.,30.,100.,300.])
    Kd, _,_ = binding_curve(guide,target,epsilon,forward_rates, concentrations,Cas, time=10*60)
    return np.log(Kd)


def calc_delta_ABA(guide, target, epsilon, forward_rates, Cas):
    on_target_ABA = calc_ABA(guide, guide, epsilon, forward_rates, Cas)
    ABA = calc_ABA(guide, target, epsilon, forward_rates, Cas)
    return ABA-on_target_ABA



def binding_curve(guide,target,epsilon,forward_rates, concentrations, Cas,time):

    # --- prepare master equation ---
    guide_length = Cas.guide_length
    mismatch_positions = active_Cas.get_mismatch_positions(guide, target, Cas)
    rate_matrix = dead_Cas.get_master_equation(epsilon,forward_rates, mismatch_positions, guide_length)


    everything_unbound = np.array([1.0] + [0.0] * (guide_length + 1))
    Pbound = []
    for c in concentrations:
        new_rate_matrix = rate_matrix.copy()
        new_rate_matrix[0][0] *= c
        new_rate_matrix[1][0] *= c
        Probability = dead_Cas.get_Probability(new_rate_matrix, everything_unbound, time)
        Pbound.append(np.sum(Probability[1:]))
    Pbound = np.array(Pbound)
    concentrations = np.array(concentrations)
    Kd, _ = curve_fit(Hill_eq, concentrations,Pbound)
    return Kd[0], Pbound, concentrations





def Hill_eq(C, Kd):
    return (1.0+Kd/C)**(-1)