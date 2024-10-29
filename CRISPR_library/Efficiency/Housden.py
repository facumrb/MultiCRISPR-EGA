from math import log10
import warnings
# Housden matrix (see function below):
# an array of 4x20=80 floats. The first 20 are for A, the next 20 for T, then C, then G
# imported from original file received from the authors: matrix_final.txt
factors = [
    0.4979, 0.7959, 0.7553, 0.6569, 0.9481, 0.7147, 0.437, 0.6212, 0.9077, 1.0, 0.1957, 0.7959, 0.6212, 0.8912, 1.0, 0.5485, 0.9942, 0.5485, 0.455, 1.0, \
    0.6699, 0.5485, 0.275, 0.5972, 0.6212, 0.7555, 1.0, 0.5131, 0.8608, 0.7553, 0.6569, 0.3417, 1.0, 0.016, 0.9146, 0.7555, 0.2906, 0.4979, 0.5485, 0.5131, 
    0.4979, 0.6869, 0.8528, 0.7643, 0.5325, 0.3417, 0.3417, 0.7643, 0.6434, 0.0092, 0.9331, 0.5325, 0.7272, 0.9708, 0.2905, 0.7272, 0.2957, 0.7918, 0.6434, 0.5062, \
    0.7918, 0.4461, 0.4851, 0.4461, 0.3417, 0.6869, 0.2417, 0.5485, 0.0947, 0.9256, 0.5325, 0.8308, 0.1255, 0.7918, 0.2544, 0.4461, 0.4979, 0.6212, 0.7918, 0.4461
]

def calcHousden(seq):
    """
    Calc housden score and round to one decimal point.
    Based on java file Crispr.java received from the authors.

    Housden: Housden et al, PMID 26350902, http://www.flyrnai.org/evaluateCrispr/
    FlyRNAi.org—the database of the Drosophila RNAi screening center and transgenic RNAi project: 2017 update

    >>> calcHousden(["ATCTGACCTCCCGGCTAATT"])
    [6.9]
    """
    seq = seq.upper()

    if len(seq) != 20:
        warnings.warn('Running Housden. The length of the guide seq is not 20. This sequence will be skipped because it may not be valid for the scoring algorithm.')
        return None
    
    if "N" in seq: # cannot do Ns
        warnings.warn('Running Housden. Guide seq contains ambiguous nucleotide "N". This sequence will be skipped as it may not be valid for the scoring algorithm.')
        return None
        

    assert(len(seq)==20)
    nuclToIndex = {"A":0,"T":1,"C":2,"G":3}

    score = 1.0
    for i in range(0, 20):
        nuclIndex = nuclToIndex[seq[i]]
        idx = (nuclIndex*20)+i
        score *= factors[idx]
    score = -1*log10(score)
    score = float("%0.1f" % score) # round to one decimal point
    return score