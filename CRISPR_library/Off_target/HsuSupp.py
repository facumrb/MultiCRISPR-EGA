import warnings,os


compTable = { "a":"t", "A":"T", "t" :"a", "T":"A", "c":"g", "C":"G", "g":"c", "G":"C", "N":"N", "n":"n", 
        "Y":"R", "R" : "Y", "M" : "K", "K" : "M", "W":"W", "S":"S",
        "H":"D", "B":"V", "V":"B", "D":"H", "y":"r", "r":"y","m":"k",
        "k":"m","w":"w","s":"s","h":"d","b":"v","d":"h","v":"b","y":"r","r":"y" }

def complRna(seq):
    " complement the sequence and translate to RNA "
    newseq = []
    for nucl in seq.upper():
        newseq.append( compTable[nucl].replace("T", "U") )
    return "".join(newseq)

def outMats():
    """
    write normalized matrices to out 
    """
    prettyMat("all", "out/hsuAll.tsv")
    prettyMat("row", "out/hsuRow.tsv")
    prettyMat("col", "out/hsuCol.tsv")
    prettyMat("none", "out/hsuNone.tsv")

def prettyMat(strat, outFname):
    """
    wrote normalized Hsu matrix to outFname
    """
    ofh = open(outFname, "w")
    row = ["pos"]
    row.extend(range(2,21))
    row = [str(x) for x in row]
    ofh.write( "\t".join(row)+"\n")

    normMat, normAvg  = parseHsuMat("./hsu2013/fig2cData.txt", strat=strat)
    for nucl, freqs in normMat.iteritems():
        row = [":".join(nucl)]
        row.extend(freqs)
        row = [str(x) for x in row]
        ofh.write( "\t".join(row)+"\n")

    row = ["avg"]
    row.extend(normAvg)
    row = [str(x) for x in row]
    ofh.write( "\t".join(row)+ "\n")

hsuMat = None # dict with (fromNucl, toNucl) -> list of 19 scores
avgFreqs = None # list of 19 scores
hsuStrat = None # loaded matrix has strat

# see hsuMat.py
nuclFreqs= {('A', 'A'): 0.4819440141370613, ('G', 'G'): 0.6571297543038187, ('U', 'T'): 0.4663759334374003, ('U', 'C'): 0.18755561795635842, ('C', 'T'): 0.3917125484856841, ('G', 'A'): 0.472948896301865, ('G', 'T'): 1.0, ('A', 'G'): 0.2796160896968995, ('U', 'G'): 0.787929779020387, ('C', 'C'): 0.0, ('A', 'C'): 0.6804018833297995, ('C', 'A'): 0.5931243444910334}
posFreqs = [0.294369386, 0.29164666, 0.190210984, 0.306896763, 0.167251773, 0.219909422, 0.169797251, 0.406475982, 0.628680509, 0.588183598, 0.907111342, 0.522909141, 1.256594863, 0.693851359, 0.552888666, 1.158572718, 1.662766602, 1.01548686, 1.428913939]


def loadHsuMat(strat):
    if strat in ["avgs", "raw"]:
        return
    global hsuMat
    global avgFreqs
    #global hsuStrat
    matFname = "./hsu2013/fig2cData.txt"
    hsuMat, avgFreqs = parseHsuMat(matFname, strat)
    #hsuStrat = strat

def parseHsuMat(fname, strat="col"):
    """ return the hsu 2013 matrix as a dict rnaNucl -> dnaNucl -> list of scores and a list of 19 averages, one per position
    >> parseHsuMat("./hsu2013/fig2cData.txt")
    """
    hsuMat = {}
    minMat = 99999.0
    maxMat = 0.0
    for line in open(fname):
        if line.startswith("nucl"):
            continue
        fs = line.rstrip("\n").split()
        # the values are in the order 19-1 3'-5' in the file, but our sequences are always 1-19, 5'-3'
        freqs = list(reversed([float(x) for x in fs[1:]]))
        if line.startswith("avg"):
            avgs = freqs
            continue
        nuclComb = fs[0]
        rnaNucl, dnaNucl = nuclComb.split(":")
        hsuMat[ (rnaNucl, dnaNucl)] = freqs
        minMat = min(min(freqs), minMat)
        maxMat = max(max(freqs), maxMat)

    minCols = []
    maxCols = []
    for i in range(0, 19):
        colVals = []
        for freqs in hsuMat.values():
            colVals.append(freqs[i])
        minCols.append(min(colVals))
        maxCols.append(max(colVals))

    pCount = 0.0001 # must use pseudo counts
    # normalize
    # "Each frequency was normalized to range from 0 to 1, such that f = (f-fmin) / (fmax-fmin)"
    normMat = {}
    for key, freqs in hsuMat.items():
        if strat=="all":
            normVals = [( f - minMat) / (maxMat - minMat) for f in freqs]
        elif strat=="row":
            minFreq, maxFreq = min(freqs), max(freqs)
            normVals = [( f - minFreq) / (maxFreq - minFreq) for f in freqs]
        elif strat=="col":
            normVals = [( f - minCols[i]) / (maxCols[i] - minCols[i]) for i, f in enumerate(freqs)]
        elif strat.startswith("none") or strat=="onlyAvgs":
            normVals = freqs
        elif strat.startswith("limit"):
            normVals = [min(f, 1.0) for f in freqs]
        else:
            assert(False)

        if not strat.startswith("none"):
            normVals = [pCount+n for n in normVals]
        normMat[key] = normVals
    hsuMat = normMat

    if strat in ["all", "row", "onlyAvgs", "limit"]:
        normAvgs = [(a-min(avgs)) / (max(avgs)-min(avgs)) for a in avgs]
    else:
        normAvgs = avgs

    if not strat.startswith("none"):
        normAvgs = [n+pCount for n in normAvgs]
    assert(min(normAvgs)!=0.0)
    avgs = normAvgs

    #print "loaded hsu mat", strat
    #print hsuMat
    #print avgs
    return hsuMat, avgs

def calcHsuSuppScore(guideSeq, otSeq, strat, dfh=None):
    """ calculate the score described on page 17 of the Hsu et al 2013 supplement PDF
    Will ignore position 0 of both the ot and the guide as the Hsu score is only
    defined for positions 1-20
    

    # mismatch in 5' part -> 

    >> calcHsuSuppScore("GAGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAG")
    0.4509132855355929
    >> calcHsuSuppScore("GTGTCCGAGCAGAAGAAGAA","GAGTCCGAGCAGAAGAAGAA")
    0.007929899123452079
    >> calcHsuSuppScore("GAGTCCGAGCAGAAGAAGAA","GAGTCAGAACAGAAGAACAA")
    3.69683024017458e-08
    >> calcHsuSuppScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
    1.0
    >> calcHsuSuppScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGT")
    0.4907323091938597
    >> calcHsuSuppScore("GGGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGC")
    0.00880386936389443
    >> calcHsuSuppScore("GAGGGGGGGGGGGGGGGGGG","GGGGGGGGGGGGGGGGGGGG")
    0.0016461094044462768
    """

    if len(guideSeq)==23:
        guideSeq = guideSeq[:20]
        otSeq = otSeq[:20]
    elif len(guideSeq) < 20:
        warnings.warn('Running HsuSupp. The length of the target sequence is less than 20. This sequence will be skipped because it may not be valid for the scoring algorithm.')
        return None
    else:
        guideSeq = guideSeq[:20]
        otSeq = otSeq[:20]
        warnings.warn('Running HsuSupp. The length of the target sequence is not 23. Maybe your target sequence is not in the format of 20nt guide+3nt PAM, which may cause accuracy problems.')

    guideSeq = guideSeq[1:20]
    otSeq = otSeq[1:20]
    #print guideSeq, len(guideSeq), otSeq, len(otSeq)
    assert(len(guideSeq)==19)
    assert(len(otSeq)==19)# Hsu ignores pos 0

    global hsuMat
    global avgFreqs
    global hsuStrat
    if hsuMat is None or hsuStrat!=strat:
        matFname = os.path.join(os.path.dirname(__file__),"hsu2013/fig2cData.txt")
        hsuMat, avgFreqs = parseHsuMat(matFname, strat)
        hsuStrat = strat

    rnaSeq = complRna(guideSeq)
    if dfh:
        dfh.write("guideDna=%s, guideRna=%s, otSeq=%s\n" % (guideSeq, rnaSeq, otSeq))
    # "Predicted cutting frequencies for genome-wide targets were calculated by
    # multiplying, in series: fest = f(1) * g(N1,N1') * f(2) * g(N2,N2') * ... * h
    # with values f(i) and g(Ni, Ni')
    # at position i corresponding, respectively, to the aggregate
    # position- and base-mismatch cutting frequencies for positions and pairings indicated in Fig. 2c"
    mismatchPosList = []
    score = 1.0
    for i in range(0, 19):
        rnaNucl, dnaNucl = rnaSeq[i].upper(), otSeq[i].upper()
        # "In case of a match, both were set equal to 1."
        if (rnaNucl, dnaNucl) in [('C', 'G'), ('U', 'A'), ('A', 'T'), ('G', 'C')]:
            f = 1.0
            g = 1.0
        else:
            f = avgFreqs[i]
            g = hsuMat[(rnaNucl, dnaNucl)][i]
            mismatchPosList.append(i)
        if strat.endswith("_Sum"):
            score += f*g
        elif strat.endswith("_allSum"):
            score += f+g
        else:
            score *= f * g
        if dfh:
            dfh.write("pos: %d, RNA-nucl: %s, DNA-nucl: %s -> f=%f, g=%f, score=%f\n" % ( i, rnaNucl, dnaNucl, f, g, score))
        #score *= g

    # "The value h meanwhile re-weighted the estimated
    # frequency by the minimum pairwise distance between consecutive mismatches in the target
    # sequence. This distance value, in base-pairs, was divided by 18 to give a maximum value of 1 (in
    # cases where fewer than 2 mismatches existed, or where mismatches occurred on opposite ends of
    # the 19 bp target-window"
    if len(mismatchPosList)<2:
        h = 1.0
    else:
        dists = []
        for left, right in zip(mismatchPosList[:-1], mismatchPosList[1:]):
            dists.append(right-left)
        minDist = min(dists)
        h = minDist / 18.0
    score *= h

    if dfh:
        dfh.write("h=%f\n" % h)
        dfh.write("score=%f\n\n" % score)
    return score