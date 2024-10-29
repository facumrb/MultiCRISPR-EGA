from subprocess import  CalledProcessError, call
import warnings
import platform,tempfile, bisect
import os
from os.path import  join,dirname,isdir,isfile


global binDir

# by default bindir is relative to the location of this library
binDir = join(dirname(__file__), "bin")

def setBinDir(path):
    global binDir
    binDir = path

def getBinPath(name, isDir=False):
    """
    get the full pathname of a platform-specific binary, in the bin/ directory relative to this directory
    """
    currPlatform = platform.system()
    binPath = join(binDir, currPlatform, name+'.exe')
    if isDir and not isdir(binPath):
        raise Exception("Could not find directory %s" % binPath)
    if not isDir and not isfile(binPath):
        raise Exception("Could not find file %s" % binPath)
    return binPath

def seqsToChariSvml(seqs):
    """ partially copied from generateSVMFile.FASTA.py in the Chari et al source code
    >>> seqsToChariSvml(["CTTCTTCAAGGTAACTGCAGA", "CTTCTTCAAGGTAACTGGGGG"])
    '0 13:1 22:1 32:1 43:1 52:1 62:1 73:1 84:1 94:1 101:1 111:1 122:1 134:1 144:1 153:1 162:1 171:1 183:1 194:1 201:1 214:1\\n0 13:1 22:1 32:1 43:1 52:1 62:1 73:1 84:1 94:1 101:1 111:1 122:1 134:1 144:1 153:1 162:1 171:1 181:1 191:1 201:1 211:1'
    """
    vecs = []
    for seq in seqs:
        assert(len(seq)==21)
        vec = []
        # end index
        for pos in range(0, 21):
            for nuclIdx, char in enumerate("GTCA"):
                val = int(seq[pos]==char)
                if val!=0:
                    vec.append( ("%d%d" % (pos+1, nuclIdx+1), val) )
        vecs.append( vec )

    lines = []
    for vec in vecs:
        vec = ["%s:%s" % (x,y) for x,y in vec]
        lines.append("0 "+" ".join(vec))
    return "\n".join(lines)

chariRanges = None

def calcChariScore(seq):
    """ return dict with chari 2015 scores, returns two lists (rawScores, rankPercent)
    Chari: Chari et al, PMID 26167643 http://crispr.med.harvard.edu/sgRNAScorer
    Unraveling CRISPR-Cas9 genome engineering parameters via a library-on-library approach. Nature Methods. 2015;12(9):823-826

    input seqs have lengths 21bp: 20 bp guide + 1bp first from PAM
    >>> calcChariScores(["CTTCTTCAAGGTAACTGCAGA", "CTTCTTCAAGGTAACTGGGGG"])
    ([0.54947621, 0.58604487], [80, 81])
    >>> calcChariScores(["CTTCTTCAAGGNAACTGCAGA"])
    ([0.9025848], [88])
    """

    if len(seq) != 21:
        warnings.warn(f'Running Chari. The length of the sequence is not 21bp. This sequence will be skipped because it may not be valid for the scoring algorithm.')
        return None
        
    # this is a rewritten version of scoreMySites.py in the Chari2015 suppl files
    chariDir = join(binDir, "src", "sgRNA.Scorer.1.0")
    modelFname = join(chariDir, '293T.HiSeq.SP.Nuclease.100.SVM.Model.txt')
    dataIn = seqsToChariSvml([seq])

    tempFname = tempfile.mktemp()

    with open(tempFname, "w") as tempFh:
        tempFh.write(dataIn + "\n")
    
    outName = tempfile.mktemp()

    svmlPath = getBinPath("svm_classify")
    cmd = [svmlPath, "-v", "0", tempFname, modelFname, outName]
    
    try:
        proc = call(cmd)
    except CalledProcessError:
        raise Exception("Could not run command '%s'" % (" ".join(cmd)))

    dataOut = open(outName).read().strip()
    
    os.remove(outName)
    os.remove(tempFname)
    
    try:
        score = float(dataOut)
    except ValueError:
        score = None  # Handle any issues with the score formatting

    return score
