function PL = getPathLoss (carrFreq, dist, PLOS, avgLossLOS, avgLossNLOS)
    PL = 20.*log10(4*pi.*carrFreq.*dist./299792458) + PLOS.*avgLossLOS + (1 - PLOS).*avgLossNLOS;
    