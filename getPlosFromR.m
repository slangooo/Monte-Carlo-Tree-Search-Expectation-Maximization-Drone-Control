function PLOS = getPlosFromR(PathLoss,avgLossLOS,avgLossNLOS,r,carrierFreq,h)
PLOS = (PathLoss-avgLossNLOS-20.*log10(4.*pi.*carrierFreq.*sqrt(r.^2+h.^2)./299792458))./(avgLossLOS - avgLossNLOS);