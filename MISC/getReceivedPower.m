function Pr = getReceivedPower(PathLoss, Pt, userBandwidth, totalBandwidth)
    Pr = userBandwidth./totalBandwidth .* Pt .* 10.^(-PathLoss/10);