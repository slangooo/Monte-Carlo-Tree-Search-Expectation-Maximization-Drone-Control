function SNR = getSNR (receivedPower, totalNoisePowerDb)
    SNR = receivedPower / 10.^(totalNoisePowerDb./10);