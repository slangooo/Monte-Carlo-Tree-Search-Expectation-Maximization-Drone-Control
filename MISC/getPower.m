function P=getPower(R, pathLossExponent,referPower,Noise,h, SNR_Thresh)
d=sqrt(R.^2+h^2);
P=db2mag(SNR_Thresh+10*log10(d.^pathLossExponent)-referPower+Noise);
end