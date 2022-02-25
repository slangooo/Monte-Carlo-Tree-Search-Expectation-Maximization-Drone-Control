function R = getCoverageRadius(sigma,SNRThresh,totalNoisePowerDb, h, userBandwidth, totalBandwidth, avgLossNLOS,avgLossLOS, carrFreq,a,b) 
Pt = sigma;  
receivedPower =10.^(SNRThresh/10).* 10.^(totalNoisePowerDb./10);
%   receivedPower = userBandwidth./totalBandwidth .* Pt .* db2mag(-PathLoss);
  PathLoss = -10.*log10(receivedPower./userBandwidth.*totalBandwidth./Pt);
%   PathLoss = 20.*log(4*pi.*carrFreq.*dist./ physconst('LightSpeed')) +...
%       PLOS.*avgLossLOS + (1 - PLOS).*avgLossNLOS;
  
%   r =sqrt((10.^((PathLoss - PLOS.*avgLossLOS + (1 - PLOS).*avgLossNLOS)./20).*physconst('LightSpeed')./(4*pi.*carrFreq)).^2-h.^2)
  
init_guess = 20;
PLOS = getProbLOS (h,init_guess,a,b);

r =sqrt((10.^((PathLoss - PLOS.*avgLossLOS - (1 - PLOS).*avgLossNLOS)./20)...
         .*299792458./(4*pi.*carrFreq)).^2-h.^2);
error = abs(init_guess - r);
iterations =0;
while any(error > 0.000001)%&& iterations < 1000
    iterations = iterations +1;
    if iterations > 10000
        
    end
    rPrev = r;
    PLOS = getProbLOS (h,r,a,b);
    r =real(sqrt((10.^((PathLoss - PLOS.*avgLossLOS - (1 - PLOS).*avgLossNLOS)./20)...
         .*299792458./(4*pi.*carrFreq)).^2-h.^2));
    error = abs(rPrev - r);
end
r(imag(r)~=0)=0;
R = r;