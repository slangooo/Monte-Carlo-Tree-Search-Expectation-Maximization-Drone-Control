% % function R = getCoverageRadius2(sigma,SNRThresh,totalNoisePowerDb, h, userBandwidth, totalBandwidth, avgLossNLOS,avgLossLOS, carrFreq,a,b) 
% % Pt = sigma;  
% % receivedPower =10.^(SNRThresh/10).* 10.^(totalNoisePowerDb./10);
% % %   receivedPower = userBandwidth./totalBandwidth .* Pt .* db2mag(-PathLoss);
% %   PathLoss = -10.*log10(receivedPower./userBandwidth.*totalBandwidth./Pt);
% % %   PathLoss = 20.*log(4*pi.*carrFreq.*dist./ physconst('LightSpeed')) +...
% % %       PLOS.*avgLossLOS + (1 - PLOS).*avgLossNLOS;
% %   
% % %   r =sqrt((10.^((PathLoss - PLOS.*avgLossLOS + (1 - PLOS).*avgLossNLOS)./20).*physconst('LightSpeed')./(4*pi.*carrFreq)).^2-h.^2)
% %   
% % init_guess = 0.5;
% % % PathLoss = 20.*log10(4*pi*carrFreq.*sqrt(r.^2+h.^2)./299792458) +...
% % %     1./(1+a.*exp(-b.*(180./pi.*atan(h./r)-a))).*(avgLossLOS - avgLossNLOS) + avgLossNLOS;
% % 
% % PLOS = init_guess*ones(1,size(h,2));
% % r = getrFromPlos(h,a,b,PLOS);
% % PLOS = getPlosFromR(PathLoss,avgLossLOS,avgLossNLOS,r,carrFreq,h);
% % % r =sqrt((10.^((PathLoss - PLOS.*avgLossLOS - (1 - PLOS).*avgLossNLOS)./20)...
% % %          .*299792458./(4*pi.*carrFreq)).^2-h.^2);
% % error = abs(init_guess - r);
% % iterations =0;
% % while any(error > 0.000001)%&& iterations < 1000
% %     iterations = iterations +1;
% %     rPrev = r;
% %     r = real(getrFromPlos(h,a,b,PLOS));
% % %     r = abs(r);
% % %     r(isnan(r))=0;
% %     PLOS = real(getPlosFromR(PathLoss,avgLossLOS,avgLossNLOS,r,carrFreq,h));
% % %     PLOS(PLOS<0) =0;
% % %     PLOS(PLOS>1) =1;
% %     error = abs(rPrev - r);
% % end
% % r(imag(r)~=0)=0;
% % R = r;

function R = getCoverageRadius2(sigma,SNRThresh,totalNoisePowerDb, h, userBandwidth, totalBandwidth, avgLossNLOS,avgLossLOS, carrFreq,a,b) 
Pt = sigma;  
receivedPower =10.^(SNRThresh/10).* 10.^(totalNoisePowerDb./10);
%   receivedPower = userBandwidth./totalBandwidth .* Pt .* db2mag(-PathLoss);
  PathLoss = -10.*log10(receivedPower./userBandwidth.*totalBandwidth./Pt);
%   PathLoss = 20.*log(4*pi.*carrFreq.*dist./ physconst('LightSpeed')) +...
%       PLOS.*avgLossLOS + (1 - PLOS).*avgLossNLOS;
r = 1:5:4000;
R = zeros(1,length(h));
for i=1:length(h)
    H= h(i);
    P = 20.*log10(4*pi*carrFreq.*sqrt(r.^2+H.^2)./299792458) +...
        1./(1+a.*exp(-b.*(180./pi.*atan(H./r)-a))).*(avgLossLOS - avgLossNLOS) + avgLossNLOS;
    [v,jj]=min(abs(P-PathLoss(i)));
    R(i)=r(jj);
end
r= R;