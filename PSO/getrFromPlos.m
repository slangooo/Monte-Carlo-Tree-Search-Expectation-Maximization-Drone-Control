function r = getrFromPlos(h,a,b,Plos)
r = h./(tan((((log((1./Plos-1)./a))./-b)+a)./(180/pi)));