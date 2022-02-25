function PLOS = getProbLOS (h,r,a,b)
PLOS = 1./(1+a.*exp(-b.*(180/pi.*atan(h./r)-a)));