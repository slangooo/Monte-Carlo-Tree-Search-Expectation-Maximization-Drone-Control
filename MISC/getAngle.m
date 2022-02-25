function angle = getAngle (u,v)
CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
angle = real(acosd(CosTheta)); %degrees
end