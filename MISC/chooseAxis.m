function Pts = chooseAxis (x,M,gran)
m=M-1;
center = mean(x);
if (m==0)
    Pts=1;
end
X= x-center;
dtheta=2*pi/m;
y=zeros(2,m);
maxSet=y;
for i=0:m-1
    maxSet(:,i+1)=[cos(i*dtheta);sin(i*dtheta)]/sqrt(cos(i*dtheta)^2+sin(i*dtheta)^2);
end
MagMax=X*maxSet;
MagMax=sum(max((X*maxSet)'));
for j=1:gran-1
    for i= 0:m-1
        y(:,i+1)=[cos(i*dtheta+j*(dtheta/gran));sin(i*dtheta+j*(dtheta/gran))]/sqrt(cos(i*dtheta+j*(dtheta/gran))^2+ sin(i*dtheta+j*(dtheta/gran))^2);
    end
    mult=X*y;
    value = sum(max(mult'));
    if value>MagMax
        MagMax=value;
        maxSet=y;
    end
end
Pts=zeros(M,1);
[V,I]=min(sqrt(X(:,1).^2+X(:,2).^2));
Pts(1)=I;
for i=2:m+1
    [V,Pts1]=max(max(X*maxSet));
    [V,Pts2]=max(X*maxSet);
    Pts(i)=Pts2(Pts1);
    for j=1: size(x,1)/(m+1)
        [V,TT]=max(X*maxSet(:,Pts1));
        X(TT,:)=nan;
    end
    %X(Pts(i),:)=nan;
    maxSet(:,Pts1)=nan;
end