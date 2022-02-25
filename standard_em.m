x = [mvnrnd([5 10],reshape([2 2],1,2).*eye(2),30); mvnrnd([0 -15],reshape([8 1],1,2).*eye(2),50);mvnrnd([-10 -10],reshape([4 8],1,2).*eye(2),50)];
%x=[1 5; 1 2; 3 4; 3 5];
m =3; % number of components
n=size(x,1); %number of test points
mu =[0 5; -10 -15; -10 -10];
sigma = var(x).*eye(2);
sigma = repmat(sigma,1,1,m);
%sigma = reshape(1:5:10,1,1,m).*sigma;

Pj = ones(1,m)/m; %initializing Pj
T=zeros(n,m); %=Pjt
while 1
%     input('')
    hold off;
    scatter(x(:,1),x(:,2));
    xlim([-40 20])
    ylim([-40 20])
    hold on;
    for j=1:m
        a=sigma(1,1,j); % horizontal radius
        b=sigma(2,2,j); % vertical radius
        x0=mu(j,1); % x0,y0 ellipse centre coordinates
        y0=mu(j,2);
        t=-pi:0.01:pi;
        xp=x0+a*cos(t);
        yp=y0+b*sin(t);
        plot(xp,yp)
    end
    for j = 1:m
        T(:,j)=Pj(j).*mvnpdf(x,mu(j,:),sigma(:,:,j));
    end
    T= T./sum(T,2);
    Probs = T;
    entro = -log2(Probs).*Probs;
    entro(isnan(entro))=0;
    sum(sum(entro'))
    Pj =sum(T,1);
    for j=1:m
        mu(j,:)=[T(:,j)'*x(:,1) T(:,j)'*x(:,2)]./Pj(j);
        sigma(:,:,j)= [T(:,j)'*(x(:,1)-mu(j,1)).^2 0;0 T(:,j)'*(x(:,2)-mu(j,2)).^2]./Pj(j);
    end
    Pj = Pj./n;
end
% % % SINR = getSINR(m,x,mu,sigma);
% % % plotNetwork(x,sigma,mu,[min(x(:,1)) max(x(:,1))],[min(x(:,2)) max(x(:,2))]);