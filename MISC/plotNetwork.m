function []=plotNetwork(x,sigma,mu,xLimit,yLimit)
m=size(mu,1);
scatter(x(:,1),x(:,2));
xlim(xLimit)
ylim(yLimit)
hold on;
for j=1:m
    a=sigma(j); % horizontal radius
    b=sigma(j); % vertical radius
    x0=mu(j,1); % x0,y0 ellipse centre coordinates
    y0=mu(j,2);
    t=-pi:0.01:pi;
    xp=x0+a*cos(t);
    yp=y0+b*sin(t);
    plot(xp,yp)
end
end