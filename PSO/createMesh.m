function [subX, subY, X, Y] = createMesh (xBoundary, yBoundary, smallSteps, largeSteps)

xSegs = linspace(0, xBoundary(3), largeSteps);
ySegs = linspace(0, yBoundary(2), largeSteps);

xSubSegs = linspace(0, xBoundary(3), smallSteps);
ySubSegs = linspace(0, yBoundary(2), smallSteps);

[X,Y] = meshgrid(xSegs, ySegs);
[subX,subY] = meshgrid(xSubSegs, ySubSegs);
% 
% figure(1)
% plot(X,Y,'k')
% hold on
% plot(Y,X,'k')
% % hold off
% plot(subX,subY,':y')
% % hold on
% plot(subY,subX,':y')
end
