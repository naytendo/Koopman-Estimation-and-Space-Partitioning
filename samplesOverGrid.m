function [xSamp,ySamp] = samplesOverGrid(numSamp, maxGrid)

% f1 = @(x)(x(1)^2*exp(-x(2)^2)-x(2)^3*exp(x(1)));
f2 = @(x)(-sin(x(2))+1/9*x(1)^2-1/2);
% Create the circle

xSamp = maxGrid*rand(2,numSamp);

ySamp = zeros(1,length(xSamp));
for ii = 1:length(xSamp)
    ySamp(ii) = f2(xSamp(:,ii));
end


x1Range = min(xSamp(1,:))-1:0.5:max(xSamp(1,:))+1;
x2Range = min(xSamp(2,:))-1:0.5:max(xSamp(2,:))+1;

[X1,X2] = meshgrid(x1Range,x2Range);

hold on
surf(X1,X2,(-sin(X2)+1/9*X1.^2-1/2),'EdgeColor','none','FaceAlpha',0.7);
plot3(xSamp(1,:), xSamp(2,:),ySamp, 'rx', 'LineWidth', 1, 'MarkerSize', 5);

xSamp = xSamp';
ySamp = ySamp';

end