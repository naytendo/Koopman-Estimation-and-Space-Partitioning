function [xSamp,ySamp] = samplesOnCircle(numSamp, radius, center)

% f1 = @(x)(x(1)^2*exp(-x(2)^2)-x(2)^3*exp(x(1)));
G = @(x)(-sin(2*x(2))+1/9*x(1)^2-1/2)+2;
% Create the circle
angles = linspace(0, 2*pi, 1080); % 1080 is the total number of points

x_1Center = center(1);
x_2Center = center(2);
x_1 = radius * cos(angles) + x_1Center; 
x_2 = radius * sin(angles) + x_2Center;

% Plot circle.
base = zeros(1,length(x_1));
plot3(x_1, x_2,base, 'k-', 'LineWidth', 2);
% % Plot center.
% hold on;
% plot(x_1Center, x_2Center, 'k+', 'LineWidth', 2, 'MarkerSize', 16);
grid on;
axis equal;
xlabel('$X_1$','fontsize',20,'interpreter','latex')
ylabel('$X_2$','fontsize',20,'interpreter','latex')

text(1.4,0,0,'$Q$','fontsize',16,'interpreter','latex')

% zlabel('G', 'FontSize', 16, 'interpreter','latex');
% % Now get random locations along the circle.
s1 = numSamp; % Number of random points to get.
Indexes = 1:floor(length(x_1)/s1):length(x_1);
x_1Samp = x_1(Indexes);
x_2Samp = x_2(Indexes);
xSamp = [x_1Samp;x_2Samp];
ySamp = zeros(1,length(xSamp));
for ii = 1:length(xSamp)
    ySamp(ii) = G(stateEvolution(xSamp(:,ii)));
end


x1Range = min(x_1)-1:0.1:max(x_1)+1;
x2Range = min(x_2)-1:0.1:max(x_2)+1;

[X1,X2] = meshgrid(x1Range,x2Range);

hold on


[X1,X2] = meshgrid(x1Range,x2Range);

X1f = zeros(size(X1));
X2f = zeros(size(X2));

for ii = 1:size(X1f,1)
    for jj = 1:size(X1f,2)
        XNew = stateEvolution([X1(ii,jj);X2(ii,jj)]);
        X1f(ii,jj) = XNew(1);
        X2f(ii,jj) = XNew(2);
    end
end

surf(X1,X2,(-sin(2*X2)+1/9*X1.^2-1/2)+2,'EdgeColor','none','FaceAlpha',0.5,'FaceColor','yellow');
plot3(xSamp(1,:), xSamp(2,:),ySamp, 'r--', 'LineWidth', 2, 'MarkerSize', 5);
text(-1,-1,3,'$G$','fontsize',16,'interpreter','latex')
xSamp = xSamp';
ySamp = ySamp';


surf(X1,X2,(-sin(2*X2f)+1/9*X1f.^2-1/2)+2,'EdgeColor','none','FaceAlpha',0.5,'FaceColor','green');
% plot3(xSamp(1,:), xSamp(2,:),ySamp, 'r--', 'LineWidth', 2, 'MarkerSize', 5);
text(-1,1.2,3,'$U_f G$','fontsize',16,'interpreter','latex')


end