numSamp = 2^8;


x = zeros(length(numSamp),2);

x(1,1) = 0;
x(1,2) = pi/3;
delta = 0.1;
for ii =2:numSamp
    x(ii,:) = pendulumDynamics(x(ii-1,:),delta);
end



xSamp = x';
ySamp = zeros(length(xSamp),1);
G = @(x)(-sin(x(2))+1/9*x(1)^2-1/2)+2;

for ii = 1:length(xSamp)
    ySamp(ii) = G(pendulumDynamics(xSamp(:,ii),delta));
end


depth = 4;
numChild = 2;
maxGrid = [2,6];
minGrid = [-2,0];



approx = zeros(1,length(xSamp));


centers = zeros(2,length(xSamp));
M = 1;
centers(:,M) = xSamp(:,1);
separation = 0.25;
indexes = zeros(1,length(xSamp));
indexes(1) = 1;
for mm = 1:length(xSamp)
    check = 0;
    for jj = 1:M
        if norm((centers(:,jj)-xSamp(:,mm))') > separation
            check = check +1;
        end
    end
    if check == M
        centers(:,M+1) = xSamp(:,mm);
        M = M+1;
        indexes(M) = mm;
    end
end
% centers = xSamp(:,1:3:end)';
centers = centers(:,1:M)';
indexes = indexes(1:M);
beta = 0.2;
K = zeros(M,M);
Kf = zeros(M,M);
type = 'wendland32';

    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
            Kf(pp,jj) = kernel(type,centers(jj,:),pendulumDynamics(centers(pp,:),delta),beta);
        end
    end

regressionCoefs = (K\K\Kf\ySamp(indexes))';
res = 0.1;
x1Range = minGrid(1):res:maxGrid(1);
x2Range = minGrid(2):res:maxGrid(2);

[X1,X2] = meshgrid(x1Range,x2Range);
X1f = zeros(size(X1));
X2f = zeros(size(X2));

for ii = 1:size(X1f,1)
    for jj = 1:size(X1f,2)
        XNew = pendulumDynamics([X1(ii,jj);X2(ii,jj)],delta);
        X1f(ii,jj) = XNew(1);
        X2f(ii,jj) = XNew(2);
    end
end
kernEstimate = zeros(size(X1));

for ii = 1:size(kernEstimate,1)
    for jj = 1:size(kernEstimate,2)
        kernVector = zeros(length(centers),1);
        for kk = 1:length(centers)
            kernVector(kk) = kernel(type,centers(kk,:),[X1f(ii,jj),X2f(ii,jj)],beta);
        end
        kernEstimate(ii,jj) = regressionCoefs*kernVector;
    end
end


% surf(X1,X2,(-sin(2*X2)+1/9*X1.^2-1/2)+2,'EdgeColor','none','FaceAlpha',0.5,'FaceColor','yellow');

% plot3(xSamp(1,:), xSamp(2,:),ySamp, 'r-', 'LineWidth', 2, 'MarkerSize', 5);
hold on
plot3(xSamp(1,:), xSamp(2,:),zeros(length(xSamp)), 'k', 'LineWidth', 2, 'MarkerSize', 5);
grid on

% text(-1,-1,3,'$G$','fontsize',16,'interpreter','latex')


% surf(X1,X2,(-sin(X2f)+1/9*X1f.^2-1/2)+2,'EdgeColor','none','FaceAlpha',0.5,'FaceColor','green');
plot3(xSamp(1,:), xSamp(2,:),ySamp, 'r-', 'LineWidth', 2, 'MarkerSize', 5);
% text(-1.75,5.75,3,'$U_f G$','fontsize',16,'interpreter','latex')
text(-1.75,4.5,0,'$Q$','fontsize',16,'interpreter','latex')
surf(X1,X2,kernEstimate,'FaceAlpha',0.9,'FaceColor','blue')



% plot3(xSamp(:,1), xSamp(:,2),ySamp, 'r.', 'LineWidth', 1, 'MarkerSize', 5);
xlabel('$x^1$','fontsize',20,'interpreter','latex')
ylabel('$x^2$','fontsize',20,'interpreter','latex')
zlabel('$U_f^M G$','fontsize',20,'interpreter','latex')
view(-67,30)

set(gca,'FontSize',16);