

radius =1;
center = [0;0];
numSamp = 2^6;



% [xSamp,ySamp] = samplesOnCircle(numSamp,radius,center);

% [xSamp,ySamp] = samplesOverGrid(numSamp,7);
[xSamp,ySamp] = samplesOnCircle(numSamp,radius,center);

depth = 4;
numChild = 2;
maxGrid = [1.5,1.5];
minGrid = [-1.5,-1.5];

% [coefs_n,numSamps_n] = averageOverGridNodes(xSamp,ySamp,depth,numChild,maxGrid,minGrid);
% 
% c = coefs_n';
% n = length(c);
% centers = zeros(n,2);
% 
% 
% % for index = 1:length(c)
% %     plot2DCell(index-1,depth,numChild,maxGrid,minGrid,coefs_n(index))
% % end
% 
% for index = 1:n 
%     centers(index,:) = get2DKernelCenter(index-1,depth,numChild,maxGrid,minGrid);
%     
% %     plot3(centers(index,1),centers(index,2),0,'kx', 'LineWidth', 1, 'MarkerSize', 5);
%     
% end
% 
% beta = 0.13;
% K = zeros(n,n);
% type = 'matern32';
% for ii = 1:n
%     for jj = 1:n
%         K(ii,jj) = kernel(type,centers(jj,:),centers(ii,:),beta);
%     end
% end
% 
% regressionCoefs = (c'*pinv(K))';

approx = zeros(1,length(xSamp));
% c = coefs_n';
% n = length(c);
% centers = zeros(n,2);
% for index = 1:n 
%     centers(index,:) = get2DKernelCenter(index-1,depth,numChild,maxGrid,minGrid);
%     
% %     plot3(centers(index,1),centers(index,2),0,'kx', 'LineWidth', 1, 'MarkerSize', 5);
%     
% end

centers = xSamp;
M = length(xSamp);
beta = 1;
K = zeros(M,M);
Kf = zeros(M,M);
type = 'matern32';

    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
            Kf(pp,jj) = kernel(type,centers(jj,:),stateEvolution(centers(pp,:)')',beta);
        end
    end

regressionCoefs = (ySamp'*inv(K)*inv(K)*Kf)';
res = 0.1;
x1Range = minGrid(1):res:maxGrid(1);
x2Range = minGrid(2):res:maxGrid(2);

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
kernEstimate = zeros(size(X1));

for ii = 1:size(kernEstimate,1)
    for jj = 1:size(kernEstimate,2)
        kernVector = zeros(length(centers),1);
        for kk = 1:length(centers)
            kernVector(kk) = kernel(type,centers(kk,:),[X1f(ii,jj),X2f(ii,jj)],beta);
        end
        kernEstimate(ii,jj) = regressionCoefs'*kernVector;
    end
end

surf(X1,X2,kernEstimate,'FaceAlpha',0.7,'FaceColor','blue')
grid on
hold on
view(135,30)
text(-1.2,1.5,0.25,'$U_f^M G$','fontsize',16,'interpreter','latex')
% plot3(xSamp(:,1), xSamp(:,2),ySamp, 'r.', 'LineWidth', 1, 'MarkerSize', 5);

