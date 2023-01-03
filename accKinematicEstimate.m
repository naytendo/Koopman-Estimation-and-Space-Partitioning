%% 
load('HumanTreadmillDataACC.mat')

hip = markerData(:,1:3);
knee = markerData(:,4:6);
ankle = markerData(:,7:9);
theta1 = zeros(1,length(hip));

theta2 = zeros(1,length(hip));
first = 1;
last = 10000;
for ii = 1:length(theta1)
    [theta1(ii),theta2(ii)] = getJointAngles(hip(ii,:),knee(ii,:),ankle(ii,:));
end
% figure()
% plot(xSamp(:,1),xSamp(:,2));
xSamp = [theta1(first:last)',theta2(first:last)']; %input data
ySamp = Output_Q(first:last,2); %output data
depth = 5;
numChild = 2;

maxGrid = [65,70];
minGrid = [-35,-5];
% 
% [coefs_n,numSamps_n] = averageOverGridNodes(xSamp,ySamp,depth,numChild,maxGrid,minGrid);
% 
% 
% 
% c = coefs_n';
% n = length(c);
% centers = zeros(n,2);
% figure()
% hold on
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
% beta = 30;
% K = zeros(n,n);
% type = 'matern32';
% for ii = 1:n
%     for jj = 1:n
%         K(ii,jj) = kernel(type,centers(ii,:),centers(jj,:),beta);
%     end
% end
% 
% regressionCoefs = (ySamp'*pinv(K))';

centers = zeros(2,length(xSamp));
Fcenters = zeros(2,length(xSamp));
M = 1;
centers(:,M) = xSamp(1,:)';
Fcenters(:,M) = xSamp(2,:)';
separation = 5;
indexes = zeros(1,length(xSamp));
indexes(1) = 1;

Findexes = zeros(1,length(xSamp));
Findexes(1) = 1;
for mm = 1:length(xSamp)
    check = 0;
    for jj = 1:M
        if norm((centers(:,jj)-xSamp(mm,1))') > separation
            check = check +1;
        end
    end
    if check == M
        centers(:,M+1) = xSamp(mm,:)';
        M = M+1;
        indexes(M) = mm;
        Findexes(M) = mm +1;
        if mm ~= length(xSamp)
            Fcenters(:,M+1) = xSamp(mm+1,:)';
        else
            Fcenters(:,M+1) = xSamp(1,:)';
        end
    end
end
% centers = xSamp(:,1:3:end)';
centers = centers(:,1:M)';
Fcenters = Fcenters(:,1:M)';
indexes = indexes(1:M);
Findexes = Findexes(1:M);
beta = 0.5;
K = zeros(M,M);

type = 'matern';

    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
        end
    end
%%
regressionCoefs =ySamp(indexes)'*pinv(K);

%%

x1Range = minGrid(1):2:maxGrid(1);
x2Range = minGrid(2):2:maxGrid(2);

[X1,X2] = meshgrid(x1Range,x2Range);



kernEstimate = zeros(size(X1));

for ii = 1:size(kernEstimate,1)
    for jj = 1:size(kernEstimate,2)
        kernVector = zeros(length(centers),1);
        for kk = 1:length(centers)
            kernVector(kk) = kernel(type,centers(kk,:),[X1(ii,jj),X2(ii,jj)],beta);
        end
        kernEstimate(ii,jj) = regressionCoefs*kernVector;
    end
end

surf(X1,X2,kernEstimate,'Edgecolor','none','FaceAlpha',0.8)
hold on
% spacing = 5;  % play around so it fits the size of your data set
% for i = 1 : spacing : length(X1(1,:))
%     plot3(X1(:,i), X2(:,i), kernEstimate(:,i),'-','color',[0.4 0.4 0.4]);
% end

% for i = 1 : spacing : length(X1(:,1))
%     plot3(X1(i,:), X2(i,:), kernEstimate(i,:),'-','color',[0.4 0.4 0.4]);
% end

plot3(xSamp(:,1),xSamp(:,2),zeros(length(xSamp),1),'k','linewidth',2)

grid on
view(3)
plot3(theta1(first:last), theta2(first:last),ySamp, 'r.', 'LineWidth', 1, 'MarkerSize', 5);