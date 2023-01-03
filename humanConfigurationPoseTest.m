load('HumanConfigurationInputs.mat')

% leftHip = HumanConfigurationInputs(:,1:3);
% leftKnee = HumanConfigurationInputs(:,4:6);
% leftAnkle = HumanConfigurationInputs(:,7:9);
% rightHip = HumanConfigurationInputs(:,10:12);
% rightKnee = HumanConfigurationInputs(:,13:15);
% rightAnkle = HumanConfigurationInputs(:,16:18);

xSamp = HumanConfigurationInputs;


centers = zeros(18,length(xSamp));
Fcenters = zeros(18,length(xSamp));
M = 1;
centers(:,M) = xSamp(1,:)';
Fcenters(:,M) = xSamp(2,:)';
separation = 10;
indexes = zeros(1,length(xSamp));
indexes(1) = 1;

Findexes = zeros(1,length(xSamp));
Findexes(1) = 1;
for mm = 1:length(xSamp)
    check = 0;
    for jj = 1:M
        if poseMetric(centers(:,jj),xSamp(mm,:)') > separation
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
beta = 1;
K = zeros(M,M);
Kf = zeros(M,M);
type = 'matern32';

    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = poseKernel(type,centers(jj,:),centers(pp,:),beta);
        end
    end

    
sigma = svd(K);

plot(log(sigma))
ylabel('log(\sigma)')


% regressionCoefs =ySamp(indexes)'*pinv(K);
% 
% 
% 
% x1Range = minGrid(1):2:maxGrid(1);
% x2Range = minGrid(2):2:maxGrid(2);
% 
% [X1,X2] = meshgrid(x1Range,x2Range);
% 
% 
% 
% kernEstimate = zeros(size(X1));
% 
% for ii = 1:size(kernEstimate,1)
%     for jj = 1:size(kernEstimate,2)
%         kernVector = zeros(length(centers),1);
%         for kk = 1:length(centers)
%             kernVector(kk) = kernel(type,centers(kk,:),[X1(ii,jj),X2(ii,jj)],beta);
%         end
%         kernEstimate(ii,jj) = regressionCoefs*kernVector;
%     end
% end
% 
% surf(X1,X2,kernEstimate)
% grid on
% view(3)
% % plot3(xSamp(:,1), xSamp(:,2),ySamp, 'r.', 'LineWidth', 1, 'MarkerSize', 5);
