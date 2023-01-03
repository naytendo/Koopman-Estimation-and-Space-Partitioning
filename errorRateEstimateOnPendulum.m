figure()
hold on

for depth = 3

radius =1;
center = [0;0];
numSamp = 2^11;

fillDist = [10 1 0.5 0.25 0.12 0.08 0.04 0.01 0.005 0.002];
Error =  zeros(1,length(fillDist));
OneOverNumCenters =  zeros(1,length(fillDist));
conditionNums=  zeros(1,length(fillDist));

maxGrid = [2,6];
minGrid = [-2,0];


for ii = 1:length(fillDist)

% [xSamp,ySamp] = samplesOnCircle(numSamp,radius,center);

% [xSamp,ySamp] = samplesOverGrid(numSamp,7);
x = zeros(length(numSamp),2);

x(1,1) = 0.5;
x(1,2) = pi/2;
dt = 0.001;
delta = 1;
for jj =2:numSamp
    x(jj,:) = pendulumDynamics(x(jj-1,:),dt);
end



xSamp = x';
ySamp = zeros(length(xSamp),1);
G = @(x)(-sin(x(2))+1/9*x(1)^2-1/2)+2;

for jj = 1:length(xSamp)
    ySamp(jj) = G(pendulumDynamics(xSamp(:,jj),dt));
end


centers = zeros(2,length(xSamp));
M = 1;
centers(:,M) = xSamp(:,1);
separation = fillDist(ii);
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
centers = centers(:,1:M)';
indexes = indexes(1:M);
beta = 1;
K = zeros(M,M);

type = 'matern52';
    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
        end
    end

regressionCoefs = K\ySamp(indexes);



approx = zeros(1,length(xSamp));
for jj = 1:length(xSamp)
    approx(jj) = getKernelEstimate(xSamp(:,jj)',regressionCoefs,centers,type,beta);
end


error = (ySamp - approx');
conditionNums(ii) = cond(K);
Error(ii) = abs(max(error));
OneOverNumCenters(ii) = 1./size(centers,1);

end




loglog(-log10(OneOverNumCenters),log10(Error),'-o')
loglog(-log10(OneOverNumCenters),2.5+2*log10(fillDist),'--')
loglog(-log10(OneOverNumCenters),log10(conditionNums),'-.')


end

ylabel('log$\| U_f G - U^M_f G\|_Y$', 'interpreter','latex');
xlabel('-log(1/N)','interpreter','latex')

legend('error','theoretical error rate','condition number')
