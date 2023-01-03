figure()
hold on



for depth = 3

radius =1;
center = [0;0];
numSamp = 2^8;

x = zeros(length(numSamp),2);

x(1,1) = 0;
x(1,2) = pi/3;
delta = 0.01;
for jj =2:numSamp
    x(jj,:) = pendulumDynamics(x(jj-1,:),delta);
end



xSamp = x';
ySamp = zeros(length(xSamp),1);

TFCoefs = [4 2 2 3 2]';
TFCenters = xSamp(:,[25  75 115 160  200])';

typeTF = 'matern52';
betaTF = 1;
for ii = 1:length(xSamp)
    xSampNew =  pendulumDynamics(xSamp(:,ii),delta);
    ySamp(ii) = xSampNew(1) + xSampNew(2)+ getKernelEstimate(xSampNew,TFCoefs,TFCenters,typeTF,betaTF);
end

fillDist = [1.25 1 0.8 0.67 0.5];
Error =  zeros(1,length(fillDist));
numCenters =  zeros(1,length(fillDist));



numChild = 2;



maxGrid = [2,6];
minGrid = [-2,0];


for ii = 1:length(fillDist)




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
% centers = xSamp(:,1:3:end)';
centers = centers(:,1:M)';
indexes = indexes(1:M);
beta = 1;
K = zeros(M,M);
Kf = zeros(M,M);
type = 'matern32';

    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
            Kf(pp,jj) = kernel(type,centers(jj,:),pendulumDynamics(centers(pp,:),delta),beta);
        end
    end

regressionCoefs = K\K\Kf\ySamp(indexes);



approx = zeros(1,length(xSamp));
for jj = 1:length(xSamp)
    approx(jj) = getKernelEstimate(xSamp(:,jj)',regressionCoefs,centers,type,beta);
end


error = (ySamp - approx');

Error(ii) = abs(max(error));
numCenters(ii) = size(centers,1);

end


order = 3/2;

t = 2*order-1-0.2;
m = floor(t);

plot(-log10(fillDist),log10(Error),'-o','linewidth',4)
plot(-log10(fillDist),0.5+2.5*(t-m)*log10(fillDist),'--','linewidth',4,'color',[0.900 0.3850 0.1280])



end

ylabel('log$\| U_f G - U^M_f G\|_Y$', 'interpreter','latex','fontsize',34);
xlabel('-log($h_{\Xi_M,Q}$)','interpreter','latex','fontsize',34)
set(gca,'fontsize',34,'color',[0.6350 0.0780 0.3])
set(gcf,'color',[1 1 1])
xlim([0 0.32])
text(.18,0.18,'$C + (t-m)log(h_{\Xi_n,Q})$','interpreter','latex','fontsize',30,'color',[1 1 1])


