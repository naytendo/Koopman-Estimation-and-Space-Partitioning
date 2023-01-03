figure()
hold on

for depth = 3

radius =1;
center = [0;0];
numSamp = [2^1 2^2 2^3 2^4 2^5];

Error =  zeros(1,length(numSamp));
fillDist = zeros(1,length(numSamp));

numChild = 2;


a = 5;
b = [1,0]';


maxGrid = [1.25,1.25];
minGrid = [-1.25,-1.25];


for ii = 1:length(numSamp)

% [xSamp,ySamp] = samplesOnCircle(numSamp,radius,center);

% [xSamp,ySamp] = samplesOverGrid(numSamp,7);

[xSamp,ySamp] = samplesOnCircle(numSamp(ii),radius,center);

% [coefs_n,numSamps_n] = averageOverGridNodes(xSamp,ySamp,depth,numChild,maxGrid,minGrid);

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
type = 'matern32';
    for pp = 1:M
        for jj = 1:M
            K(pp,jj) = kernel(type,centers(jj,:),centers(pp,:),beta);
        end
    end

regressionCoefs = (ySamp'*pinv(K))';
[xSamp2,ySamp2] = samplesOnCircle(2^6,radius,center);

approx = zeros(1,length(xSamp2));
for jj = 1:length(xSamp2)
    approx(jj) = getKernelEstimate(xSamp2(jj,:),regressionCoefs,centers,type,beta);
end


error = (ySamp2 - approx');

Error(ii) = abs(max(error));
fillDist(ii) = 2*pi/numSamp(ii);


end


order = 3/2;

t = 2*order-1-0.25;
m = floor(t);
plot(log10(numSamp),log10(Error),'-o')
plot(log10(numSamp),(t-m)*log10(fillDist),'--')



end

ylabel('log$\| U_f G - U^M_f G\|_Y$', 'interpreter','latex');
xlabel('log(Number of Samples)')

text(1,0,'$C h_{\Xi_n,Q}^{(t-m)}$','interpreter','latex')

