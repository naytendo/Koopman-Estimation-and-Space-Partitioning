for depth = 3

radius =1;
center = [0;0];
numSamp = 2^8;

fillDist = [0.2 0.15 0.1 0.06 0.03 0.01 0.003];
condNum31 =  zeros(1,length(fillDist));
condNum32 =  zeros(1,length(fillDist));
condNum33  =  zeros(1,length(fillDist));

for ii = 1:length(fillDist)

numChild = 2;



maxGrid = [2,6];
minGrid = [-2,0];




% [xSamp,ySamp] = samplesOnCircle(numSamp,radius,center);

% [xSamp,ySamp] = samplesOverGrid(numSamp,7);
x = zeros(length(numSamp),2);

x(1,1) = 0;
x(1,2) = pi/3;
delta = 0.001;
for jj =2:numSamp
    x(jj,:) = pendulumDynamics(x(jj-1,:),delta);
end



xSamp = x';
ySamp = zeros(length(xSamp),1);
G = @(x)(-sin(x(2))+1/9*x(1)^2-1/2)+2;

for jj = 1:length(xSamp)
    ySamp(jj) = G(pendulumDynamics(xSamp(:,jj),delta));
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
beta = 1;

K31 = zeros(M,M);
K32 = zeros(M,M);
K33 = zeros(M,M);
type31 = 'wendland31';
type32 = 'wendland32';
type33 = 'wendland33';
    for pp = 1:M
        for jj = 1:M
            K31(pp,jj) = kernel(type31,centers(jj,:),centers(pp,:),beta);
            K32(pp,jj) = kernel(type32,centers(jj,:),centers(pp,:),beta);
            K33(pp,jj) = kernel(type33,centers(jj,:),centers(pp,:),beta);
        end
    end
condNum31(ii) = cond(K31);
condNum32(ii) = cond(K32);
condNum33(ii) = cond(K33);
% e = eig(K31);
% lambdaMin(ii) = min(e);
% lambdaMax(ii) = max(e);

end
figure()
plot(-log(0.5*fillDist),log(condNum31),'linewidth',2)
hold on
plot(-log(0.5*fillDist),log(condNum32),'--','linewidth',2)
plot(-log(0.5*fillDist),log(condNum33),'-.','linewidth',2)
ylabel('log(cond(K))');
xlabel('-log($q_x$)','interpreter','latex')
grid on
legend('$C^2$','$C^4$','$C^6$','interpreter','latex','location','southeast');
% 
% figure(2)
% plot(-log(0.5*fillDist),log(1./lambdaMin),'b','linewidth',2)
% hold on
% plot(-log(0.5*fillDist),log(lambdaMax),'r','linewidth',2)
% xlabel('-log(q_x)')
% grid on

end