figure()
hold on
type = ["exp", "matern32", "wendland32"];
lineType = ["-", "-.", "--"];
for tt = 1:length(type)

    radius =1;
    center = [0;0];
    numSamp = 2^10;

    fillDist =  [0.5 0.2 0.15 0.1 0.06 0.05 0.03 0.01];
    Error =  zeros(1,length(fillDist));
    numCenters =  zeros(1,length(fillDist));



    numChild = 2;



    maxGrid = [2,6];
    minGrid = [-2,0];

    for ii = 1:length(fillDist)


        x = zeros(length(numSamp),2);

        x(1,1) = 0;
        x(1,2) = pi/3;
        delta = 0.0001;
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
        indexes = indexes(1:M);
        beta = 1;
        K = zeros(M,M);
        

            for pp = 1:M
              for jj = 1:M
                    K(pp,jj) = kernel(type(tt),centers(jj,:),centers(pp,:),beta);
          
              end
            end
        
        Error(ii) =cond(K);
    end

plot(-log10(0.5*fillDist),log10(Error),lineType(tt),'linewidth',2)

end


legend('Exponential','Sobolev-Matern','Wendland','location','northwest')
set(gcf,'color','w');
set(gca,'fontsize',20)
% xlim([1 2.6])
grid on