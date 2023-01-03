function kernEstimate = getKernelEstimate(xSamp,coefs,centers,type,beta)

% determines what the estimate is is for a particular sample so we can
% compare it to the actual function
kernVector = zeros(size(centers,1),1);
for kk = 1:size(centers,1)
       kernVector(kk) = kernel(type,centers(kk,:),xSamp,beta);
end

kernEstimate = coefs'*kernVector;



end