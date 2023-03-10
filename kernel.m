function kern = kernel(type,x,y,beta)
% Calculate value of kernel given the type of kernel k(.,.),two elements
% and parameters.
% Input:
%       type    -   string, type of kernel 
%       x,y     -   column vectors
%       beta    -   column vectors, parameters of kernel (indications to be added)
% Output:
%       KK      -   k(x,y)
r = norm(x-y);
% Examine the inputs (to be added)
switch type
    case 'exp'
        Dxy = r^2;
        kern = exp(-beta * Dxy);
    case 'wendland31'
        if 1-r <= 0
            kern = 0;
        else
            kern = (1-r)^4*(4*r+1);
        end
     case 'wendland32'
        if 1-r <= 0
            kern = 0;
        else
            kern = (1-r)^6*(35*r^2+18*r+3);
        end
     case 'wendland33'
        if 1-r <= 0
            kern = 0;
        else
            kern = (1-r)^8*(32*r^3 + 25*r^2+8*r+1);
        end
    case 'imq'
        kern = 1/(5^2+r^2)^beta;
    case 'matern32'
        kern = (1 + sqrt(3)*r/beta)*exp(-sqrt(3)*r/beta);
    case 'matern52'
        kern = (1 + sqrt(5)*r/beta + 5*r^2/beta^2)*exp(-sqrt(5)*r/beta);
    case 'matern'
        kern = 1/(8*beta^3)*(1+beta*r)*exp(-beta*r);
    otherwise 
   
        error('Kernel type is not supported!');
end