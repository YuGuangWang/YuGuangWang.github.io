function [w,y] = pntset(quadr, deg)
% [w,y] = pntset(ip, deg)
% returns the weights w and nodes y (size:[N 3]) of quadrature rule ip
%
% Inputs:
% ip -- name of quadrature rule
% deg -- i) degree which the quadrature rule is exact at or associated with
%
% Outputs:
% w -- weights of the quadrature rule
% y -- nodes of the quadrature rule


switch quadr
    case 'SD'
        [w,y] = SD(deg);
    case 'GL'
        [w,y] = GL(deg);
    case 'SP'
        N = 2*deg^2;
        [w,y] = SP(N);
end
w = w';
y = y';