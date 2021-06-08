function phi = Wendland_nr(r, k, delta, n)
% phi = Wendland_nr(r, k, delta, n)
% Evaluate the normalised equal area RBF 
% (2^(k-1)*gamma(k)*gamma(ell+2*k+1))/(gamma(ell+1)*gamma(2*k)) phi_k(r/delta_ellka)
% where r is the Euclidean distance for 0 <= r <= 2, and
% delta_ellka = (ell+2k+1)*gamma(k+1/2)/(2*sqrt(a)*gamma(k+1)).
% Index k determines which RBF
% Each RBK is normalized so phi(0) = 1 (affects k = 3, 5)
% k = -1 : phi(r) = max(1-r, 0)
% k = 0 : phi(r) = max(1-r, 0)^2 in C^0
% k = 1 : phi(r) = max(1-r, 0)^4 * (4*r+1) in C^2
% k = 2 : phi(r) = max(1-r, 0))^6 * (35*r^2+18*r+3)/3 in C^4
% k = 3 : phi(r) = max(1-r, 0))^8 * (32*r^3+25*r^2+8*r+1) in C^6
% k = 4 : phi(r) = max(1-r, 0))^10 * (429*r^4+450*r^3+210*r^2+50*r+5) in C^8
% k = 5 : phi(r) = max(1-r, 0))^12 * (2048*r^5+2697*r^4+1644*r^3+566*r^2+108*r+9) in C^10
% 
% Inputs:
% r     : variable
% delta : Scaling factor delta > 0 (Default delta = 1)
% n     : Dimension of points x, y: r = ||x - y||
%       : Default value n = 2 for S^2 in R^3
%
% Outputs:
% phi   : function value of normalised RBF, an array of the same size to input r
%
% See the Maple worksheet rbf_wendland
%
% Ref1: Holger Wendland "Piecewise polynomial, positive definite and
% compactly supported radial function of minimal degree",
% Advances in Computational Mathemaitcs 4 (1995) 389-396.
%
% Ref2: A. Chernih, I. H. Sloan, R. S. Womersley. "Wendland functions with
% increasing smoothness converge to a Gaussian",
% Advances in Computational Mathematics 40 (2014): 185-200.

% Default to S^2 in R^3
if nargin < 4
    n = 2;
end;
if nargin < 3
    delta = 1;
end;
if delta ~= 1
    r = r/delta;
end;

% normalise the RBF for k larger than 1
if k>0||k==0
    a = 1;
    ell = floor(k + n/2)+1;
%     coe_w = 2^(k-1)*gamma(k)*gamma(ell+2*k+1)/(gamma(ell+1)*gamma(2*k));
    delta_ellka = (ell+2*k+1)*gamma(k+1/2)/(2*sqrt(a)*gamma(k+1));
    r = r/delta_ellka;
end

% Positive par of 1 - r
rp = max(1-r, 0);

% Select which RBF
switch k
    
    case {-1}
        
        % Hat function
        phi = max(1 - r, 0);
        
    case {0}
        
        % Wendland v = 2, k = 0 function in C^0 and H_s(S^2) for s = 3/2
        phi = rp.^2;
               
        
    case {1}
        
        % Wendland v = 3, k = 1 function in C^2 and H_s(S^2) for s = 5/2
        phi = rp.^4 .* (4*r+1);
        
        
    case {2}
        
        % Wendland v = 4, k = 2 function in C^4 and H_s(S^2) for s = 7/2
        phi = rp.^6 .* ((35*r+18).*r+3)/3;
        
                
    case {3}
        
        % Wendlandv = 5, k = 3 function in C^6 and H_s(S^2) for s = 9/2
        phi = rp.^8 .* (((32*r+25).*r+8).*r+1);
              
        
    case {4}
        
        % Wendland v = 6, k = 4 function in C^8
        phi = rp.^10 .* ((((429*r+450).*r+210).*r+50).*r+5)/5;
        
        
    case {5}
        
        % Wendland v = 7, k = 5 function in C^10
        phi = rp.^12 .* (((((2048*r+2697).*r+1644).*r+566).*r+108).*r+9)/9;
        
        
    case {6}
        
        % Wendland v = 8, k = 6 function in C^12
        phi = rp.^14 .* ((((((46189*r+73206).*r+54915).*r+24500).*r+6755).*r+1078).*r+77)/77;
        
        
    otherwise
        
        fprintf('RBF warning: Unknown case k = %d\n', k);
        phi = [];
        return;
        
end
% if k > 0
%     phi = phi*coe_w;
% end