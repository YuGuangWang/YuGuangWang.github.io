function F = NeedletPsi(j, X, w, h, Y)
% F = NeedletPsi(j, X, w, h, Y)
% Evaluate needlet psi_{jk} at points in Y on S^d, d >= 2
% -- Input arguments -- 
% j = order of needlet: j = 0, 1, 2, ....
%     needlet includes Gegenbauer polynomials of degree deg: 2^{j-2} < deg < 2^j
% X = d+1 by Nj array of nodes of needlet quadrature rule
% w = 1 by Nj array of weights ofr needlet quadrature rule: w > 0
%     Needlet quadrature rule exact for polynomials of degee <= 2^{j+1}
% h = needlet filter with support (1/2, 2): h(t)^2 + h(2*t)^2 = 1
% Y = d+1 by N array of points where needlet is to be evaluated
% -- Output arguments --
% F = Nj by N array of values of the Nj needlets of order j
%     evaluated at the N points on S^d specified as columns of Y

% Needlet filter h has support [1/2, 2]: continuity ==> h(1/2) = h(2) = 0
% Thus order j needlet only includes Gengenbauer polynomials of degree deg:
% where 2^{j-2} < deg < 2^j, so
% Order j = 0  <==>  deg =  0
% Order j = 1  <==>  deg =  1
% Order j = 2  <==>  deg =  2, 3
% Order j = 3  <==>  deg =  3, 4, 5, 6, 7
% Order j = 4  <==>  deg =  5, ..., 15
% Order j = 5  <==>  deg =  9, ..., 31
% Order j = 6  <==>  deg = 17, ..., 63
% etc

% Dimension of spherical harmonics of degree L on S^d
Zdl = @(d, L) (2*L+d-1).*gamma(L+d-1)./(gamma(d).*gamma(L+1));

[d1, Nj] = size(X);
[d1a, N] = size(Y);
if d1a ~= d1
    fprintf('Error in NeedletPsi: Input points must be on S^d\n');
    fprintf('size(X) = %d, %d, size(Y) = %d, %d\n', size(X), size(Y));
    psi = NaN;
    return
end
% Dimension of sphere
d = d1 - 1;
% Parameter of Gegenbauer polynomial for S^d
% Jacobi: alpha = beta = lambda - 1/2 = (d-2)/2
lam = (d-1)/2;

% Force needlet quadrature weights to be a column vector
w = w(:);
wsqrt = sqrt(w);

% Vector of ones of size N
eN = ones(1,N);

% Order j = 0 <==> degree 0 polynomial
% Normalized Gegenbauer polynomial of degree 0 is C(z) = 1
% needlet filter has h(1) = 1
if j == 0
    F = wsqrt(:,eN).*ones(Nj,N);
    return;
end

% Order j = 1 term <==> degree 1 polynomial
% Nornalized Gegenbauer polynomial of degree 1 is C(z) = z
if j == 1
    coeff = Zdl(d,1)*wsqrt;
    %F0 = zeros(Nj, N);
    %for i = 1:N
    %    yi = Y(:,i);
    %    z = X'*yi;
    %    F0(:,i) = coeff.*z;
    %end;
    % Array syntax
    Z = X'*Y;
    F = coeff(:,eN).*Z;
    %err1 = norm(F-F0, inf)
    return;
end

% Order j > 1
   
% Degrees and coefficients in three-term recurrence
% for Gegenbauer polynomials
Deg = [1:2^j-1];
An = 2*(Deg+lam-1)./Deg;
Cn = -(Deg+2*lam-2)./Deg;
%
Z = X'*Y;
P0 = ones(size(Z));
P1 = (2*lam)*Z;

% Gegenbauer polynomials not included in needlet of order j
for n = 2:2^(j-2)
    P = An(n)*Z.*P1 + Cn(n)*P0;
    P0 = P1;
    P1 = P;
end;

% Gegenbauer polynomials included in needlet of order j
Deg = [2^(j-2)+1:2^j-1];
hDeg = h(Deg/2^(j-1));
% Yuguang's simplificattion of Z(d,ell) / C_ell^lam(1)
aDeg = (2/(d-1))*Deg + 1;
coeff = hDeg.*Zdl(d,Deg);
F = zeros(Nj, N);
for i = 1:numel(Deg)
    
    n = Deg(i);
    P = An(n)*Z.*P1 + Cn(n)*P0;
    P0 = P1;
    P1 = P;
    
    % Normalization with check of Yuguang's simplification
    Pn1 = gamma(n+2*lam)/gamma(n+1)/gamma(2*lam);
    const0 = coeff(i)/Pn1;
    const1 = hDeg(i)*aDeg(i);
    chk = norm(const0-const1, inf);
    if chk > 10*n*eps
        keyboard
    end
    F = F + const1*P;  
    
end;
    
F = wsqrt(:,eN).*F;


