function m = polyfit(x,y,lambda,n)
%Polynomial fitting
%   input: x, regressors
%          y, output values
%          lambda, regularization parameter
%          n, degree of polynomial
% output: estimate, same structure as LinRegressRegul

% Note that we are considering a multi-dimensional case

% new regressor matrix
x2 = polyinput(x, n);

m = LinRegressRegul(x2,y,lambda);
m.n = n;
if n == 1
    if lambda == 0
        m.label = 'LR';
    else
        m.label = sprintf('LRR (.3f)', lambda);
    end
else
    if lambda == 0
        m.label = sprintf('Poly %d', n);
    else
        m.label = sprintf('Poly %d (%.3f)', n, lambda);
    end
end
end