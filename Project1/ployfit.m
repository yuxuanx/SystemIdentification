function m = ployfit(x,y,lambda,n)
%Polynomial fitting
%   input: x, regressors
%          y, output values
%          lambda, regularization parameter
%          n, degree of polynomial
% output: estimate, same structure as LinRegressRegul

% Note that we are considering a multi-dimensional case

% new regressor matrix

% number of samples
numSamples = size(x,1);

% dimension
dim = size(x,2);

x2 = ones(numSamples,1);
for p = 1:n
    C = unique(nchoosek(rectpulse(1:dim,p),p),'rows');
    % for each possible combination of order p (n>=1)
    for i = 1:size(C,1)
        x2 = [x2 prod(x(:,C(i,:)),2)];
    end
end
m = LinRegressRegul(x2,y,lambda);
m.n = n;

end