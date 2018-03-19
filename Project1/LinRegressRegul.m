function m = LinRegressRegul(x,y,lambda)
%Linear regression function with regularization
% input: x, regressor matrix, one row for each sample
%        y, output values, one row for each sample
%        lambda, regularization parameter
% output: estimation structure
%         m.theta, parameters
%         m.variance, parameter uncertainties

% Regularization
x2 = [x;sqrt(lambda*eye(size(x)))];
y2 = [y;zeros(size(y))];

% Do linear regression
m = LinRegress(x2,y2);
m.model = 'LRR';

end

