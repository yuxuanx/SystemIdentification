function m = LinRegressRegul(x,y,lambda)
%Linear regression function with regularization
% input: x, regressor matrix, one row for each sample
%        y, output values, one row for each sample
%        lambda, regularization parameter
% output: estimation structure
%         m.theta, parameters
%         m.variance, parameter uncertainties

% Regularization
if lambda == 0
    x2 = x;
    y2 = y;
else
    x2 = [x;sqrt(lambda*eye(size(x)))];
    y2 = [y;zeros(size(y))];
end

% Do linear regression
m = LinRegress(x2,y2);
m.x = x;
m.lambda = lambda;
m.label = sprintf('LRR %.3f', lambda);

end

