function m = LinRegress(x,y)
%Linear regression function
% input: x, regressor matrix, one row for each sample
%        y, output values, one row for each sample
% output: estimation structure
%         m.theta, parameters
%         m.variance, parameter uncertainties

% Check if x and y have the same number of rows
numSamples = size(x,1);
if numSamples == size(y,1)
    % This is a linear regression model
    m.model = 'LR';
    % Least mean squares estimate
    xSquare = x'*x;
    m.theta = xSquare\x'*y;
    % Check if y has more than one output
    if size(y,2) > 1
        m.variance = 'None';
    else
        estimateError = (y-x*m.theta);
        estimateSquareErrorSum = estimateError'*estimateError;
        m.variance = inv(xSquare)/(numSamples-size(x,2))*estimateSquareErrorSum;
    end
else
    error('Matrix x and y should have the same number of rows!')
end

end

