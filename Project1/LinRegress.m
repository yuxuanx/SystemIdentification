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
    inv_xSquare = inv(x'*x);
    m.theta = x\y;
%     m.theta = xSquare\x'*y;
    % Check if y has more than one output
    if size(y,2) > 1
        m.variance = 'None';
    else
        estimateError = (y-x*m.theta);
        
%         temp = estimateError'*(eye(numSamples) - x*inv_xSquare*x')*estimateError;
%         m.variance = inv_xSquare*temp/(numSamples-size(x,2));
        
        estimateSquareErrorSum = estimateError'*estimateError;
        m.variance = inv_xSquare/(numSamples-size(x,2))*estimateSquareErrorSum;
    end
    m.x = x;
    m.y = y;
elseif numSamples ~= size(y,1)
    error('Matrix x and y should have the same number of rows!')
elseif numSamples < size(x,2)
    error('Number of samples should be no smaller than the number of parameters')
end

end

