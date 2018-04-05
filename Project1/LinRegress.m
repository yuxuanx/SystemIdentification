function m = LinRegress(x,y)
%Linear regression function
% input: x, regressor matrix, one row for each sample
%        y, output values, one row for each sample
% output: estimation structure
%         m.theta, parameters
%         m.variance, parameter uncertainties

if isempty(y)
    m.x = x;
    return;
end

% Check if x and y have the same number of rows

[numSamplesX, numDimsX] = size(x);
[numSamplesY, numDimsY] = size(y);
if numSamplesX == numSamplesY
    % This is a linear regression model
    m.model = 'LR';
    m.label = m.model;
    % Least mean squares estimate
    m.theta = x\y;
%     m.theta = (x'*x)\x'*y;
    % Check if y has more than one output
    if numDimsY > 1
        m.variance = 'None';
    else
        estimateError = (y-x*m.theta);  % Should use the evalModel function
        %estimateError = (y-evalModel(m, x));
        inv_xSquare = inv(x'*x);
        estimateSquareErrorSum = estimateError'*estimateError;
        m.variance = inv_xSquare/(numSamplesX-numDimsX)*estimateSquareErrorSum;
    end
    m.x = x;
    m.y = y;
elseif numSamplesX ~= numSamplesY
    error('Matrix x and y should have the same number of rows!')
elseif numSamplesX < numDimsX
    error('Number of samples should be no smaller than the number of parameters')
end

end

