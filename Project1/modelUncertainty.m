function [ymin, ymax] = modelUncertainty(m, x, n, alpha)
%Model estimation uncertainty
% input: m, regression model
%        x, regressors
%        n, number of random pertubations
%        alpha, confidence level
% output: ymin, minimum predicted values
%         ymax, maximum predicted values

d = size(m.theta, 1);
thetaRand = randn(d, n);
ys = [];
m_tmp = m; % a temporary regression model (clone of input model)
for i = 1:n 
  thetaRand_i = thetaRand(:, i);
  thetaDelta_i = thetaRand_i/norm(thetaRand_i)*chi2inv(alpha, d);
  m_tmp.theta = m.theta + m.variance^0.5*thetaDelta_i;
  ys = [ys, evalModel(m_tmp, x)];  % append pred. values with altered theta
end
ymin = min(ys, [], 2); 
ymax = max(ys, [], 2);




if 0 == 1
    predictions = evalModel(m,x);
    theta_hat = m.theta;  % parameter estimation
    variance_hat = m.variance;
    % calculate the upper and lower bound
    n = 20; % number of samples, also try {5, 10}
    confidenceInterval = 0.95;
    d = size(theta_hat,1);  % dimension of theta
    exampleThetaUp = zeros(d,n);
    exampleThetaLow = zeros(d,n);
    examplePredictionUp = zeros(size(y,1),n);
    examplePredictionLow = zeros(size(y,1),n);
    sumSquarePredictionUp = zeros(n,1);
    sumSquarePredictionLow = zeros(n,1);

    for k = 1:n
        % generate random vector
        DeltaTheta = rand(d,1);
        % scaling
        DeltaTheta = DeltaTheta*sqrt(chi2inv(confidenceInterval,d));
        % generate example vector for calculating upper/lower bound
        exampleThetaUp(:,k) = theta_hat+sqrtm(variance_hat)*DeltaTheta;
        exampleThetaLow(:,k) = theta_hat-sqrtm(variance_hat)*DeltaTheta;
        examplePredictionUp(:,k) = m.x*exampleThetaUp(:,k);
        examplePredictionLow(:,k) = m.x*exampleThetaLow(:,k);
        sumSquarePredictionUp(k) = ...
            (examplePredictionUp(:,k)-predictions)'*(examplePredictionUp(:,k)-predictions);
        sumSquarePredictionLow(k) = ...
            (examplePredictionLow(:,k)-predictions)'*(examplePredictionLow(:,k)-predictions);
    end
    
    % find upper&lower bounds
    [~,maxIndex] = max(sumSquarePredictionUp);
    [~,minIndex] = max(sumSquarePredictionLow);
    ymax = examplePredictionUp(:,maxIndex);
    ymin = examplePredictionLow(:,minIndex);
    
end