% Remember to test the model on a second noise-free data set
%% generate linear data
N = 1000; % sample size, try 10,100,1000,10000
var = 1;
[x,y] =linearData(N,var);

% linear regression
m = ployfit(x,y,0,1);
estimateLinear = m.x*m.theta;
% here, mean square error is used assess model quality
mseLinear = immse(y,estimateLinear);

% polynomial regression
m = ployfit(x,y,0,5);
estimatePoly = m.x*m.theta;
msePoly = immse(y,estimatePoly);

% Monte Carlo simulation
numTrial = 100;
% assume theta is a vector with length 1; this corresponds to case 1 with
% no constant term
parameterEstimates = zeros(1,numTrial);
for i = 1:numTrial
    % maybe add a time seed here
    [x,y] =linearData(N,var);
    m = LinRegress(x,y);
    parameterEstimates(:,i) = m.theta;
end

figure
histogram(parameterEstimates)
xlabel('Parameter estimate'); ylabel('Frequency')

% KNN model
k = 10; % number of neighbors
m = knnRegressor(x,y,k);
predictions = evalModel(m,x);

%% generate polynomial data
N = 100; % sample size, try 10,100,1000,10000
var = 1;
[x,y] = polyData(N,var);
% linear regression
m = ployfit(x,y,0,1);
estimateLinear = m.x*m.theta;
% Finding: parameter goes to zero for the linear model with increasing
% sample size; Why? recall the variance-bias tradeoff and the
% interpretation of model variance -- a measure of how much variability the
% estimate would have if we change the data used for model fitting. If we
% fit an unflexible model to a polynomial data, the bias would be large,
% which implies a small variance. The larger the sample size we have, the
% less model uncertainty we will have. This property holds for all the
% models.

% try polynomial models of various of data
n = 5; % polynomial degree
lambda = 1; % regularization
m = ployfit(x,y,0,1);

% unsymmetrical data
N = 1000; % sample size, try 10,100,1000,10000
var = 1;
[x,y] = polyData(N,var,1);
% linear regression
m = ployfit(x,y,0,1);
estimateLinear = m.x*m.theta;
% polynomial regression
m = ployfit(x,y,0,5);
estimatePoly = m.x*m.theta;
% Finding: linear model can not give a good estimate to the second half of
% the unsymmetric data due to the nonlinearity of the model
% HOW TO FIX IT??
m = ployfit(x,y,1,1);
estimateRegu = m.x*m.theta;

% KNN model
k = 10; % number of neighbors
m = knnRegressor(x,y,k);
predictions = evalModel(m,x);

%% generate chirp data
N = 1000; % sample size, try 50,1000
var = 1; % noise level, try 1, 0.1
[x,y] = chirpData(N,var);
% polynomial regression
m = ployfit(x,y,0,5);
estimatePoly = m.x*m.theta;
% KNN model
k = 10; % number of neighbors
m = knnRegressor(x,y,k);
predictions = evalModel(m,x);

%% Estimating two or higher dimensional functions

