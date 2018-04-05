%% Validate problem 1.a
clear;clc; close all;
dbstop if error


% TEST OF VARIANCE
x = [1;2;3];
y = [1.1;1.8;3.1];
variance_calculated = 0.06*[7/3, -1; -1, 1/2];  % calculated by hand
m = polyfit(x,y,0,1);
assert(m.variance, variance_calculated, 1e-5);

numSamples = 50;
numVariables = 1;
dimension = 1;
x = 10*rand(numSamples,numVariables);

% add some noise
stdNoise = 1;
noise = stdNoise*randn(numSamples,dimension);
y = 2+2*x + noise;

% Validation estimate and its standard deviation
% m_LR = LinRegress(x,y);
m_LR = polyfit(x,y,0,1);
thetaEstimate = m_LR.theta;

if dimension == 1
    standardDeviation = sqrtm(inv(x'*x))*stdNoise;
    standardDeviationEstimate = sqrtm(m_LR.variance);
end




% plot a linear regression model with uncertainty region
plotModel(y,x,m_LR);
% plot a linear regression and a KNN model
k = round(sqrt(numSamples));
m_KNN = knnRegressor(x,y,k);
plotModel(y,x,m_LR,m_KNN);

%% Validate problem 1.d
clear;clc
dbstop if error
% quadratic function
numSamples = 25;
numVariables = 1;
dimension = 1;
lambda = 0;
x = 10*rand(numSamples,numVariables);

% add some noise
stdNoise = 1;
noise = stdNoise*randn(numSamples,dimension);

n = 2;
y_2ndOrder = 1 + 5*x + 0.5*x.^2 + noise;
m = polyfit(x,y_2ndOrder,lambda,n);
thetaEstimate2nd = m.theta;
plotModel(y_2ndOrder,x,m)

n = 3;
y_3rdOrder = 1 + 5*x + 0.5*x.^2 - 0.1*x.^3 + noise;
m2 = polyfit(x,y_3rdOrder,lambda,2);
m3 = polyfit(x,y_3rdOrder,lambda,n);
thetaEstimate3rd = m3.theta;
plotModel(y_3rdOrder,x,m2, m3)
