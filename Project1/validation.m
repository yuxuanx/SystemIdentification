%% Validate problem 1.a
clear;clc
dbstop if error

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
m_LR = ployfit(x,y,0,1);
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
numSamples = 1000;
numVariables = 1;
dimension = 1;
lambda = 0;
x = 10*rand(numSamples,numVariables);

% add some noise
stdNoise = 1;
noise = stdNoise*randn(numSamples,dimension);

n = 2;
y_2ndOrder = 1 + 0.5*x + 0.25*x.^2 + noise;
m = ployfit(x,y_2ndOrder,lambda,n);
thetaEstimate2nd = m.theta;

n = 3;
y_3rdOrder = 1 + 0.5*x + 0.25*x.^2 + 0.125*x.^3 + noise;
m = ployfit(x,y_3rdOrder,lambda,n);
thetaEstimate3rd = m.theta;
