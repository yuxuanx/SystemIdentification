clear;clc
%% Validate problem 1.a
numSamples = 10;
numVariables = 3;
dimension = 1;
x = 10*rand(numSamples,numVariables);
% theta = rand(numVariables,dimension);
theta = [-3;4;2];

% add some noise
stdNoise = 1;
noise = stdNoise*randn(numSamples,dimension);
y = x*theta + noise; 

% Validation estimate and its standard deviation
m = LinRegress(x,y);
thetaEstimate = m.theta;

if dimension == 1
    standardDeviation = sqrtm(inv(x'*x))*stdNoise;
    standardDeviationEstimate = sqrtm(m.variance);
end

% plot a linear regression model with uncertainty region
plotModel(y,x,m);
% plot a linear regression and a KNN model
m2 = knnRegressor(x,y,3);
plotModel(y,x,m,m2);

%% Validate problem 1.d
% quadratic function
numSamples = 100;
numVariables = 3;
dimension = 6;
lambda = 0;
n = 3;
x = 100*rand(numSamples,numVariables);

x2 = ones(numSamples,1);
for p = 1:n
    C = unique(nchoosek(rectpulse(1:numVariables,p),p),'rows');
    % for each possible combination of order p (n>=1)
    for i = 1:size(C,1)
        x2 = [x2 prod(x(:,C(i,:)),2)];
    end
end

theta = rand(size(x2,2),dimension);

% add some noise
stdNoise = 0.1;
noise = stdNoise*randn(numSamples,dimension);
y = x2*theta + noise;

m = ployfit(x,y,lambda,n);
thetaEstimate = m.theta;
