%% generate linear data
N = 1000; % sample size, try 10,100,1000,10000
var = 1; % noise variance
% Remember to test the model on a second noise-free data set
[x,y] =linearData(N,0);
[x_noise,y_noise] =linearData(N,var);

% linear regression
m = ployfit(x_noise,y_noise,0,1);
m_x = ployfit(x,y,0,1); % this is one is only used for extracting x2
estimateLinear = m_x.x*m.theta;
% here, mean square error is used assess model quality
mseLinear = immse(y,estimateLinear);

% Yes, the estimate converge to the true function when the number of data
% goes to infinity

% polynomial regression
m = ployfit(x_noise,y_noise,0,5);
m_x = ployfit(x,y,0,5);
estimatePoly = m_x.x*m.theta;
msePoly = immse(y,estimatePoly);

% plot with the model quality versus number of data
N = 100*(0.5:0.5:10);
msePoly = zeros(length(N),1);
mseLinear = zeros(length(N),1);
for i = 1:length(N)
    [x,y] =linearData(N(i),0);
    [x_noise,y_noise] =linearData(N(i),var);
    
    m = ployfit(x_noise,y_noise,0,5);
    m_x = ployfit(x,y,0,5);
    estimatePoly = m_x.x*m.theta;
    msePoly(i) = immse(y,estimatePoly);
    
    m = ployfit(x_noise,y_noise,0,1);
    m_x = ployfit(x,y,0,1);
    estimateLinear = m_x.x*m.theta;
    mseLinear(i) = immse(y,estimateLinear);
end
figure; hold on
plot(N,mseLinear); plot(N,msePoly)

%% Monte Carlo simulation
numTrial = 10000;
N = 1000;
var = 1;
% assume theta is a vector with length 1; this corresponds to case 1 with
% no constant term
parameterEstimates = zeros(2,numTrial);
for i = 1:numTrial
    % maybe add a time seed here
    [x,y] =linearData(N,var);
    m = ployfit(x,y,0,1);
    parameterEstimates(:,i) = m.theta;
end

figure
subplot(2,1,1)
histogram(parameterEstimates(1,:))
xlabel('Parameter estimate'); ylabel('Frequency')
subplot(2,1,2)
histogram(parameterEstimates(2,:))
xlabel('Parameter estimate'); ylabel('Frequency')

% Question:, the parameter variance is related to x, if x changes every MC
% trial, hwo to validate? keep the x the same in every trial??

% KNN model, evaluate K with different var/N
n = 2:10; % number of neighbors
var = 0.5:0.5:2;
N = 100*(1:10);

mseKNN_nk = zeros(length(n),length(N));
mseKNN_vark = zeros(length(n),length(var));

var = 1;
for k = 1:length(n)
    for i = 1:length(N)
        [x_noise,y_noise] =linearData(N(i),var);
        m = knnRegressor(x_noise,y_noise,n(k));
        predictions = evalModel(m,x);
        mseKNN_nk(k,i) = immse(y,predictions);
    end
end

figure
mesh(100*(1:10),2:10,mseKNN_nk)
ylabel('number of neighbors'); xlabel('number of samples'); zlabel('mse')

N = 1000;
var = 0.5:0.5:2;
for k = 1:length(n)
    for i = 1:length(var)
        [x_noise,y_noise] =linearData(N,var(i));
        m = knnRegressor(x_noise,y_noise,n(k));
        predictions = evalModel(m,x);
        mseKNN_vark(k,i) = immse(y,predictions);
    end
end

figure
mesh(0.5:0.5:2,2:10,mseKNN_vark)
ylabel('number of neighbors'); xlabel('noise variance'); zlabel('mse')

%% generate polynomial data
N = 1000; % sample size, try 10,100,1000,10000
[x,y] = polyData(N,0);
var = 1;
[x_noise,y_noise] = polyData(N,var);
% linear regression
m = ployfit(x_noise,y_noise,0,1);
variance  = m.variance;
% Finding: recall the expression for the model uncertainty,
% var(noise)*inv(x'*x), the larger the sample size, the smaller inv(x'*x)

% try polynomial models of various degrees on data
n = 2:10; % polynomial degree
lambda = 0; % regularization
msePoly = zeros(length(n),1);
for i = 1:length(n)
    m = ployfit(x_noise,y_noise,0,n(i));
    m_x = ployfit(x,y,0,n(i));
    estimatePoly = m_x.x*m.theta;
    msePoly(i) = immse(y,estimatePoly);
end

% try degree 10 on data size 15, and different regularization
% if no regularization, overfitting
N = 15;
n = 10;
var = 0.1;
[x,y] = polyData(N,0);
[x_noise,y_noise] = polyData(N,var);
lambda = [0,0.01,0.1,1,10,100];

msePoly = zeros(length(lambda),1);
for i = 1:length(lambda)
    m_x = ployfit(x,y,lambda(i),n);
    m = ployfit(x_noise,y_noise,lambda(i),n);
    estimatePoly = m_x.x*m.theta;
    msePoly(i) = immse(y,estimatePoly);
end

%% unsymmetrical data
N = 100; % sample size, try 10,100,1000,10000
var = 1;
[x,y] = polyData(N,0,1);
[x_noise,y_noise] = polyData(N,var,1);
% linear regression
m_x = ployfit(x,y,0,1);
m = ployfit(x_noise,y_noise,0,1);
estimateLinear = m_x.x*m.theta;

figure; hold on
plot(y); plot(estimateLinear);

% Finding: linear model can not give a good estimate to the second half of
% the unsymmetric data due to the nonlinearity of the model
% HOW TO FIX IT? -- Use polynomial model

% KNN model, evaluate K with different var/N
n = 2:10; % number of neighbors
var = 0.5:0.5:2;
N = 100*(1:10);

mseKNN_nk = zeros(length(n),length(N));
mseKNN_vark = zeros(length(n),length(var));

var = 1;
for k = 1:length(n)
    for i = 1:length(N)
        [x_noise,y_noise] =polyData(N(i),var,1);
        m = knnRegressor(x_noise,y_noise,n(k));
        predictions = evalModel(m,x);
        mseKNN_nk(k,i) = immse(y,predictions);
    end
end

figure
mesh(100*(1:10),2:10,mseKNN_nk)
ylabel('number of neighbors'); xlabel('number of samples'); zlabel('mse')

N = 1000;
var = 0.5:0.5:2;
for k = 1:length(n)
    for i = 1:length(var)
        [x_noise,y_noise] =polyData(N,var(i),1);
        m = knnRegressor(x_noise,y_noise,n(k));
        predictions = evalModel(m,x);
        mseKNN_vark(k,i) = immse(y,predictions);
    end
end

figure
mesh(0.5:0.5:2,2:10,mseKNN_vark)
ylabel('number of neighbors'); xlabel('noise variance'); zlabel('mse')

%% generate chirp data
N = 50; % sample size, try 50,1000
var = 0.05; % noise level, try 0.05, 0.2
[x,y] = chirpData(N,0);
[x_noise,y_noise] = chirpData(N,var);
% polynomial regression (try different degrees)
n = 1:10;
msePoly = zeros(length(n),1);
for i = 1:length(n)
    m_x = ployfit(x,y,0,n(i));
    m = ployfit(x_noise,y_noise,0,n(i));
    estimatePoly = m_x.x*m.theta;
    msePoly(i) = immse(y,estimatePoly);
end

% Is it better to have a high-degree polynomial + some regularization than 
% having a lower degree polynomial without regularization?

% Do the same thing for KNN

%% Estimating two dimensional functions
% [xq,yq] = meshgrid(0:.2:10, 0:.2:10);
% vq = griddata(x(:,1),x(:,2),y,xq,yq);
% mesh(xq,yq,vq);
% hold on
% plot3(x(:,1),x(:,2),y,'o');

% try both twoDimData1&2
N = 1000;
[x,y] = twoDimData1(N,0);
N = 100;
[x_noise,y_noise] = twoDimData1(N,1);
n = 3;
m_x = ployfit(x,y,0,n); % try 2 to 5
m = ployfit(x_noise,y_noise,0,n);
estimatePoly = m_x.x*m.theta;
msePoly = immse(y,estimatePoly);

figure
[xq,yq] = meshgrid(0:.2:10, 0:.2:10);
% vq = griddata(x(:,1),x(:,2),y,xq,yq);
ve = griddata(x(:,1),x(:,2),estimatePoly,xq,yq);
% mesh(xq,yq,vq);
mesh(xq,yq,ve);
hold on
plot3(x(:,1),x(:,2),y,'o');

figure
k = 20;
m = knnRegressor(x_noise,y_noise,k);
predictions = evalModel(m,x);
mseKNN = immse(y,predictions);
ve = griddata(x(:,1),x(:,2),predictions,xq,yq);
mesh(xq,yq,ve);
hold on
plot3(x(:,1),x(:,2),y,'o');

%% Estimating ten dimensional functions