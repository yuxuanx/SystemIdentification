%% generate linear data
N = 1000; % sample size, try 10,100,1000,10000
var = 1; % noise variance
% Remember to test the model on a second noise-free data set
[x,y] =linearData(10*N,0);
[x_noise,y_noise] =linearData(N,var);

% linear regression
m = ployfit(x_noise,y_noise,0,1);
estimateLinear = evalModel(m,x);
% here, mean square error is used assess model quality
mseLinear = immse(y,estimateLinear);

% Yes, the estimate converge to the true function when the number of data
% goes to infinity

% polynomial regression
m = ployfit(x_noise,y_noise,0,5);
estimatePoly = evalModel(m,x);
msePoly = immse(y,estimatePoly);

% plot with the model quality versus number of data
N = 100*(0.5:0.5:10);
msePoly = zeros(length(N),1);
mseLinear = zeros(length(N),1);
for i = 1:length(N)
    [x,y] =linearData(N(i),0);
    [x_noise,y_noise] =linearData(N(i),var);
    
    m = ployfit(x_noise,y_noise,0,5);
    estimatePoly = evalModel(m,x);
    msePoly(i) = immse(y,estimatePoly);
    
    m = ployfit(x_noise,y_noise,0,1);
    estimateLinear = evalModel(m,x);
    mseLinear(i) = immse(y,estimateLinear);
end
figure; hold on
plot(N,mseLinear,'-o','Linewidth',2); plot(N,msePoly,'-o','Linewidth',2)
xlabel('Sample size'); ylabel('Mean square error')
legend('constant + linear term','5th order polynomial')

%% Monte Carlo simulation
numTrial = 1000;
N = 1000;
var = 1;
% assume theta is a vector with length 1; this corresponds to case 1 with
% no constant term
parameterEstimates = zeros(2,numTrial);
averagevariance = zeros(2,2,numTrial);
for i = 1:numTrial
    % maybe add a time seed here
    [x,y] =linearData(N,var);
    m = ployfit(x,y,0,1);
    parameterEstimates(:,i) = m.theta;
    averagevariance(:,:,i) = inv(m.x'*m.x);
end
averagevariance = mean(averagevariance,3);

figure
subplot(2,1,1)
histogram(parameterEstimates(1,:))
xlabel('Parameter estimate (constant)'); ylabel('Frequency')
subplot(2,1,2)
histogram(parameterEstimates(2,:))
xlabel('Parameter estimate (linear term)'); ylabel('Frequency')

% Question:, the parameter variance is related to x, if x changes every MC
% trial, hwo to validate? keep the x the same in every trial??

% KNN model, evaluate K with different var/N
n = 2:10; % number of neighbors
var = 0.5:0.5:2;
N = 100*(1:10);

mseKNN_nk = zeros(length(n),length(N));
mseKNN_vark = zeros(length(n),length(var));

var = 1;
for i = 1:length(N)
    [x_noise,y_noise] =linearData(N(i),var);
    for k = 1:length(n)
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
for i = 1:length(var)
    [x_noise,y_noise] =linearData(N,var(i));
    for k = 1:length(n)
        m = knnRegressor(x_noise,y_noise,n(k));
        predictions = evalModel(m,x);
        mseKNN_vark(k,i) = immse(y,predictions);
    end
end

figure
mesh(0.5:0.5:2,2:10,mseKNN_vark)
ylabel('number of neighbors'); xlabel('noise variance'); zlabel('mse')

%% generate polynomial data
N = 100; % sample size, try 10,100,1000,10000
[x,y] = polyData(10*N,0);
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
    m = ployfit(x_noise,y_noise,lambda,n(i));
    estimatePoly = evalModel(m,x);
    msePoly(i) = immse(y,estimatePoly);
end

figure
plot(n,msePoly,'Linewidth',2)
xlabel('Model degree'); ylabel('Mean square error')

% try degree 10 on data size 15, and different regularization
% if no regularization, overfitting
N = 15;
n = 10;
var = 1;
[x,y] = polyData(10*N,0);
[x_noise,y_noise] = polyData(N,var);
lambda = [0.01,0.1,1,10,100];

msePoly = zeros(length(lambda),1);
for i = 1:length(lambda)
    m = ployfit(x_noise,y_noise,lambda(i),n);
    estimatePoly = evalModel(m,x);
    msePoly(i) = immse(y,estimatePoly);
end
figure
semilogx(lambda,msePoly,'Linewidth',2)
xlabel('Regularization term'); ylabel('Mean square error')

%% unsymmetrical data
N = 100; % sample size, try 10,100,1000,10000
var = 1;
[x_noise,y_noise] = polyData(N,var,1);
[x,y] = polyData(10*N,0,1);
% linear regression
m_LR = ployfit(x_noise,y_noise,0,1);
m_Poly = ployfit(x_noise,y_noise,0,4);
m_KNN = knnRegressor(x_noise,y_noise,10);
plotModel(y,x,m_LR,m_Poly,m_KNN)

% Finding: linear model can not give a good estimate to the second half of
% the unsymmetric data due to the nonlinearity of the model
% HOW TO FIX IT? -- Use polynomial model, or unparametric model, e.g., KNN

% KNN model, evaluate K with different var/N
n = 2:10; % number of neighbors
var = 0.5:0.5:2;
N = 100*(1:10);

mseKNN_nk = zeros(length(n),length(N));
mseKNN_vark = zeros(length(n),length(var));

var = 1;
for i = 1:length(N)
    [x_noise,y_noise] =polyData(N(i),var,1);
    for k = 1:length(n)
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
for i = 1:length(var)
    [x_noise,y_noise] =polyData(N,var(i),1);
    for k = 1:length(n)
        m = knnRegressor(x_noise,y_noise,n(k));
        predictions = evalModel(m,x);
        mseKNN_vark(k,i) = immse(y,predictions);
    end
end

figure
mesh(0.5:0.5:2,2:10,mseKNN_vark)
ylabel('number of neighbors'); xlabel('noise variance'); zlabel('mse')

%% generate chirp data
N = 100; % sample size, try 50,1000
var = 0.2; % noise level, try 0.05, 0.2
[x,y] = chirpData(10*N,0);
[x_noise,y_noise] = chirpData(N,var);
% polynomial regression (try different degrees)
n = 2:20;
msePoly = zeros(length(n),1);
for i = 1:length(n)
    m = ployfit(x_noise,y_noise,10,n(i));
    estimatePoly = evalModel(m,x);
    msePoly(i) = immse(y,estimatePoly);
end

% Is it better to have a high-degree polynomial + some regularization than 
% having a lower degree polynomial without regularization?

% This is an ambuguious question, if there is a lower degree polynomial 
% model that we used caused large bias, and we increased the degree and got
% a giher degree polynomial that casued smaller bias. If the increasing
% variance brought by the increasing of the degree can be reduced by
% introducing regularization, there is no reason why we should not use the
% high degree one with regularization

% Do the same thing for KNN

% subplot for var = 0.05, 0.2, N = 50, 1000, k = 2:10
n = 2:20;
mseKNN = zeros(length(n),1);
for k = 1:length(n)
    m = knnRegressor(x_noise,y_noise,n(k));
    predictions = evalModel(m,x);
    mseKNN(k) = immse(y,predictions);
end
figure; hold on
plot(n,mseKNN)
plot(n,msePoly)

%% Estimating two dimensional functions
% [xq,yq] = meshgrid(0:.2:10, 0:.2:10);
% vq = griddata(x(:,1),x(:,2),y,xq,yq);
% mesh(xq,yq,vq);
% hold on
% plot3(x(:,1),x(:,2),y,'o');

% Try to fix the testing data, that is when degree, number of neighbors
% changes, the test, validation data remains the same

% try both twoDimData1&2
N = 100;
[x_noise,y_noise] = twoDimData2(N,1);
[x,y] = twoDimData2(10*N,0);
n = 1;
m = ployfit(x_noise,y_noise,0,n);
estimatePoly = evalModel(m,x);
msePoly = immse(y,estimatePoly);

figure
[xq,yq] = meshgrid(0:.2:10, 0:.2:10);
vq = griddata(x(:,1),x(:,2),y,xq,yq);
mesh(xq,yq,vq);
hold on

plot3(x(:,1),x(:,2),estimatePoly,'o');
legend('True function','Estimation')
xlabel('Dimension 1');ylabel('Dimension 2'),zlabel('y')

figure
k = 4;
m = knnRegressor(x_noise,y_noise,k);
predictions = evalModel(m,x);
mseKNN = immse(y,predictions);
ve = griddata(x(:,1),x(:,2),y,xq,yq);
mesh(xq,yq,ve);
hold on
plot3(x(:,1),x(:,2),predictions,'o');
legend('True function','Estimation')
xlabel('Dimension 1');ylabel('Dimension 2'),zlabel('y')

%% Estimating ten dimensional functions
N = 500;
[x_noise,y_noise] = tenDimData(N,1);
[x,y] = tenDimData(10*N,0);
n = 3;
m = ployfit(x_noise,y_noise,10,n);
estimateLinear = evalModel(m,x);
mseLinear = immse(y,estimateLinear);
k = ceil(sqrt(N));
m = knnRegressor(x_noise,y_noise,k);
predictions = evalModel(m,x);
mseKNN = immse(y,predictions);