%% generate linear data
close all;
isOctave = exist('OCTAVE_VERSION', 'builtin') ~= 0;
if isOctave
    pkg load statistics
end

% 3.1 a and b

N = 1000; % sample size, try 10,100,1000,10000
var = 1; % noise variance
% Remember to test the model on a second noise-free data set
[x,y] = linearData(N,0);
[x_noise,y_noise] = linearData(N,var);

% linear regression
m = polyfit(x_noise, y_noise, 0, 1);
estimateLinear = evalModel(m,x);
% here, mean square error is used assess model quality

mseLinear = mse(y,estimateLinear);

% Yes, the estimate converge to the true function when the number of data
% goes to infinity

% polynomial regression
m = polyfit(x_noise,y_noise,0,5);
estimatePoly = evalModel(m,x);
msePoly = mse(y,estimatePoly);

% plot with the model quality versus number of data
N = 10.^[1:4];
msePoly = zeros(length(N),1);
thetaPoly = zeros(length(N),6);
mseLinear = zeros(length(N),1);
thetaLinear = zeros(length(N),2);

for i = 1:length(N)
    [x,y] = linearData(N(i),0);
    [x_noise,y_noise] = linearData(N(i),var);
    
    m = polyfit(x_noise,y_noise,0,5);
    estimatePoly = evalModel(m,x);
    msePoly(i) = mse(y,estimatePoly);
    thetaPoly(i,:) = m.theta;
    
    m = polyfit(x_noise,y_noise,0,1);
    estimateLinear = evalModel(m,x);
    mseLinear(i) = mse(y,estimateLinear);
    thetaLinear(i,:) = m.theta;
end
figure; hold on
loglog(N,mseLinear,'-o','Linewidth',2); loglog(N,msePoly,'-o','Linewidth',2)
xlabel('Sample size'); ylabel('Mean square error')
legend('constant + linear term','5th order polynomial')

% 3.1 c

%% Monte Carlo simulation
numTrial = 1000;
N = 10;
var = 1;
% assume theta is a vector with length 1; this corresponds to case 1 with
% no constant term
parameterEstimates = zeros(2,numTrial);
averagevariance = zeros(2,2,numTrial);
for i = 1:numTrial
    % maybe add a time seed here
    [x,y] =linearData(N,var);
    m = polyfit(x,y,0,1);
    parameterEstimates(:,i) = m.theta;
    averagevariance(:,:,i) = inv(m.x'*m.x);
end
averagevariance = mean(averagevariance,3);

%figure
%subplot(2,1,1)
%hist(parameterEstimates(1,:))
%xlabel('Parameter estimate (constant)'); ylabel('Frequency')
%subplot(2,1,2)
%hist(parameterEstimates(2,:))
%xlabel('Parameter estimate (linear term)'); ylabel('Frequency')

function [xx, yy] = cplot_mesh(x, y)
    X = [x(:), y(:)];
    limit_min = min(X); limit_max = max(X);
    x_min = limit_min(1); y_min = limit_min(2);
    x_max = limit_max(1); y_max = limit_max(2);
    [xx, yy] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max,100));
end

function cplot2d_MC(x, y, label_x, label_y)
    hold on; box on;
    X = [x(:), y(:)];
    mu = mean(X); sigma = cov(X);
    [xx, yy] = cplot_mesh(x, y);
    F_MC = reshape(mvnpdf([xx(:), yy(:)], mu, sigma), length(xx), length(yy));
    scatter(X(:,1), X(:,2),'.'); 
    contour(xx, yy, F_MC, 5);
    title('Monte Carlo'); xlabel(label_x); ylabel(label_y);
end

function cplot2d_theo(m, x, y, label_x, label_y)
    [xx, yy] = cplot_mesh(x, y);
    F = reshape(mvnpdf([xx(:), yy(:)], m.theta, m.variance), length(xx), length(yy));
    contour(xx, yy, F, 5);
    title('Expression'); xlabel(label_x); ylabel(label_y);
end

function cplot2d(m, x, y, label_x, label_y)
    % contour plot of montecarlo vs theoretical parameter uncertainty
    subplot(1,2,1);
    cplot2d_theo(m, x, y, label_x, label_y)
    subplot(1,2,2);
    cplot2d_MC(x, y, label_x, label_y)
    
    
    %X = [x(:), y(:)];
    %mu = mean(X); sigma = cov(X); 
    %limit_min = min(X); limit_max = max(X);
    %x_min = limit_min(1); y_min = limit_min(2);
    %x_max = limit_max(1); y_max = limit_max(2);
    %[xx, yy] = meshgrid(linspace(x_min, x_max, 100), linspace(y_min, y_max,100));
    %F = reshape(mvnpdf([xx(:), yy(:)], m.theta, m.variance), length(xx), length(yy));
    %F_MC = reshape(mvnpdf([xx(:), yy(:)], mu, sigma), length(xx), length(yy));
    
    %subplot(1,2,1); 
    
    %contour(xx, yy, F, 5);
    %title('Expression')
    %xlabel(label_x); ylabel(label_y)
    
    %subplot(1,2,2); hold on; box on;
    
    %scatter(X(:,1), X(:,2),'.'); 
    %contour(xx, yy, F_MC, 5);
    %title('Monte Carlo')
    %xlabel(label_x); ylabel(label_y)
end
figure; hold on;
cplot2d(m, parameterEstimates(1,:), parameterEstimates(2,:), 
    'Parameter estimate (constant)', 'Parameter estimate (linear)')

% Question:, the parameter variance is related to x, if x changes every MC
% trial, hwo to validate? keep the x the same in every trial??

% 3.1 d

% KNN model, evaluate K with different var/N
n = 2:10; % number of neighbors
var = 0.5:0.5:2;
N = 100*(1:10);

mseKNN_nk = zeros(length(n),length(N));
mseKNN_vark = zeros(length(n),length(var));
%var = 1;
for i = 1:length(N)
    [x_noise,y_noise] = linearData(N(i),1);
    for k = 1:length(n)
        m = knnRegressor(x_noise,y_noise,n(k));
        predictions = evalModel(m,x);
        mseKNN_nk(k,i) = mse(y,predictions);
    end
end
figure
surf(N,n,mseKNN_nk)
ylabel('number of neighbors'); xlabel('number of samples'); zlabel('mse')

%N = 1000;
%var = 0.5:0.5:2;
for i = 1:length(var)
    [x_noise,y_noise] = linearData(1000,var(i));
    for k = 1:length(n)
        m = knnRegressor(x_noise,y_noise,n(k));
        predictions = evalModel(m,x);
        mseKNN_vark(k,i) = mse(y,predictions);
    end
end
figure
surf(var,n,mseKNN_vark)
ylabel('number of neighbors'); xlabel('noise variance'); zlabel('mse')

% 3.2 a

%% generate polynomial data
N = 100; % sample size, try 10,100,1000,10000
[x,y] = polyData(10*N,0);
var = 1;
[x_noise,y_noise] = polyData(N,var);
% linear regression
m = polyfit(x_noise,y_noise,0,1);
variance  = m.variance;
% Finding: recall the expression for the model uncertainty,
% var(noise)*inv(x'*x), the larger the sample size, the smaller inv(x'*x)

% plot with the model quality versus number of data
N = 10.^[1:4];
msePoly = zeros(length(N),1);
mseLinear = zeros(length(N),1);
models_linear = [];
for i = 1:length(N)
    [x,y] = polyData(N(i),0);
    [x_noise,y_noise] = polyData(N(i),var);
    
    m = polyfit(x_noise,y_noise,0,5);
    estimatePoly = evalModel(m,x);
    msePoly(i) = mse(y,estimatePoly);
    
    m = polyfit(x_noise,y_noise,0,1);
    estimateLinear = evalModel(m,x);
    mseLinear(i) = mse(y,estimateLinear);
    models_linear = [models_linear, m];
end
figure; hold on
loglog(N,mseLinear,'-o','Linewidth',2); loglog(N,msePoly,'-o','Linewidth',2)
xlabel('Sample size'); ylabel('Mean square error')
legend('constant + linear term','5th order polynomial')

figure;
x_label = 'Parameter estimate (constant)'; y_label = 'Parameter estimate (linear)';
mu = models_linear(i).theta; sigma = models_linear(i).variance;
xx = linspace(0, 3, 100);
yy = linspace(-1.5, -0.5, 100);

for i = 1:length(models_linear)
    subplot(2,2,i); cplot2d_theo(models_linear(i), xx, yy, x_label, y_label)
end

% 3.2 b
% try polynomial models of various degrees on data
n = 1:10; % polynomial degree
lambda = 0; % regularization
msePoly = zeros(length(n),1);
models = [];
[x,y] = polyData(100,0);
[x_noise,y_noise] = polyData(100,var);
for i = 1:length(n)
    m = polyfit(x_noise,y_noise,lambda,n(i));
    estimatePoly = evalModel(m,x);
    msePoly(i) = mse(y,estimatePoly);
    models = [models, m];
end
plotModel(y, x, models(1), models(2), models(3), models(4), models(10));

figure
plot(n,msePoly,'Linewidth',2)
xlabel('Model degree'); ylabel('Mean square error')

% try degree 10 on data size 15, and different regularization
% if no regularization, overfitting
N = 15;
n = 10;
var = 1;
[x_noise,y_noise] = polyData(N, var);
[x,y] = polyData(N, 0);
lambda = 10.^[-2:1:2];
models = [];
%rand('seed',1)
msePoly = zeros(length(lambda),1);
for i = 1:length(lambda)
    m = polyfit(x_noise,y_noise,lambda(i),n);
    estimatePoly = evalModel(m, x);
    msePoly(i) = mse(y, estimatePoly);
    models = [models, m];
end
figure
semilogx(lambda,msePoly,'Linewidth',2)
xlabel('Regularization term'); ylabel('Mean square error')

plotModel(y, x, models(1), models(2), models(3), models(4), models(5))

%3.2 c and d
%% unsymmetrical data
N = 100; % sample size, try 10,100,1000,10000
var = 1;
[x_noise,y_noise] = polyData(N,var,1);
[x,y] = polyData(N,0,1);
% linear regression
m_LR = polyfit(x_noise,y_noise,0,1);
m_Poly = polyfit(x_noise,y_noise,0,4);
m_KNN = knnRegressor(x_noise,y_noise,10);
plotModel(y,x,m_LR,m_Poly,m_KNN)

m_LR = polyfit(x_noise,y_noise,2,1);
m_Poly = polyfit(x_noise,y_noise,2,4);
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
        mseKNN_nk(k,i) = mse(y,predictions);
    end
end

figure
surf(100*(1:10),2:10,mseKNN_nk)
ylabel('number of neighbors'); xlabel('number of samples'); zlabel('mse')

N = 1000;
var = 0.5:0.5:2;
for i = 1:length(var)
    [x_noise,y_noise] =polyData(N,var(i),1);
    for k = 1:length(n)
        m = knnRegressor(x_noise,y_noise,n(k));
        predictions = evalModel(m,x);
        mseKNN_vark(k,i) = mse(y,predictions);
    end
end

figure
surf(0.5:0.5:2,2:10,mseKNN_vark)
ylabel('number of neighbors'); xlabel('noise variance'); zlabel('mse')

% 3.3 a
%% generate chirp data
N = [50, 1000]; % sample size, try 50,1000
var = 0.2; % noise level, try 0.05, 0.2

% polynomial regression (try different degrees)
n = 2:20;
msePoly = zeros(length(n),length(N));
mseKNN = zeros(length(n),length(N));
for j = 1:length(N)
[x,y] = chirpData(N(j),0);
[x_noise,y_noise] = chirpData(N(j),var);
    for i = 1:length(n)
        mPoly = polyfit(x_noise,y_noise,10,n(i));
        estimatePoly = evalModel(mPoly,x);
        msePoly(i, j) = mse(y,estimatePoly);
        
        mKNN = knnRegressor(x_noise,y_noise, n(i));
        estimateKNN = evalModel(mKNN, x);
        mseKNN(i, j) = mse(y, estimateKNN);
    end
end
figure; hold on; box on;
semilogy(n, msePoly); semilogy(n, mseKNN); 
xlabel('Model degree'); ylabel('mse');
legend('poly (50 samples)', 'poly (1000 samples)','knn (50 samples)', 'knn (1000 samples)',...
'location','NorthWest')

%figure; hold on;
%plot(y,x,mPoly, mKNN)

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
    mse_KNN(k) = mse(y,predictions);
end
figure; hold on
semilogy(n,mseKNN)
semilogy(n,msePoly)

%% Estimating two dimensional functions
% [xq,yq] = meshgrid(0:.2:10, 0:.2:10);
% vq = griddata(x(:,1),x(:,2),y,xq,yq);
% mesh(xq,yq,vq);
% hold on
% plot3(x(:,1),x(:,2),y,'o');

% Try to fix the testing data, that is when degree, number of neighbors
% changes, the test, validation data remains the same

% try both twoDimData1&2

function poly_surf(data_handle, var, lambda, poly_degree)
    N = 100;
    rand('seed', 1);
    [x_noise,y_noise] = data_handle(N, var);
    [x,y] = data_handle(10*N,0);
    n = poly_degree;
    m = polyfit(x_noise,y_noise,lambda,n);
    estimatePoly = evalModel(m,x);
    mse = mse(y,estimatePoly);
    
    [xq,yq] = meshgrid(0:.2:10, 0:.2:10);
    vq = griddata(x(:,1),x(:,2),y,xq,yq);
    mesh(xq,yq,vq);
    hold on;

    plot3(x(:,1),x(:,2),estimatePoly,'o');
    legend('True function',sprintf('Estimation (mse=%.3f)', mse))
    xlabel('Dimension 1');ylabel('Dimension 2'),zlabel('y')
end

function knn_surf(data_handle, var, n_neighbors)
    N = 100;
    k = n_neighbors;
    rand('seed', 1);
    [x_noise,y_noise] = data_handle(N, var);
    [x,y] = data_handle(10*N,0);
    m = knnRegressor(x_noise,y_noise,k);
    predictions = evalModel(m,x);
    mse = mse(y,predictions);
    [xq,yq] = meshgrid(0:.2:10, 0:.2:10);
    ve = griddata(x(:,1),x(:,2),y,xq,yq);
    mesh(xq,yq,ve);
    hold on;
    plot3(x(:,1),x(:,2),predictions,'o');
    legend('True function',sprintf('Estimation (mse=%.3f)', mse))
    xlabel('Dimension 1');ylabel('Dimension 2'),zlabel('y')
end

function poly_mse(data_handle, lambda, degree)
    N = 100;
    rand('seed', 1);
    [x_noise,y_noise] = data_handle(N, 1);
    [x,y] = data_handle(10*N,0);
    n = degree;
    msePoly = zeros(1, length(n));
    for i = 1:length(n)
        m = polyfit(x_noise, y_noise, lambda, n(i));
        predictions = evalModel(m,x);
        msePoly(i) = mse(y,predictions);
    end
  plot(n, msePoly); xlabel('polynomial degree'); ylabel('mse');
end
figure;
poly_mse(@twoDimData1, 0, 2:5); title('twoDimData1');
figure;
poly_mse(@twoDimData2, 0, 2:5); title('twoDimData2')

function knn_mse(data_handle, neighbors)
    N = 100;
    rand('seed', 1);
    [x_noise,y_noise] = data_handle(N, 1);
    [x,y] = data_handle(10*N,0);
    k = neighbors;
    mseKNN = zeros(1, length(k));
    for i = 1:length(k)
        m = knnRegressor(x_noise,y_noise,k(i));
        predictions = evalModel(m,x);
        mseKNN(i) = mse(y,predictions);
    end
  plot(k, mseKNN); xlabel('number of neighbors'); ylabel('mse');
end
figure;
knn_mse(@twoDimData1, 1:20)
figure;
knn_mse(@twoDimData2, 1:20)

% polynomial without regularization
figure;
poly_surf(@twoDimData1, 1, 0, 2); title('twoDimData1 (poly)')
figure;
poly_surf(@twoDimData2, 1, 0, 4); title('twoDimData2 (poly)')

figure;
knn_surf(@twoDimData1, 1, 8); title('twoDimData1 (KNN)')
figure;
knn_surf(@twoDimData2, 1, 3); title('twoDimData2 (KNN)')


%% Estimating ten dimensional functions
N = 1000;
[x_noise,y_noise] = tenDimData(N,1);
[x,y] = tenDimData(10*N,0);
lambda = 10.^[-2:4];
mse_1 = zeros(length(lambda),1);
mse_3 = zeros(length(lambda),1);
for i = 1:length(lambda)
    m_1 = polyfit(x_noise, y_noise, lambda(i), 1);
    y_pred = evalModel(m_1, x);
    mse_1(i,1) = mse(y,y_pred);

    m_3 = polyfit(x_noise, y_noise, lambda(i), 3);
    y_pred = evalModel(m_3,x);
    mse_3(i,1) = mse(y,y_pred);
end
figure; hold on; box on;
semilogx(lambda, mse_1); semilogx(lambda, mse_3);
legend('degree 1', 'degree 3'); xlabel('regularization parameter'); ylabel('mse'); 

k=10:10:100;
mseKNN = zeros(length(k), 1);
for i = 1:length(k)
    %k = ceil(sqrt(N));
    m = knnRegressor(x_noise,y_noise,k(i));
    predictions = evalModel(m,x);
    mseKNN(i,1) = mse(y,predictions);
end
figure; hold on; box on;
plot(k, mseKNN);
xlabel('number of neighbors'); ylabel('mse'); 
