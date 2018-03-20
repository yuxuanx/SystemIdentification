function [x,y] = tenDimData(N,var)
%[x,y] =tenDimData(N,var). This function generates data for project 1. You
% indicate the number of data with N, and var indicates the variance of the
% additive Gaussian white noise.


    x=rand(N,10);

    X=[ones(N,1) x];
    % Linear part
    thetaL=[-2.5 0.4 -0.5 0.2 0.6, -0.2 0.3 0.5 -0.4 -0.3 0.7];
    y=X*thetaL'-0.3*x(:,3).^2+0.6*x(:,5).^3+0.5*x(:,2).*x(:,6).*x(:,10)-0.5*x(:,2).*x(:,4).*x(:,10)+sqrt(var)*randn(N,1);
end

