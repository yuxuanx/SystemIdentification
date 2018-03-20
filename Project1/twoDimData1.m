function [x,y] = twoDimData1(N,var)
%[x,y] =twoDimData1(N,var). This function generates data for project 1. You
% indicate the number of data with N, and var indicates the variance of the
% additive Gaussian white noise.


    x=rand(N,2)*10;

    X=[ones(N,1) x x.^2 x(:,1).*x(:,2)];
    theta=[-2.5 0.4 -0.5 0.02 0.06, -0.02];
    y=X*theta'+sqrt(var)*randn(N,1);
end

