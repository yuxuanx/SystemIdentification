function [x,y] = twoDimData2(N,var)
%[x,y] =twoDimData2(N,var). This function generates data for project 1. You
% indicate the number of data with N, and var indicates the variance of the
% additive Gaussian white noise.


    x=rand(N,2)*10;


    y=5*sin(x(:,1).*x(:,2)/10)+sqrt(var)*randn(N,1);
end

