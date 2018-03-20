function [x,y] = linearData(N,var)
%[x,y] =linearData(N,var) This function generates data for project 1. You
% indicate the number of data with N, and var indicates the variance of the
% additive Gaussian white noise.

x=rand(N,1)*10;

y=1.5+0.5*x+sqrt(var)*randn(N,1);
end

