function [x,y] = chirpData(N,var)
%[x,y] =chirpData(N,var) This function generates data for project 1. You
% indicate the number of data with N, and var indicates the variance of the
% additive Gaussian white noise. 

 x=rand(N,1)*10;

 

y=10*(.5-3.2*(x/10)+2.1*(x/10).^2+1*sin(x.^2)./(1+x.^2))+sqrt(var)*randn(N,1);
end

