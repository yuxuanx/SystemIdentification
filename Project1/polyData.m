function [x,y] = polyData(N,var,sym)
%[x,y] =polyData(N,var), or polyData(N,var,sym) This function generates data for project 1. You
% indicate the number of data with N, and var indicates the variance of the
% additive Gaussian white noise. An optional argument sym can be set to 1
% to obtain unsymmetric data.

if nargin==2, sym=0; end;

if sym==1,
  x=[rand(round(0.9*N),1)*5;5+rand(round(0.1*N),1)*5];
else
 x=rand(N,1)*10;
end

y=10*(.5-3.2*(x/10)+1*(x/10).^2+5*(x/10).^3-4*(x/10).^4)+sqrt(var)*randn(N,1);
end

