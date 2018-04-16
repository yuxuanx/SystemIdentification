dbstop if error

clear;clc

N = 100;
u = rand(N,1);
std = 0.01;
e = std*randn(N,1);

a = [-1.5,0.7]; na = length(a);
b = [0.2,0.1]; nb = length(b); 
nk = 1;
y = zeros(N,1);
for i = 1:N
    ytemp = zeros(na,1); idx = i-na:i-1; ytemp(idx>0) = y(idx(idx>0));
    utemp = zeros(nb,1); idx = i-nb-nk+1:i-nk; utemp(idx>0) = u(idx(idx>0));
    y(i) = -flip(a)*ytemp + flip(b)*utemp;
end
y = y + e; % ARX model, no delay for noise term

m = arxfit([y;u],[na,nb,nk]);
tfsys = id2tf(m)
% ltiview('lsim',tf(tfsys));

prediction = idpredict(m,[y;u],1);
simulation = idsimulate(m,[y;u]);

horizon = 1;
idcompare([y;u],m,horizon)