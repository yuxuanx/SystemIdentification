dbstop if error

clear;clc

N = 50;
u = randi(2,[N,1])-1;
e = randn(N,1);

a = [1,1]; na = length(a);
b = 1; nb = length(b); 
nk = 2;
y = zeros(N,1);
for i = 1:N
    ytemp = zeros(na,1); idx = i-na:i-1; ytemp(idx>0) = y(idx(idx>0));
    utemp = zeros(nb,1); idx = i-nb-nk+1:i-nk; utemp(idx>0) = u(idx(idx>0));
    y(i) = -flip(a)*ytemp + flip(b)*utemp + e(i);
end

m = arxfit([y;u],[na,nb,nk]);
tfsys = id2tf(m)

horizon = 3;
idcompare([y;u],m,horizon)