dbstop if error

%% ARX model
clear;clc

N = 100;
u = randi(2,[N,1])-1;
e = 0.1*randn(N,1);

a = [-1.5 0.7]; na = length(a);
b = [1 0.5]; nb = length(b); 
nk = 0;
y = zeros(N,1);
for i = 1:N
    ytemp = zeros(na,1); idx = i-na:i-1; ytemp(idx>0) = y(idx(idx>0));
    utemp = zeros(nb,1); idx = i-nb-nk+1:i-nk; utemp(idx>0) = u(idx(idx>0));
    y(i) = -flip(a)*ytemp + flip(b)*utemp + e(i);
end

% w = 0;
w = 0.1*randn(N,1); % add white noise to test robustness
m = arxfit([y+w;u],[na,nb,nk]);
tfsys = id2tf(m)
ltiview(tfsys) % only works for causal system

horizon = 2;
idcompare([y;u],m,horizon)

%% OE model
clear;clc

N = 1000;
u = randi(2,[N,1])-1;
e = 0.1*randn(N,1);

f = [-1.5 0.7]; nf = length(f);
b = [1 0.5]; nb = length(b); 
nk = 0;
y = zeros(N,1);
for i = 1:N
    ytemp = zeros(nf,1); idx = i-nf:i-1; ytemp(idx>0) = y(idx(idx>0));
    utemp = zeros(nb,1); idx = i-nb-nk+1:i-nk; utemp(idx>0) = u(idx(idx>0));
    etemp = zeros(nf,1); idx = i-nf:i-1; etemp(idx>0) = e(idx(idx>0));
    y(i) = -flip(f)*ytemp + flip(b)*utemp + flip(f)*etemp;
end

% opt = 'approximate';
opt = 'optimal';
% w = 0; 
w = 0.1*randn(N,1); % add white noise to test robustness
m = oefit([y+w;u],[nf,nb,nk],opt);

tfsys = id2tf(m)
ltiview(tfsys) % only works for causal system

horizon = 2;
idcompare([y;u],m,horizon)

%% System Identification
% First system
load('exercise1.mat')

na = 1;
nb = 1;
nk = 0;
mArx = arxfit([y;u],[na,nb,nk]);

horizon = 2;
idcompare([y;u],mArx,horizon)

nf = 1;
nb = 1;
nk = 0;
mOe = oefit([y;u],[nf,nb,nk]);

horizon = 2;
idcompare([y;u],mOe,horizon)

% Second system
load('exercise2.mat')

na = 1;
nb = 1;
nk = 0;
mArx = oefit([z1(:,1);z1(:,2)],[na,nb,nk],'optimal');

horizon = 2;
idcompare([z2(:,1);z2(:,2)],mArx,horizon)

nf = 1;
nb = 1;
nk = 0;
mOe = oefit([z1(:,1);z1(:,2)],[nf,nb,nk],'optimal');

horizon = 2;
idcompare([z2(:,1);z2(:,2)],mOe,horizon)
