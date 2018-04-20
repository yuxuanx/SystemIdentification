dbstop if error

%% ARX model
clear;clc

N = 100;
u = randi(2,[N,1])-1;
e = 0.5*randn(N,1);

a = [-1.5 0.7]; na = length(a);
b = [1 0.5]; nb = length(b); 
nk = 0;
y = zeros(N,1);
for i = 1:N
    ytemp = zeros(na,1); idx = i-na:i-1; ytemp(idx>0) = y(idx(idx>0));
    utemp = zeros(nb,1); idx = i-nb-nk+1:i-nk; utemp(idx>0) = u(idx(idx>0));
    y(i) = -flip(a)*ytemp + flip(b)*utemp + e(i);
end

w = 0;
% w = 0.1*randn(N,1); % add white noise to test robustness
m = arxfit([y+w;u],[na,nb,nk]);
tfsys = id2tf(m)
ltiview(tfsys) % only works for causal system

horizon = 2;
idcompare([y;u],m,horizon);

%% OE model
clear;clc

N = 1000;
u = randi(2,[N,1])-1;
e = 0.5*randn(N,1);

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
w = 0;
% w = 0.1*randn(N,1); % add white noise to test robustness
m = oefit([y+w;u],[nf,nb,nk],opt);

tfsys = id2tf(m)
plottype = 'impulse'; % impulse response
% plottype = 'pzmap'; % pole/zero map
ltiview(plottype,tfsys) % only works for causal system

horizon = 2;
idcompare([y;u],m,horizon);

%% System Identification
% First system
clear;clc
load('exercise1.mat')

na = 7;
nb = 9;
nk = 0;
mArx = arxfit([y;u],[na,nb,nk]);

tfsys = id2tf(mArx)
% plottype = 'impulse'; % impulse response
plottype = 'pzmap'; % pole/zero map
ltiview(plottype,tfsys) % only works for causal system

horizon = 1;
[prediction, simulation, groundTruth] = idcompare([y;u],mArx,horizon);

immseArxPre = immse(prediction,groundTruth);
immseArxSim = immse(simulation,groundTruth);

nf = 2;
nb = 2;
nk = 0;
mOe = oefit([y;u],[nf,nb,nk],'approximate');

tfsys = id2tf(mOe)
% plottype = 'impulse'; % impulse response
plottype = 'pzmap'; % pole/zero map
ltiview(plottype,tfsys) % only works for causal system

horizon = 1;
[prediction, simulation, groundTruth] = idcompare([y;u],mOe,horizon);

immseOePre = immse(prediction,groundTruth);
immseOeSim = immse(simulation,groundTruth);

%% Second system
clear;clc
load('exercise2.mat')

na = 5;
nb = 7;
nk = 0;
mArx = arxfit([z1(:,1);z1(:,2)],[na,nb,nk]);

tfsys = id2tf(mArx)
plottype = 'impulse'; % impulse response
% plottype = 'pzmap'; % pole/zero map
ltiview(plottype,tfsys) % only works for causal system

horizon = 2;
[prediction, simulation, groundTruth] = idcompare([z2(:,1);z2(:,2)],mArx,horizon);

immseArxPre = immse(prediction,groundTruth);
immseArxSim = immse(simulation,groundTruth);

nf = 1;
nb = 3;
nk = 0;
mOe = oefit([z1(:,1);z1(:,2)],[nf,nb,nk],'optimal');

tfsys = id2tf(mOe)
plottype = 'impulse'; % impulse response
% plottype = 'pzmap'; % pole/zero map
ltiview(plottype,tfsys) % only works for causal system

horizon = 2;
[prediction, simulation, groundTruth] = idcompare([z2(:,1);z2(:,2)],mOe,horizon);

immseOePre = immse(prediction,groundTruth);
immseOeSim = immse(simulation,groundTruth);
