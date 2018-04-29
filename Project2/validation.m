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
% ltiview(tfsys) % only works for causal system

horizon = 1;
idcompare([y;u],m,horizon);

%% OE model
clear;clc

N = 1000;
u = randi(2,[N,1])-1;
e = 0*randn(N,1);

f = [-1.5 0.7]; nf = length(f);
b = [1 0.5]; nb = length(b); 
nk = 1;
y = zeros(N,1);
for i = 1:N
    ytemp = zeros(nf,1); idx = i-nf:i-1; ytemp(idx>0) = y(idx(idx>0));
    utemp = zeros(nb,1); idx = i-nb-nk+1:i-nk; utemp(idx>0) = u(idx(idx>0));
    etemp = zeros(nf,1); idx = i-nf:i; etemp(idx>0) = e(idx(idx>0));
    y(i) = -flip(f)*ytemp + flip(b)*utemp + flip([1 f])*etemp;
end

utest = u(1:N/2);
ytest = y(1:N/2);

uvalidation = u(N/2+1:end);
yvalidation = y(N/2+1:end);

% opt = 'approximate';
opt = 'optimal';
w = 0;
% w = 0.1*randn(N,1); % add white noise to test robustness
mOptimal = oefit([ytest+w;utest],[nf,nb,nk],opt);

% tfsys = id2tf(mOptimal)
% plottype = 'impulse'; % impulse response
% plottype = 'pzmap'; % pole/zero map
% ltiview(plottype,tfsys) % only works for causal system

[groundTruth,simulation] = idcompare([yvalidation;uvalidation],mOptimal,[]);
immseOeSimOptimal = immse(simulation,groundTruth);

opt = 'approximate';
mApprox = oefit([ytest+w;utest],[nf,nb,nk],opt);
[groundTruth,simulation] = idcompare([yvalidation;uvalidation],mApprox,[]);
immseOeSimApprox = immse(simulation,groundTruth);

%% System Identification
% First system
clear;clc
load('exercise1.mat')

dataSize = length(y);

% identify ARX model order
V = arxstruc(iddata(y(1:dataSize/2),u(1:dataSize/2)),...
    iddata(y(dataSize/2+1:end),u(dataSize/2+1:end)),struc(1:10,1:10,1:10));
nn = selstruc(V,0);
mArx = arxfit([y;u],[3 6 1]);

tfsys = id2tf(mArx)
% plottype = 'impulse'; % impulse response
plottype = 'pzmap'; % pole/zero map
% ltiview(plottype,tfsys) % only works for causal system

horizon = 1;
[groundTruth,simulation,prediction] = idcompare([y;u],mArx,horizon);

immseArxPre = immse(prediction,groundTruth);
immseArxSim = immse(simulation,groundTruth);


% identify OE model order
V = ivstruc(iddata(y(1:dataSize/2),u(1:dataSize/2)),...
    iddata(y(dataSize/2+1:end),u(dataSize/2+1:end)),struc(1:10,1:10,1:10));
nn = selstruc(V,0);
mOe = oefit([y;u],nn,'approximate');

tfsys = id2tf(mOe)
plottype = 'impulse'; % impulse response
% plottype = 'pzmap'; % pole/zero map
ltiview(plottype,tfsys) % only works for causal system

horizon = 1;
[groundTruth,simulation] = idcompare([y;u],mOe,horizon);

immseOeSim = immse(simulation,groundTruth);

%% Second system
clear;clc
load('exercise2.mat')

% identify ARX model order
V = arxstruc(iddata(z1(:,1),z1(:,2)),iddata(z1(:,1),z1(:,2)),struc(1:20,1:20,1:20));
nn = selstruc(V,0);
mArx = arxfit([z1(:,1);z1(:,2)],nn);

tfsys = id2tf(mArx)
plottype = 'impulse'; % impulse response
% plottype = 'pzmap'; % pole/zero map
% ltiview(plottype,tfsys) % only works for causal system

horizon = 2;
[groundTruth,simulation,prediction] = idcompare([z2(:,1);z2(:,2)],mArx,horizon);

immseArxPre = immse(prediction,groundTruth);
immseArxSim = immse(simulation,groundTruth);

% identify OE model order
V = ivstruc(iddata(z1(:,1),z1(:,2)),iddata(z1(:,1),z1(:,2)),struc(1:10,1:10,1:10));
nn = selstruc(V,0);
mOe = oefit([z1(:,1);z1(:,2)],nn,'optimal');

tfsys = id2tf(mOe)
plottype = 'impulse'; % impulse response
% plottype = 'pzmap'; % pole/zero map
% ltiview(plottype,tfsys) % only works for causal system

horizon = 2;
[groundTruth,simulation] = idcompare([z2(:,1);z2(:,2)],mOe,horizon);

immseOeSim = immse(simulation,groundTruth);
