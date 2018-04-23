function simulation = idsimulate(model,z)
%Simulation, equal to Inf-step predictor

% extract parameters
na = model.n(1);
nb = model.n(2);
nk = model.n(3);
theta = model.theta;
N = length(z)/2;
u = z(N+1:end);

simulation = zeros(N,1);

for i = 1:N
    ytemp = zeros(na,1); idx = i-na:i-1; ytemp(idx>0) = simulation(idx(idx>0));
    utemp = zeros(nb,1); idx = i-nb-nk+1:i-nk; utemp(idx>0) = u(idx(idx>0));
    simulation(i) = -flip(theta(1:na))'*ytemp + flip(theta(na+1:end))'*utemp;
end

end

