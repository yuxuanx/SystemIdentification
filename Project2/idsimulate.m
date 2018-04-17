function simulation = idsimulate(model,z)
%Simulation, equal to Inf-step predictor

% extract parameters
na = model.n(1);
nk = model.n(3);
theta = model.theta;

A = [1;theta(1:na)]';
B = [zeros(nk,1);theta(na+1:end)]';
m = idpoly(A,B);
simulation = lsim(m,z(length(z)/2+1:end));

end

