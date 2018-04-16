function simulation = idsimulate(model,z)
%Simulation, equal to Inf-step predictor

% simulation = idpredict(model,z,model.n(1)+1);

% extract parameters
na = model.n(1);
nb = model.n(2);
nk = model.n(3);

% size of data set
N = length(z)/2;
% extract inputs
u = z(N+1:end);

% construct regressors (output parts + input parts)
X = zeros(N,nb);
xtemp = repmat(1-nk:-1:1-nk-nb+1,N,1) + (0:N-1)';
for i = 1:N
    for j = 1:nb
        if xtemp(i,j) > 0
            X(i,j) = u(xtemp(i,j));
        end
    end
end

simulation = X*model.theta(na+1:end);

end

