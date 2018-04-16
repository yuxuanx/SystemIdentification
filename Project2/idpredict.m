function prediction = idpredict(model,z,horizon)

%K-step predictor, k = horizon

% extract parameters
na = model.n(1);
nb = model.n(2);
nk = model.n(3);

theta_hat_K = model.theta;
% K-step ahead prediction only use outputs up to time t-k, and all inputs;
% therefore, set parameters of y(t-k+1)...y(t-1) to zero
theta_hat_K(1:min(horizon-1,na)) = 0;

% size of data set
N = length(z)/2;
% extract outputs
y = z(1:N);
% extract inputs
u = z(N+1:end);

% construct regressors (output parts + input parts)
X = zeros(N,na+nb);
xtemp1 = repmat(0:-1:1-na,N,1) + (0:N-1)';
for i = 1:N
    for j = 1:na
        if xtemp1(i,j) > 0
            X(i,j) = y(xtemp1(i,j));
        end
    end
end
xtemp2 = repmat(1-nk:-1:1-nk-nb+1,N,1) + (0:N-1)';
for i = 1:N
    for j = 1:nb
        if xtemp2(i,j) > 0
            X(i,na+j) = u(xtemp2(i,j));
        end
    end
end

prediction = X*theta_hat_K;

end

