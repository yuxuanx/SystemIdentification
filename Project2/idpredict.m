function prediction = idpredict(model,z,horizon)

%K-step predictor, k = horizon

if strcmp(model.model,'OE')
    prediction = idsimulate(model,z);
    return;
end

% extract parameters
na = model.n(1);
nk = model.n(3);
theta = model.theta;

% size of data set
N = length(z)/2;
% extract outputs
y = z(1:N);
% extract inputs
u = z(N+1:end);

% K-step ahead prediction only use outputs up to time t-k, and all inputs;
% \hat{y}(t|t-k) = \bar{H}_k(q)B(q)u(t) + (1-\bar{H}_k(q)A(q))y(t)
% \bar{H}_k(q) = \sum_{\tau=0}^{k-1}h(\tau)q^{-\tau}
% H(q) = 1/A(q) = 1/(1+a_1q^-1+...+a_naq^-na) = Taylor expansion = ...
% = \sum_{i=0}^Inf(a_1q^-1+...+a_naq^-na)^i = \frac{(a_1q^-1+...+a_naq^-na)^n}...
% {a_1q^-1+...+a_naq^-na-1}

%ARX estimator
% y(t)=-a_1*y(t−1)-...-a_na*y(t−na)+b_1*u(t−nk)+...+b_nb*u(t−nb−nk+1)+e(t)
% A(q)y(t) = B(q)u(t-nk) + e(t)

% Step 1, determine polynomial coefficients of \bar{H}_k, i.e.,
% coefficients of order 0 to order k-1
polynomial = poly2sym([flip(theta(1:na));0]); % x = q^-1
HqTaylorCoeffs = zeros(horizon,horizon); % store coefficients of Hq (Taylor)
for i = 1:horizon
    temp = sym2poly((-polynomial)^(i-1));
    if length(temp) < horizon
        HqTaylorCoeffs(i,end-length(temp)+1:end) = temp;
    elseif length(temp) > horizon
        HqTaylorCoeffs(i,:) = temp(end-horizon+1:end);
    else
        HqTaylorCoeffs(i,:) = temp;
    end 
end
Hk_bar_coeffs = sum(HqTaylorCoeffs); % store coefficients of \bar{H}_k from order 0 to k-1
A_poly = poly2sym([flip(theta(1:na));1]);
B_poly = poly2sym([flip(theta(na+1:end));zeros(nk,1)]);
H_k_poly = poly2sym(Hk_bar_coeffs);

theta_y = -flip(sym2poly(1-H_k_poly*A_poly));
theta_y = theta_y(2:end);
theta_y(1:horizon-1) = 0;

theta_u = flip(sym2poly(H_k_poly*B_poly));
theta_u = theta_u(nk+1:end);

% construct regressors (output parts + input parts)

na = length(theta_y);
nb = length(theta_u);
X = zeros(N,na+nb);
xtemp1 = repmat(0:-1:1-na,N,1) + (0:N-1)';
for i = 1:N
    for j = 1:na
        if xtemp1(i,j) > 0
            X(i,j) = -y(xtemp1(i,j));
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

prediction = X*[theta_y theta_u]';

end

