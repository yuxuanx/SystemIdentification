function m = arxfit(z,n)

%ARX estimator
% y(t)=-a_1*y(t−1)-...-a_na*y(t−na)+b_1*u(t−nk)+...+b_nb*u(t−nb−nk+1)+e(t)

% Input: z -- data set, z = [y(1)...y(N),u(1)...u(N)]
%        x -- [-y(t-1)...-y(t-na),u(t-nk)...u(t−nb−nk+1)]
%        n = [na,nb,nk]
%        na -- number of past outputs
%        nb -- number of past inputs
%        nk -- delays in the input signal

% Output: m.theta -- [a_1...a_na,b_1...b_nb]
%         m.model = 'ARX'

% extract parameters
na = n(1);
nb = n(2);
nk = n(3);

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

m.theta = X\y;
m.n = n;
m.model = 'ARX';

end

