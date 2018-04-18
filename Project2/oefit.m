function m = oefit(z,n,opt)
%OE model

% extract parameters
nb = n(1);
nf = n(2);
nk = n(3);

na = 4*nf;
nb = 4*nb;

mArx = arxfit(z,[na nb nk]);

ys = idsimulate(mArx,z);

switch(opt)
    case 'approximate'
        % approximated method
        m = arxfit([ys;z(length(z)/2+1:end)],[nf nb/4 nk]);
    case 'optimal'
        % optimal method
        % First, pass the simulated output to a FIR filter with mArx.theta(1:na)
        N = length(z)/2;
        u = z(N+1:end); % original input
        yfir = zeros(N,1);
        ufir = zeros(N,1);
        for i = 1:N
            ytemp = zeros(na,1); idx = i-na:i-1; ytemp(idx>0) = ys(idx(idx>0));
            yfir(i) = flip(mArx.theta(1:na))'*ytemp;
            utemp = zeros(na,1); idx = i-na:i-1; utemp(idx>0) = u(idx(idx>0));
            ufir(i) = flip(mArx.theta(1:na))'*utemp;
        end
        % Second, estimate a new ARX model using the filtered signals and the original model order
        m = arxfit([yfir;ufir],[nf nb/4 nk]);
end

m.model = 'OE';

end

