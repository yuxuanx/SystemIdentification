function m = oefit(z,n)
%OE model

% extract parameters
nb = n(1);
nf = n(2);
nk = n(3);

na = 4*nf;
nb = 4*nb;

mArx = arxfit(z,[na nb nk]);

ys = idsimulate(mArx,z);

% approximated method
m = arxfit([ys;z(length(z)/2+1:end)],[na nb nk]);

% optimal method


m.model = 'OE';

end

