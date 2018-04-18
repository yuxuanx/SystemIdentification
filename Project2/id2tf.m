function tfsys = id2tf(m)

%Convert an ARX/OE model to transfer function

switch(m.model)
    case 'ARX'
        A = [1 m.theta(1:m.n(1))'];
        B = [zeros(1,m.n(3)) m.theta(m.n(1)+1:end)'];
        tfsys = tf(A,B,1,'Variable','q^-1');
    case 'OE'
        F = [1 m.theta(1:m.n(1))'];
        B = [zeros(1,m.n(3)) m.theta(m.n(1)+1:end)'];
        tfsys = tf(F,B,1,'Variable','q^-1');
end

end

