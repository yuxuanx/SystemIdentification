function tfsys = id2tf(arx)

%Convert an ARX model to transfer function
A = [1 arx.theta(1:arx.n(1))'];
B = [zeros(1,arx.n(3)) arx.theta(arx.n(1)+1:end)'];
tfsys = tf(A,B,1,'Variable','q^-1');

end

