function m = knnRegressor(x,y,n)
%KNN Regressor

m.x = x;
m.y = y;
m.n = n;
m.model = 'KNN';
m.label = sprintf('KNN %d', n);

end

