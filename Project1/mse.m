function value = mse(x,y)
  e2 = (x-y).^2;
  value = mean(e2);
end