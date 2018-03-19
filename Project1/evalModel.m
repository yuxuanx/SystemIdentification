function predictions = evalModel(model,regressors)
%Evaluate model
% input: model
%        a matrix of regressors
% output: prediction of the output for input regressors

predictions = regressors*model;

end

