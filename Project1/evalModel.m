function predictions = evalModel(model,regressors)
%Evaluate model
% input: model
%        a matrix of regressors
% output: prediction of the output for input regressors


switch(model.model)
    case 'KNN'
        numRegressorKNN = size(model.x,1);
        numRegressorInput = size(regressors,1);
        predictions = zeros(numRegressorInput,size(model.y,2));
        for i = 1:numRegressorInput
            dist = sum((repmat(regressors(i,:),numRegressorKNN,1)-model.x).^2,2);
            [~,sortpos] = sort(dist,'ascend');
            predictions(i,:) = mean(model.y(sortpos(1:model.n),:));
        end
    otherwise
        if isfield(model,'n')
            m = ployfit(regressors,[],model.lambda,model.n);
            regressors = m.x;
        end
        predictions = regressors*model.theta;
end


end

