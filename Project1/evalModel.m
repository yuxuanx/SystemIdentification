function y_pred = evalModel(m, x)
%Evaluate model
% input: m, regression model
%        x, a matrix of regressors
% output: y, predicted values

function D = pdist2( X, Y )
Yt = Y';
XX = sum(X.*X,2);
YY = sum(Yt.*Yt,1);
D = bsxfun(@plus,XX,YY)-2*X*Yt;
end

switch(m.model)
    case 'LR'
        if isfield(m, 'n')
            x2 = polyinput(x, m.n);
        else
            x2 = polyinput(x, 1);
        end
        y_pred = x2*m.theta;
    case 'KNN'
        d = pdist2(x, m.x);
        [ds, I] = sort(d, 2);
        y_pred = mean(m.y(I(:,1:m.n)), 2);
     otherwise
         error(sprintf('The evalModel function is not implemented for model %s.', m.model))
  end

%switch(m.model)
%    case 'KNN'
%        numRegressorKNN = size(m.x,1);
%        numRegressorInput = size(x,1);
%        y_pred = zeros(numRegressorInput,size(m.y,2));
%        for i = 1:numRegressorInput
%            dist = sum((repmat(x(i,:),numRegressorKNN,1)-m.x).^2,2);
%            [~,sortpos] = sort(dist,'ascend');
%            y_pred(i,:) = mean(m.y(sortpos(1:m.n),:));
%        end
%    otherwise
%        if isfield(m,'n')
%            m_tmp = polyfit(x,[],m.lambda,m.n);
%            x = m_tmp.x;
%        end
%        y_pred = x*m.theta;
%end


end

