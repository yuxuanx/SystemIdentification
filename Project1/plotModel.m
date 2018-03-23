function plotModel(varargin)
%Plot data together with the estimated function

% input: generated data, regressors, model, estimated model, 
% model1,model2,model3 can be multiple, e.g.,

% extract generated data
y = varargin{1};
% extract regressors
x = varargin{2};

if nargin==3 % only one model
    % include uncertainty, using e.g., confidence interval, in the plot, 
    % can be illustrated as a shaded area
    m = varargin{3};
    modelType = m.model;
    
    model_hat = modelSelection(modelType,x,y,m);
    
    % prediction
    predictions = evalModel(model_hat,x);
    % parameter estimation
    theta_hat = model_hat.theta;
    variance_hat = model_hat.variance;
    
    % if the data is one dimension, do plotting
    if size(y,2) == 1
        % calculate the upper and lower bound
        n = 20; % number of samples, also try {5, 10}
        confidenceInterval = 0.95;
        d = size(theta_hat,1);  % dimension of theta
        exampleThetaUp = zeros(d,n);
        exampleThetaLow = zeros(d,n);
        examplePredictionUp = zeros(size(y,1),n);
        examplePredictionLow = zeros(size(y,1),n);
        sumSquarePredictionUp = zeros(n,1);
        sumSquarePredictionLow = zeros(n,1);
        
        for k = 1:n
            % generate random vector
            DeltaTheta = rand(d,1);
            % scaling
            DeltaTheta = DeltaTheta*sqrt(chi2inv(confidenceInterval,d));
            % generate example vector for calculating upper/lower bound
            exampleThetaUp(:,k) = theta_hat+sqrtm(variance_hat)*DeltaTheta;
            exampleThetaLow(:,k) = theta_hat-sqrtm(variance_hat)*DeltaTheta;
            examplePredictionUp(:,k) = x*exampleThetaUp(:,k);
            examplePredictionLow(:,k) = x*exampleThetaLow(:,k);
            sumSquarePredictionUp(k) = ...
                (examplePredictionUp(:,k)-predictions)'*(examplePredictionUp(:,k)-predictions);
            sumSquarePredictionLow(k) = ...
                (examplePredictionLow(:,k)-predictions)'*(examplePredictionLow(:,k)-predictions);
        end
        
        % find upper&lower bounds
        [~,maxIndex] = max(sumSquarePredictionUp);
        [~,minIndex] = max(sumSquarePredictionLow);
        f_max = examplePredictionUp(:,maxIndex);
        f_min = examplePredictionLow(:,minIndex);
        
        figure; hold on
        plot(y,'LineWidth',2); plot(predictions,'LineWidth',2)
        plot(f_max,'--','LineWidth',2); plot(f_min,'--','LineWidth',2);
        patch([1:length(f_max) fliplr(1:length(f_max))], [f_max' fliplr(f_min')], 'y','FaceAlpha',.1)
        legend('Generated data (Truth)','Predictions','Approximate upper bound',...
            'Approximate lower bound')
    end
else % multiple models
    % if the data is one dimension, do plotting
    if size(y,2) == 1
        numModel = nargin-2;
        predictions = zeros(size(y,1),numModel);
        figure; hold on
        plot(y,'LineWidth',2);
        Legend = cell(numModel+1,1);
        Legend{1} = 'Generated data (truth)';
        for i = 1:numModel
            % extract model, one by one
            m = varargin{i+2};
            modelType = m.model;
            model_hat = modelSelection(modelType,x,y,m);
            % prediction
            predictions(:,i) = evalModel(model_hat,x);
            plot(predictions(:,i),'LineWidth',2)
            Legend{i+1} = modelType;
        end
        legend(Legend)
    end
end

end

function model = modelSelection(modelType,x,y,m)
switch(modelType)
    case 'LR'   % linear regression
        model = LinRegress(x,y);
    case 'LRR'  % linear regression with regularization
        lambda = m.lambda; % regularization factor
        model = LinRegressRegul(x,y,lambda);
    case 'PRR'  % (multidimensional) polynomial linear regression
        lambda = m.lambda;
        n = m.n;  % polynomial order
        model = ployfit(x,y,lambda,n);
    case 'KNN' % KNN regression
        k = m.n;
        model = knnRegressor(x,y,k);
end
end

