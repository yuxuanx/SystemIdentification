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
    
    % prediction
    predictions = evalModel(m,x);
    % parameter estimation
    theta_hat = m.theta;
    variance_hat = m.variance;
    
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
            examplePredictionUp(:,k) = m.x*exampleThetaUp(:,k);
            examplePredictionLow(:,k) = m.x*exampleThetaLow(:,k);
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
        
        [~,I] = sort(x);
        
        figure; hold on
        plot(x,y,'o','Markersize',8); plot(x(I),predictions(I),'LineWidth',2)
        plot(x(I),f_max(I),'--','LineWidth',2); plot(x(I),f_min(I),'--','LineWidth',2);
        patch([x(I)' fliplr(x(I)')], [f_max(I)' fliplr(f_min(I)')], 'y','FaceAlpha',.1)
        xlabel('x'); ylabel('y')
        legend('Generated data (Truth)','Estimation','Approximate upper bound',...
            'Approximate lower bound')
    end
else % multiple models
    % if the data is one dimension, do plotting
    [~,I] = sort(x);
    if size(y,2) == 1
        numModel = nargin-2;
        predictions = zeros(size(y,1),numModel);
        figure; hold on
        plot(x,y,'o','markersize',8);
        xlabel('x'); ylabel('y')
        Legend = cell(numModel+1,1);
        Legend{1} = 'Generated data (Truth)';
        for i = 1:numModel
            % extract model, one by one
            m = varargin{i+2};
            % prediction
            switch(m.model)
                case 'LR'  % (multidimensional) polynomial linear regression
                    predictions(:,i) = evalModel(m,x);
                case 'KNN' % KNN regression
                    predictions(:,i) = evalModel(m,x);
            end
            
            plot(x(I),predictions(I,i),'LineWidth',2)
            Legend{i+1} = m.model;
        end
        legend(Legend)
    end
end

end

