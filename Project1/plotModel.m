function plotModel(varargin)
%Plot data together with the estimated function

% input: generated data, regressors, model, estimated model, 
% model1,model2,model3 can be multiple, e.g.,

% extract generated data
y = varargin{1};
% extract regressors
x = varargin{2};
[numSamplesX, numDimensionsX] = size(x);
[numSamplesY, numDimensionsY] = size(y);

if nargin==3 % only one model
    % include uncertainty, using e.g., confidence interval, in the plot, 
    % can be illustrated as a shaded area
    m = varargin{3};
    
    % if the data is one dimension, do plotting
    if numDimensionsY == 1
        alpha = 0.95;
        n = 100;
        [f_min, f_max] = modelUncertainty(m, x, n, alpha);
        
        [~,I] = sort(x);
        
        figure; hold on
        plot(x,y,'o','Markersize',8);
        predictions = evalModel(m,x);
        plot(x(I),predictions(I),'LineWidth',2)
        plot(x(I),f_max(I),'--','LineWidth',2); 
        plot(x(I),f_min(I),'--','LineWidth',2);
        %patch([x(I)', fliplr(x(I)')], [f_max(I)', fliplr(f_min(I)')], 'y', 'FaceAlpha',0.5);
        xlabel('x'); ylabel('y')
        legend('Generated data (Truth)',m.label, ...
        'Approximate upper bound', 'Approximate lower bound')
    end
else % multiple models
    % if the data is one dimension, do plotting
    [~,I] = sort(x);
    if numDimensionsY == 1
        x_model = linspace(x(I(1)), x(I(end)), 100)';
        numModel = nargin-2;
        predictions = zeros(size(x_model,1),numModel);
        %predictions = zeros(size(y,1),numModel);
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
                    predictions(:,i) = evalModel(m,x_model);
                case 'KNN' % KNN regression
                    predictions(:,i) = evalModel(m,x_model);
            end
            plot(x_model,predictions(:,i),'LineWidth',2)
            %plot(x(I),predictions(I,i),'LineWidth',2)
            Legend{i+1} = m.label;
            %Legend{i+1} = m.model;
        end
        legend(Legend)
    end
end

end

