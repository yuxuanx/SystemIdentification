function varargout = idcompare(z,model,horizon)
%Plot function, compare the simulated/predicted output together with the
% true output with uncertainty region

simulation = idsimulate(model,z);
groundTruth = z(1:length(z)/2);

if strcmp(model.model,'ARX')
    prediction = idpredict(model,z,horizon);
    varargout{1} = groundTruth;
    varargout{2} = simulation;
    varargout{3} = prediction;
    
    figure
    hold on
    plot(groundTruth,'Linewidth',2);
    plot(prediction,'--','Linewidth',2);
    plot(simulation,'-.','Linewidth',2);
    str = strcat('Predicted output, K=', num2str(horizon));
    xlabel('Time step'); ylabel('Output')
    legend('True output',str,'Simulated Output')
elseif strcmp(model.model,'OE')
    varargout{1} = groundTruth;
    varargout{2} = simulation;
    
    figure
    hold on
    plot(groundTruth,'Linewidth',2);
    plot(simulation,'-.','Linewidth',2);
    xlabel('Time step'); ylabel('Output')
    legend('True output','Simulated/Predicted Output')
end




end

