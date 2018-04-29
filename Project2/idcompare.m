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

theta = model.theta;
na = model.n(1);
nk = model.n(3);
A = [1;theta(1:na)]';
B = [zeros(nk,1);theta(na+1:end)]';
sys = idpoly(A,B);
sys = setcov(sys,model.variance);

% figure
% h = iopzplot(sys);
% showConfidence(h,2);
% 
% figure
% h = bodeplot(sys);
% showConfidence(h,2);
% 
% figure
% h = nyquistplot(sys);
% showConfidence(h,2);



end

