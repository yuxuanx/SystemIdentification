function idcompare(z,model,horizon)
%Plot function, compare the simulated/predicted output together with the
% true output with uncertainty region

prediction = idpredict(model,z,horizon);

simulation = idsimulate(model,z);

trueOutput = z(1:length(z)/2);

figure
hold on
plot(trueOutput,'Linewidth',2);
plot(prediction,'Linewidth',2);
plot(simulation,'Linewidth',2);
str = strcat('Predicted output, K=', num2str(horizon));
xlabel('Time step'); ylabel('Output')
legend('True output',str,'Simulated Output')


end

