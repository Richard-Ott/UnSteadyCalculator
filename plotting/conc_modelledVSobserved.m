function H = conc_modelledVSobserved(best_pred,N10,dN10,N14,dN14)
% plots the compairson of modelled vs observed N10 and N14 concentrations
% as output from the sampler against the measured
% concentrations

hold on; box on; grid on;

nSamp = length(best_pred)/2;
N10m  = best_pred(1:nSamp);
N14m  = best_pred(nSamp+1:end);
    
cc = lines(nSamp);
H=figure();
hold on
for i = 1:nSamp
    % plot measure data
    errorbar(N10(i),N14(i),dN14(i),dN14(i),dN10(i),dN10(i),'Color',cc(i,:));

    % plot best model
    plot(N10m(i),N14m(i),'p',"Color",cc(i,:),'MarkerSize',15);

end


legend({'measured concentrations','modelled concentrations'})
ylabel('N14 at/g')
xlabel("N10 at/g")
end