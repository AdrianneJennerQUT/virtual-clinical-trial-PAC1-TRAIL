% plotting Gompertzian fit to data

data2 = [10000	10000	10000	10000	10000	10000
20000	13300	20000	15600	6670	13300
26700	17800	24400	26700	20000	17800
68900	44400	75600	66700	20000	24400
84400	102000	91100	57800	24400	26700
136000	191000	218000	182000	109000	133000
338000	309000	358000	358000	253000	200000
378000	511000	407000	438000	324000	367000];


r = 0.052696;
K = 1745423909.2602;
C01 = 15000;
C02 = 10000;
C03 = 5000;

tvec = linspace(0,7,100);
G2 = K*exp(log(C02/K)*exp(-r*tvec));


figure
hold on 
errorbar(0:7,mean(data2'),std(data2'),'LineWidth',2,'CapSize',12,'Color',[7 59 76]/255)
plot(tvec, G2,':','Color',[17 138 178]/255,'LineWidth',3)
xlabel('Time (days)')
ylabel('KGN cell count')
set(gca,'FontSize',18)
legend('Cell count data','Gompertz growth fit')
xlim([0 7])

