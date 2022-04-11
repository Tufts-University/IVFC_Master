voltage_405socket=[.03,.02,.03,.02,.03,.03;...
                   .07,.04,.06,.04,.07,.06;...
                   .14,.09,.11,.07,.14,.13;...
                   .28,.16,.21,.13,.25,.24;...
                   .5,.29,.38,.23,.44,.44;...
                   .83,.47,.64,.39,.73,.73];
gain=[2:0.2:3];
pmt=1:6;
%voltage_405socket=voltage_405socket./max(max(voltage_405socket));
figure; plot(gain,voltage_405socket,'.','MarkerSize',23)
boldify
xticks([2:.2,3]);
xticks([2:.2:3]);
ylabel('Normalized Voltage')
xlabel('Gain')
legend('PMT1','PMT2','PMT3','PMT4','PMT5','PMT6','location','southoutside','numcolumns',6)
set(gca,'FontSize',20)
set(gca,'yscale','log')

yyaxis('right')
 x=500:50:750;
y=0.0045*x+2.6364;
y=10.^y;
%plot(gain,y)
plot(gain,y./max(y))
set(gca,'yscale','log')
legend('PMT1','PMT2','PMT3','PMT4','PMT5','PMT6','PMT Gain Graph (Hamamatsu)','location','southoutside','numcolumns',7)