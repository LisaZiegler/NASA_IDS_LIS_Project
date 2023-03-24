subplot(5,1,1)
%plot((11).date, EE61_variables(11).value, '-k') %chla
%hold on
plot(meteordata1.VarName1, meteordata1.VarName6, '-','linewidth',2)
ylabel('Wind Speed (m/s)')
xlim([0 350]);
hold on
subplot(5,1,2)
plot(stations.variables(1).date, stations.variables(1).value, 'o','markersize',10) %chla
%hold on
%plot(mean_EE61_variables(1).dates, mean_EE61_variables(1).values, 'or', 'linewidth',1)
ylabel('chla (ug chl/L)')
hold on
subplot(5,1,3)
plot(stations.variables(18).date, stations.variables(18).value, 'o','markersize',10) %chla
hold on
%plot(mean_EE61_variables(1).dates, mean_EE61_variables(1).values, 'or', 'linewidth',1)
ylabel('DO2 (mg O2/l)')
hold on
subplot(5,1,4)
plot(stations.variables(7).date, stations.variables(7).value, 'o','markersize',10) %chla
%hold on
%plot(mean_EE61_variables(1).dates, mean_EE61_variables(1).values, 'or', 'linewidth',1)
%plot(mean_EE61_variables(5).dates, mean_EE61_variables(5).values, 'or','linewidth',1)
ylabel('NH4 (mg N/l)')
xlim([0 350])
hold on
subplot(5,1,5)
plot(stations.variables(9).date, stations.variables(9).value, 'o','markersize',10) %chla
hold on
%plot(mean_EE61_variables(1).dates, mean_EE61_variables(1).values, 'or', 'linewidth',1)
ylabel('NO3 (mg N/l)')
hold on
subplot(9,1,5)
plot(EE61_variables(7).date, EE61_variables(7).value, '-k') % PO4
hold on
plot(mean_EE61_variables(7).dates, mean_EE61_variables(7).values, 'or','linewidth',1)
hold on
ylabel('PO4 (mg P/l)')
hold on
subplot(9,1,7)
plot(EE61_variables(4).date, EE61_variables(4).value, '-k') % DOP
hold on
plot(mean_EE61_variables(4).dates, mean_EE61_variables(4).values, 'or','linewidth',1)
ylim([0 0.15])
ylabel('DOP (mg P/l)')
hold on
subplot(9,1,8)
plot(EE61_variables(3).date, EE61_variables(3).value, '-k') % DON
hold on
plot(mean_EE61_variables(3).dates, mean_EE61_variables(3).values, 'or','linewidth',1)
ylabel('DON (mg N/l)')
hold on
subplot(9,1,9)
plot(EE61_variables(8).date, EE61_variables(8).value, '-k') % PON
hold on
plot(mean_EE61_variables(8).dates, mean_EE61_variables(8).values, 'or','linewidth',1)
ylabel('PON (mg N/l)')

legend('obs','mean','Orientation','horizontal','location','southoutside')