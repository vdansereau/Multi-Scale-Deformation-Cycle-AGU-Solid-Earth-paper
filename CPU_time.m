%% This code plots the CPU time calculated for simulations using 
% We = 0.001, alpha = 4, ddam = 0.1, th = 10^10 s and different t-steps
% Used to produce figure 6: CPU_time

path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_001/CPU_time/';


dt = [10^4 10^5 10^6 10^7 10^8]./10^15; %s
CPU = [90133.8 10592.1 1382.21 328.46 74.5621]; %s 


figure1 = figure;
axes1 = axes('Parent',figure1, 'FontSize', 14);
box(axes1,'on');
set(axes1,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log', 'XGrid', 'on', 'YGrid','on');
xlabel('\it \Delta t', 'FontSize', 14);
ylabel('CPU time','FontSize', 14);

hold(axes1,'on');

loglog(dt,CPU, 'Color', 'k', 'LineWidth', 1.0, 'Marker', '.', 'MarkerSize', 15);

saveas(figure1, [path, 'CPU_time.pdf']);