%%This script plots 
%(i) the macro damage increment vs strain (or time) 
%(ii) the macro shear stress vs strain (or time)
%(iii) the macro Ebrit vs strain (or time), normalized by the time step
%(iiii) the minimum De over the domain vs strain (or time)
%for simulations using different Th values.
%It is used to produce figure 8: T_h


clear all;
%path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_001/comp_th/dt_10_5/';
%path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_1/comp_th/dt_10_4/';
path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_10/comp_th/dt_10_3/';


%% Set the range of (dimensional) healing times
th = [10^(8) 10^(9) 10^(10) 10^(11)]; %s

%% Set the colormap
filenb = [1 2 3 4];
%cc=hsv(length(th));
%Custom map for 5 different data, to avoid having a yellow curve
cc = [1      0.800000011920929   0; ...  
      0      0.800000011920929   0.400000005960464; ...
      0      0.4                 1.0000; ...    
      0.8000 0                   1.0000];


%% Set the figures and figure axes
figure1 = figure;
    axes1 = axes('Parent', figure1, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio', [2 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'FontSize', 10); 
    xlim([0 0.0005]);    
    %ylim([1e-12 1e2]);
    box('on');
    hold('all');
ylabel('macro damage increment', 'FontSize', 12);
xlabel('time', 'FontSize', 12);
 

figure2 = figure;
    axes2 = axes('Parent',figure2, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio',[2 1 1], ...
    'XMinorTick', 'on', 'xtick', [0 1*10^(-4) 2*10^(-4) 3*10^(-4) 4*10^(-4) 5*10^(-4)], 'ytick', [0 1*10^(-7) 2*10^(-7) 3*10^(-7)], 'YMinorTick', 'on', ...
    'FontSize', 10); 
    xlim([0 0.0005]);    
    %ylim([1e-12 1e2]);
    box('on');
    hold('all');
ylabel('\it \Sigma \rm', 'FontSize', 12); 
xlabel('time', 'FontSize', 12);


figure3 = figure;
    axes3 = axes('Parent',figure3, 'XGrid', 'on', 'YGrid', 'on', 'YScale', 'log', 'PlotBoxAspectRatio',[2 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'FontSize', 10); 
    xlim([0 0.0005]);    
    %ylim([1e-20 1e-10]);
    box('on');
    hold('all');
    
ylabel('\it E_{brit}/\Delta t \rm', 'FontSize', 12);
xlabel('time', 'FontSize', 12); 


figure4 = figure;
    axes4 = axes('Parent',figure4, 'XGrid', 'on', 'YGrid', 'on', 'YScale', 'log', 'PlotBoxAspectRatio',[2 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'FontSize', 10); 
    xlim([0 0.0005]);    
    %ylim([1e-20 1e-10]);
    box('on');
    hold('all');
    
ylabel('De_{min}', 'FontSize', 12);
xlabel('time', 'FontSize', 12); 


%% Upload output files
for i = 1 : 1 : length(filenb)
data = load([path, 'test_', int2str(filenb(i)), '.txt'])';

global macro element

[macro.time,ia,ic]= unique(data(1,:)); %model iteration
macro.time        = macro.time - data(1,1);
macro.dam_rate    = data(2,ia);        %avalanche count
macro.k           = data(3,ia);        %nb of sub-iterations
macro.eps_x_top   = data(4,ia);        %x-displacement of the top boundary
macro.eps_x_mid   = data(5,ia);        %x-displacement of the mid boundary
macro.eps_x_bot   = data(6,ia);        %x-displacement of the bottom boundary
macro.De_min      = data(7,ia);        %Minimum local Deborah number
macro.sigma_x_top = data(8,ia);        %tangential stress integrated on the top boundary
macro.sigma_x_mid = data(9,ia);        %tangential stress integrated on the mid boundary
macro.sigma_x_bot = data(10,ia);       %tangential stress integrated on the bottom boundary
macro.E_brit      = data(11,ia);       %Energy dissipated in brittle failure
macro.area        = data(12,ia);       %Cumulated area of each avalanche
element.first_dam = find(data(2,ia), 1, 'first');
macro.dt          = macro.time(2)-macro.time(1); %model time step (s)
f                 = 1/(macro.dt);
data = [];




%% Calculate d(sigma)/dt
dsigma = diff(macro.sigma_x_top(element.first_dam:end))/macro.dt;

% Calculate skewness
%  Measure of the symmetry of a distribution
%  skewness = sum i = 1, N (Yi - mean(Y))^3/N * 1/std^3
%  zero for a normal distribution, <0 for a left-skewed dist, and >0 for right-skewed dist. 
skness = sum((dsigma - mean(dsigma)).^3.0)/length(dsigma)*(1.0/std(dsigma,1)^3.0)

% Calculate kurtosis
%  Measure of the flateness/peakedness of a distribution
%  kurtosis = sum i = 1, N (Yi - mean(Y))^4/N * 1/std^4
%  ktosis = 3 for normal dist, >3 for a peaked distribution, < 3 for a flat dist
ktosis = sum((dsigma - mean(dsigma)).^4.0)/length(dsigma)*(1.0/std(dsigma,1)^4.0)



%%
plot(macro.time, macro.dam_rate,'Color', cc(i,:), 'LineWidth', 1.0, 'Displayname', ['\it T_h \rm = ', num2str(th(i))], 'Parent', axes1);
hold on;



%%   
plot(macro.time, macro.sigma_x_top, 'Color', cc(i,:), 'LineWidth', 1.0, 'Displayname', ['\it T_h \rm = ', num2str(th(i))], 'Parent', axes2);
hold on;


%%
plot(macro.time, macro.E_brit/macro.dt, 'Color', cc(i,:), 'LineWidth', 1.0, 'Displayname', ['\it T_h \rm = ', num2str(th(i))], 'Parent', axes3);
hold on;


%%
plot(macro.time, macro.De_min, 'Color', cc(i,:), 'LineWidth', 1.0, 'Displayname', ['\it T_h \rm = ', num2str(th(i))], 'Parent', axes4);
hold on;



end

lgd = legend(axes1,'show');
set(lgd, 'Location', 'southeast', 'FontSize', 10);
hold on;
saveas(figure1, [path, 'dam_rate.pdf']);

%lgd = legend(axes2,'show');
set(lgd, 'Location', 'northeast', 'FontSize', 10);
hold on;
saveas(figure2, [path, 'stress_strain.pdf']); 

lgd = legend(axes3,'show');
set(lgd, 'Location', 'southeast', 'FontSize', 10);
hold on;
saveas(figure3, [path, 'E_brit.pdf']);

lgd = legend(axes4,'show');
set(lgd, 'Location', 'southeast', 'FontSize', 10);
hold on;
saveas(figure4, [path, 'De_min.pdf']);