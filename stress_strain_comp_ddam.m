%%This script plots 
%(i) the macro damage increment vs strain (or time) 
%(ii) the macro shear stress vs strain (or time)
%(iii) the macro Ebrit vs strain (or time), normalized by the time step
%(iiii) the PDF of the macro Ebrit (normalized by the time step)
%(iiiii) the PDF of the macro damage increment

%for simulations using different delta d values.
%It is used to produce figure 9 and 10: comp_ddam_We_XX.pdf

clear all;
%% Set the path to the output files: determines the value of De that is used
%path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_001/comp_ddam/alpha_8/';
path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_1/comp_ddam/alpha_8/';

%% Set the range of delta d values
ddam = [0.1 0.3 0.5 0.7 0.9];
filenb = [1 2 3 4 5];

%% Set the colormap
%cc=parula(length(filenb)); 
%Custom colormap for 5 different data
cc = [1      0                   0; ...
      1      0.800000011920929   0; ...  
      0      0.800000011920929   0.400000005960464; ...
      0      0.4                 1.0000; ...    
      0.8000 0                   1.0000];


%% Set the figures and figure axes
figure1 = figure;
    axes1 = axes('Parent', figure1, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio', [3 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'FontSize', 10); 
    xlim([0 0.000005]);    
    ylim([0 3e-3]);
    box('on');
    hold('all');
ylabel('damage rate', 'FontSize', 10);
xlabel('time', 'FontSize', 10);
 

figure2 = figure;
    axes2 = axes('Parent',figure2, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio',[3 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'FontSize', 10); 
    xlim([0 0.000005]);    
    ylim([0 6e-7]);
    box('on');
    hold('all');
ylabel('\it \Sigma \rm', 'FontSize', 10); 
xlabel('time', 'FontSize', 10);


figure3 = figure;
    axes3 = axes('Parent',figure3, 'XGrid', 'on', 'YGrid', 'on', 'YScale', 'log', 'PlotBoxAspectRatio',[3 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'FontSize', 10); 
    xlim([0 0.000005]);    
    %ylim([1e-20 1e-10]);
    box('on');
    hold('all');    
ylabel('\it E_{brit}/ \Delta t \rm', 'FontSize', 10);
xlabel('time', 'FontSize', 10);


figure4 = figure;
axes4 = axes('Parent', figure4, 'FontSize', 16);
%edges = [-10 : 0.02 : -6];
mini_E_brit = 10^(-8); %minimum dam_rate value
maxi_E_brit = 10^(-6); %maximum dam_rate value
nbin_E_brit = (log10(maxi_E_brit)-log10(mini_E_brit))*30;
set(axes4,'XMinorTick','on','XScale','log', 'YScale', 'log', 'Xlim', [10^-8 10^-6], 'Ylim', [10^2 10^8]);
ylabel('PDF', 'FontSize', 16); 
xlabel('\it E_{brit}/ \Delta t \rm', 'FontSize', 16);
hold('all');


figure5 = figure;
axes5 = axes('Parent', figure5, 'FontSize', 20);
%edges2 = [-8 : 0.1 : -2];
mini_dam_rate = 10^(-8); %minimum dam_rate value
maxi_dam_rate = 10^(-2); %maximum dam_rate value
nbin_dam_rate = (log10(maxi_dam_rate)-log10(mini_dam_rate))*10;
set(axes5,'XMinorTick','on','YScale','log', 'XScale', 'log', 'Xlim', [10^(-8) 10^(-2)], 'xtick', [10^(-8) 10^(-6) 10^(-4) 10^(-2)], 'Ylim', [10^(-3) 10^(8)]);
ylabel('PDF', 'FontSize', 20); 
xlabel('dam rate', 'FontSize', 20);
hold('all');


%% Upload output files
for i = length(filenb) : -1 : 1
data = load([path, 'test_', int2str(filenb(i)), '.txt'])';

global macro element

[macro.t,ia,ic]= unique(data(1,:));  %model iteration
macro.dt          = macro.t(2)-macro.t(1); %model time step (s)
ib = find((macro.t(2:end)-macro.t(1:end-1)) < macro.dt*10);
macro.time        = data(1,ib);
macro.dam_rate    = data(2,ib);         %avalanche count
macro.k           = data(3,ib);         %nb of sub-iterations
macro.eps_x_top   = data(4,ib);         %x-displacement of the top boundary
macro.eps_x_mid   = data(5,ib);         %x-displacement of the mid boundary
macro.eps_x_bot   = data(6,ib);         %x-displacement of the bottom boundary
macro.De_min      = data(7,ia);        %Minimum local Deborah number
macro.sigma_x_top = data(8,ib);         %tangential stress integrated on the top boundary
macro.sigma_x_mid = data(9,ib);         %tangential stress integrated on the mid boundary
macro.sigma_x_bot = data(10,ib);        %tangential stress integrated on the bottom boundary
macro.E_brit      = data(11,ib);        %Energy dissipated in brittle failure
macro.area        = data(12,ib);        %Cumulated area of each avalanche
element.first_dam = find(data(2,ib), 1, 'first');
f                 = 1/(macro.dt);
data = [];


%%
plot(macro.time, macro.dam_rate,'Color', cc(i,:), 'LineWidth', 1.0, 'Displayname', ['\delta\it d \rm = ', num2str(ddam(i))], 'Parent', axes1);
hold on;


%%   
plot(macro.time, macro.sigma_x_top, 'Color', cc(i,:), 'LineWidth', 1.0, 'Displayname', ['\delta\it d \rm = ', num2str(ddam(i))], 'Parent', axes2);
hold on;


%%
plot(macro.time, macro.E_brit./macro.dt, 'Color', cc(i,:), 'LineWidth', 1.5, 'Displayname', ['\delta\it d \rm = ', num2str(ddam(i))], 'Parent', axes3);
hold on;



%% Cut t-series after the first macro-rupture (first few cycles)
[pks, locs] = findpeaks(macro.sigma_x_top);
macro.dam_rate(1:locs(2)) = [];
macro.E_brit(1:locs(2)) = [];

macro.dam_rate(find(macro.dam_rate == 0)) = NaN;
macro.E_brit(find(macro.E_brit == 0)) = NaN;

%%
%histogram(macro.E_brit./macro.dt, 10.^edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', cc(i,:), 'LineWidth', 1.5, 'Displayname', ['\delta\it d \rm = ', num2str(ddam(i))], 'Parent', axes4);
%hold on;

PowerPDF3(macro.E_brit./macro.dt, nbin_E_brit, mini_E_brit, maxi_E_brit);

plot(macro.bin_center, macro.pdf, 'LineWidth', 1.5, 'Color', cc(i,:), 'Displayname', ['\delta\it d \rm = ', num2str(ddam(i))], 'Parent', axes4);
hold on;
macro.pdf = [];
macro.bin_center = [];

%% Not used: PowerPDF3 function used to get a curve instead of a (bar) histogram
%histogram(macro.dam_rate, 10.^edges2, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', cc(i,:), 'LineWidth', 1.5, 'Displayname', ['\delta\it d \rm = ', num2str(ddam(i))], 'Parent', axes5);
%hold on;

PowerPDF3(macro.dam_rate, nbin_dam_rate, mini_dam_rate, maxi_dam_rate);

plot(macro.bin_center, macro.pdf, 'LineWidth', 1.5, 'Color', cc(i,:), 'Displayname', ['\delta\it d \rm = ', num2str(ddam(i))], 'Parent', axes5);
hold on;
macro.pdf = [];
macro.bin_center = [];

end


lgd = legend(axes1,'show');
set(lgd, 'Location', 'northeast', 'FontSize', 10);
hold on;
saveas(figure1, [path, 'time_strain.pdf']);

%lgd = legend(axes2,'show');
%set(lgd, 'Location', 'northeast', 'FontSize', 10);
%hold on;
saveas(figure2, [path, 'stress_strain.pdf']); 

lgd = legend(axes3,'show');
set(lgd, 'Location', 'northeast', 'FontSize', 10);
hold on;
saveas(figure3, [path, 'E_brit.pdf']);

lgd = legend(axes4,'show');
set(lgd, 'Location', 'northwest', 'FontSize', 12);
hold on;
saveas(figure4, [path, 'Ebrit_hist.pdf']);

%lgd = legend(axes5,'show');
%set(lgd, 'Location', 'southwest', 'FontSize', 12);
%hold on;
saveas(figure5, [path, 'damrate_hist.pdf']);


