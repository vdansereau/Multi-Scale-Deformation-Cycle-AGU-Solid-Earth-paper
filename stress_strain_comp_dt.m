%%This script plots, for different simulations with different time steps
(i) the macro stress vs strain (or time), 
(ii) the macro damage increment vs strain (or time),
(iii) the integrated brittle energy release vs strain (or time), normalized by the time step,
(iiii) the pdf of the brittle energy release, normalized by the time step.
It is used to produce figure 6 and figure C1: convergence.pdf and convergence2.pdf
 

clear all;
path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_1/comp_dt/alpha_4/th_10_9/';


%% Choose the time step (dt) range according to the value of De used in the simulations
%dt = [10^4 10^5 10^6 10^7 10^8];   % We = 0.001: Delta t = 10^4 s, 10^5 s, 10^6 s, 10^7 s, 10^8 s
dt = [10^3 10^4 10^5 10^6 10^7];   % We = 0.1  : Delta t = 10^3 s, 10^4 s, 10^5 s, 10^6 s, 10^7 s
%dt = [10^2 10^3 10^4 10^5 10^6];    % We = 10   : Delta t = 10^2 s, 10^3 s, 10^4 s, 10^5 s, 10^6 s


%% Set the colormap
filenb = [1 2 3 4 5];
%cc=hsv(length(filenb));

%Custom map for 5 different data, to avoid having a yellow curve
cc = [1      0                   0; ...
      1      0.800000011920929   0; ...  
      0      0.800000011920929   0.400000005960464; ...
      0      0.4                 1.0000; ...    
      0.8000 0                   1.0000];



%% Set the figures and figure axes
figure1 = figure;
    axes1 = axes('Parent', figure1, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio', [3 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'FontSize', 12); 
    xlim([0 0.00005]);    
    ylim([0 3e-3]);
    box('on');
    hold('all');
ylabel('damage rate', 'FontSize', 12);
xlabel('time', 'FontSize', 12);
 

figure2 = figure;
    axes2 = axes('Parent',figure2, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio',[3 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ... %'YTick', [0.0e-7 0.5e-7 1.0e-7 1.5e-7 2.0e-7],...
    'FontSize', 12); 
    xlim([0 0.00005]);    
    %ylim([1e-12 1e2]);
    ytickformat('%.1f');
    box('on');
    hold('all');
ylabel('\it \Sigma \rm', 'FontSize', 12); 
xlabel('time', 'FontSize', 12);


figure3 = figure;
    axes3 = axes('Parent',figure3, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio',[3 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'FontSize', 12); 
    xlim([0 0.00005]);    
    ylim([0 2e-7]);
    box('on');
    hold('all');    
ylabel('\it E_{brit} \rm', 'FontSize', 12);
xlabel('time', 'FontSize', 12); 


figure4 = figure;
    axes4 = axes('Parent',figure4, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio',[3 1 1], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'FontSize', 12); 
    xlim([0 0.00005]);    
    %ylim([1e-20 1e-10]);
    box('on');
    hold('all');  
ylabel('\it \epsilon_{x_top} \rm', 'FontSize', 12);
xlabel('time', 'FontSize', 12); 


figure5 = figure;
axes5 = axes('Parent',figure5);
pbaspect([1 2 1]); %For We = 0.001 complete figure
%pbaspect([1 1 1]); %For We = 0.1, 10 Appendix figure
set(axes5,'XMinorTick','on', 'Xlim', [10^(-8) 2*10^(-7)], 'FontSize', 8);
ylabel('PDF', 'FontSize', 8); 
xlabel('\it E_{brit} \rm', 'FontSize', 8);
hold('all');
edges = [10^(-8) : 0.2*10^(-8) : 10^(-6)];
%mini_E_brit = 10^(-8); %minimum E_brit value
%maxi_E_brit = 10^(-6); %maximum E_brit
%nbin_E_brit = (log10(maxi_E_brit)-log10(mini_E_brit))*50;



% figure6 = figure;
% axes6 = axes('Parent', figure6, 'FontSize', 12);
% edges2 = [-8 : 0.1 : -2];
% set(axes6,'XMinorTick','on','YScale','log', 'Xlim', [10^(-8) 10^(-2)], 'xtick', [10^(-8) 10^(-6) 10^(-4) 10^(-2)], 'Ylim', [0 0.5]);
% ylabel('Probability', 'FontSize', 12); 
% xlabel('dam rate', 'FontSize', 12);
% hold('all');
% 
% 
% figure7 = figure;
% axes7 = axes('Parent', figure7, 'FontSize', 12);
% edges3 = [-2 : 0.01 : 1];
% set(axes7,'XMinorTick','on','YScale','log', 'XScale', 'log', 'Xlim', [10^(-2) 10^(1)], 'Ylim', [0 0.5]);
% ylabel('Probability', 'FontSize', 12); 
% xlabel('avalanche area', 'FontSize', 12);
% hold('all');



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
macro.E_brit(find(macro.E_brit == 0)) = NaN;
macro.area        = data(12,ib);        %Cumulated area of each avalanche
element.first_dam = find(data(2,ib), 1, 'first');
f                 = 1/(macro.dt);
data = [];


%% Find first peak in stress to cut t-series after the first macro-rupture
[pks, locs] = findpeaks(macro.sigma_x_top);



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
plot(macro.time, macro.eps_x_top, 'Color', cc(i,:), 'LineWidth', 1.5, 'Displayname', ['\Delta\it t \rm = ', num2str(dt(i))], 'Parent', axes4);
hold on;


%%
plot(macro.time, macro.dam_rate,'Color', cc(i,:), 'LineWidth', 1.5, 'Displayname', ['\Delta\it t \rm = ', num2str(dt(i))], 'Parent', axes1);
hold on;


%%   
plot(macro.time, macro.sigma_x_top, 'Color', cc(i,:), 'LineWidth', 1.5, 'Displayname', ['\Delta\it t \rm = ', num2str(dt(i))], 'Parent', axes2);
hold on;


%%
plot(macro.time, macro.E_brit./macro.dt, 'Color', cc(i,:), 'LineWidth', 1.5, 'Displayname', ['\Delta\it t \rm = ', num2str(dt(i))], 'Parent', axes3);
hold on;


%% For PDF calculations, cut t-series after the first macro-rupture (first few cycles)
macro.dam_rate(1:locs(2)) = [];
macro.E_brit(1:locs(2)) = [];

macro.dam_rate(find(macro.dam_rate == 0)) = NaN;




histogram(macro.E_brit./macro.dt, edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', cc(i,:), 'LineWidth', 1.5, 'Displayname', ['\Delta\it t \rm = ', num2str(dt(i))], 'Parent', axes5);
hold on;


%%
%histogram(macro.dam_rate, 10.^edges2, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', cc(i,:), 'LineWidth', 1.0, 'Displayname', ['\Delta\it t \rm = ', num2str(dt(i))], 'Parent', axes6);
%hold on;


%%
%histogram(macro.area, 10.^edges3, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', cc(i,:), 'LineWidth', 1.0, 'Displayname', ['\Delta\it t \rm = ', num2str(dt(i))], 'Parent', axes7);
%hold on;



end

lgd = legend(axes1,'show');
set(lgd, 'Location', 'northeast', 'FontSize', 10);
hold on;
saveas(figure1, [path, 'dam_rate.pdf']);

%lgd = legend(axes2,'show');
%set(lgd, 'Location', 'northeast', 'FontSize', 10);
%hold on;
saveas(figure2, [path, 'stress_strain.pdf']); 

%lgd = legend(axes3,'show');
%set(lgd, 'Location', 'northeast', 'FontSize', 10);
%hold on;
saveas(figure3, [path, 'E_brit.pdf']);

% lgd = legend(axes4,'show');
% set(lgd, 'Location', 'northeast', 'FontSize', 10);
% hold on;
saveas(figure4, [path, 'eps_x_top.pdf']);

% lgd = legend(axes5,'show');
% set(lgd, 'Location', 'northwest', 'FontSize', 12);
% hold on;
saveas(figure5, [path, 'Ebrit_pdf.pdf']);

% lgd = legend(axes6,'show');
% set(lgd, 'Location', 'southwest', 'FontSize', 12);
% hold on;
% saveas(figure6, [path, 'damrate_hist.pdf']);

