%%This script plots
%(i) the time series of x-surface displacement and damage
%(ii) the time series of macro shear stress and x-surface velocity of the upper left corner of the domain
%(iii) the time series of the post-seismic x-surface velocity of the upper left corner of the domain (log-log scale)
%It is used to produce figure 11: discussion.pdf

clear all;
%path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_001/dt_10_5/C_10_4/th_10_10/alpha_4/ddam_10/';
path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_1/dt_10_4/C_10_4/th_10_9/alpha_4/ddam_50/';
filenb = [1]; %We = 0.1
%filenb = [2]; %We = 0.001

%U0 = 10^(-9); %prescribed velocity on bottom boundary for We = 0.001
U0= 10^(-8); %prescribed velocity on bottom boundary for We = 0.1

%for i = 1 : 1 : length(filenb)
i = 1;
data = load([path, 'test_', int2str(filenb(i)), '.txt'])';

global macro element

[macro.time,ia,ic]= unique(data(1,:)); %model iteration
macro.time        = macro.time - data(1,1); %if run started from a restart
macro.dam_rate    = data(2,ia);        %avalanche count3
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


data = load([path, 'geod_', int2str(filenb(i)), '.txt'])';

element.time = unique(data(1,:));

for k = 1:1:length(element.time)
    ind = find(data(1,:)== element.time(k));
    element.top_right(k,1:7)= data(:,ind(2));             
    element.top_mid(k,1:7)  = data(:,ind(3));
    element.top_left(k,1:7) = data(:,ind(4));
    element.int_right(k,1:7)= data(:,ind(1));
    element.int_mid(k,1:7)  = data(:,ind(6));
    element.int_left(k,1:7) = data(:,ind(5));
    element.int_up(k,1:7)   = data(:,ind(7));              %upper side of the interface, at x = 0.5
    element.int_down(k,1:7) = data(:,ind(8));              %lower side of the interface, at x = 0.5
end
element.time = [];
%each element.() array contains : time(1), x(2), y(3), ux(4), uy(5), Ux(6), %Uy(7)


%% Calculate d sigma dt and find the local maxima in stress based on this

macro.dsigma = [0 diff(macro.sigma_x_top)];


% Identify threshold for the stress drop using the histogram of dsigma
%threshold = 1.0*10^(-8); %alpha = 4, ddam = 0.1, We = 0.1
%threshold = 0.5*10^(-8); %alpha = 4, ddam = 0.3, We = 0.1
threshold = 1.25*10^(-9); %alpha = 4, ddam = 0.5, We = 0.1
%threshold = 0.3*10^(-9); %alpha = 4, ddam = 0.5, We = 0.001
%threshold = 1.0*10^(-9); %alpha = 4, ddam = 0.2, 0.3 We = 0.001
%threshold = 0.5*10^(-8); %alpha = 4, ddam = 0.1, We = 0.001
%threshold = 1.0*10^(-8); %alpha = 6, ddam = 0.1, We = 0.001


pks_sigma = find((-macro.dsigma) >= threshold);

%eliminate first macro-rupture
%pks_sigma(1:3) = []; %for ddam = 0.1, eliminating 3 first cycles
pks_sigma(1) = []; %for ddam = 0.5, eliminating 1st cycle is enough



%% Plot time series of surface displacement and damage

% figure1 = figure;
%     axes1 = axes('Parent',figure1, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio',[2 1 1], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YScale', 'log', ...
%     'FontSize', 12); 
%     xlim([0 0.00001]);    
%     box('on');
%     hold('all');
%     ylabel('\it u_{sfc} \rm', 'FontSize', 12); %path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/test/';
%     xlabel('time', 'FontSize', 12);
%     plot(macro.dt*[1:1:length(element.top_left(:,4))], element.top_left(:,4), 'Color', 'k', 'LineWidth', 1.5, 'Parent', axes1);
% hold on;
%        
% saveas(figure1, [path, 'u_x_time.pdf']); 


figure1 = figure;
[AX, H1, H2] = plotyy(macro.dt*[pks_sigma(1):1:length(element.top_left(:,6))], element.top_left(pks_sigma(1):length(element.top_left(:,6)),6), macro.dt*[pks_sigma(1):1:length(element.top_left(:,6))], macro.dam_rate(pks_sigma(1):1:length(element.top_left(:,6))));
set(H1, 'Color', 'k', 'LineWidth', 1.0, 'Displayname', ['\it \epsilon_{sfc} \rm']);
set(H2, 'Color', 'c', 'LineWidth', 1.0, 'Displayname', ['macro damage incr.']);
set(AX(1), 'FontSize', 12, 'xlim', [macro.dt*pks_sigma(1) 1e-5], 'ylim', [0 1.2e-6], 'ytick', [0 0.2e-6 0.4e-6 0.6e-6 0.8e-6 1.0e-6 1.2e-6], 'PlotBoxAspectRatio',[2 1 1]);
set(AX(2), 'FontSize', 12, 'xlim', [macro.dt*pks_sigma(1) 1e-5], 'ylim', [0 0.2e-3], 'ytick', [0 0.2e-3], 'PlotBoxAspectRatio',[2 1 1]);

set(get(AX(1),'Ylabel'),'String','\it \epsilon_{sfc} \rm', 'Color', [0 0 0], 'FontSize', 12);      set(AX(1), 'YColor', [0 0 0]);
set(get(AX(2),'Ylabel'),'String','macro damage incr.', 'Color', 'c', 'FontSize', 12);       set(AX(2), 'YColor', 'c');
hold on;

xlabel('time', 'FontSize', 12);

%legend(AX(1),'show', 'Location', 'northeast');

hold off;
saveas(figure1, [path, 'disp_dam.pdf']); 




%% Scatter plot of u_sfc vs damage

figure5 = figure;

axes5 = axes('Parent',figure5, 'XScale', 'Log', 'YScale', 'Log', 'FontSize', 14, 'PlotBoxAspectRatio',[1 1 1]);

ylabel('\it u_{sfc} \rm', 'FontSize', 14); 
xlabel('macro damage incr.', 'FontSize', 14);
xlim([10^(-8) 10^(-2)]);
ylim([-10^(2) -10^(-1)]);
box('on');
hold('all');

scatter(macro.dam_rate(pks_sigma(1):1:length(element.top_left(:,4))), element.top_left(pks_sigma(1):1:length(element.top_left(:,4)),4)-U0/U0, 'MarkerEdgeColor', [0 0 0], 'Marker', '.');


saveas(figure5, [path, 'scatter_usfc_damage.pdf']);



%% Plot time series of surface displacement

% figure4 = figure;
%     axes4 = axes('Parent',figure4, 'PlotBoxAspectRatio',[2 1 1], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', ...
%     'FontSize', 12); 
%     xlim([0 0.00001]);    
%     box('on');
%     hold('all');
%     ylabel('\it \epsilon_{sfc} \rm', 'FontSize', 12); %path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/test/';
%     xlabel('time', 'FontSize', 12);
%     plot(macro.dt*[1:1:length(element.top_left(:,6))], element.top_left(:,6), 'Color', 'k', 'LineWidth', 1.5, 'Parent', axes4);
% hold on;
%        
% saveas(figure4, [path, 'eps_x_time.pdf']); 



%% Find the local minima in stress: each local min value of sigma between 2 peaks

for i = 1:length(pks_sigma)-1
    [mins, ind] = min(macro.sigma_x_top(pks_sigma(i):pks_sigma(i+1)));
    mins_sigma(i) = pks_sigma(i)+ind;
end



%% Plot stress time series with the x-velocity and an indication of each stress local maxima and u_sfc

% figure2 = figure;
%     axes2 = axes('Parent',figure2, 'PlotBoxAspectRatio',[2 1 1], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', ...
%     'FontSize', 12); 
%     xlim([0 0.00001]);    
%     ylim([0 0.0000003]);
%     box('on');
%     hold('all');
%     ylabel('\it \sigma_{N} \rm', 'FontSize', 12); 
%     xlabel('time', 'FontSize', 12);
%     plot(macro.time, macro.sigma_x_top, 'Color', 'k', 'LineWidth', 1.5, 'Parent', axes2);
% hold on;
% 
% %for i = 1: length(pks_sigma)
% %    plot([macro.time(pks_sigma(i)) macro.time(pks_sigma(i))],[0 max(macro.sigma_x_top)], 'Color', 'k', 'LineStyle', '--');
% %end
%         
% saveas(figure2, [path, 'stress.pdf']); 


figure2 = figure;
[AX, H1, H2] = plotyy(macro.dt*[pks_sigma(1):1:length(element.top_left(:,4))], macro.sigma_x_top(pks_sigma(1):1:length(element.top_left(:,4))), macro.dt*[pks_sigma(1):1:length(element.top_left(:,4))], element.top_left(pks_sigma(1):length(element.top_left(:,4)),4)-U0/U0 );
set(H1, 'Color', 'k', 'LineWidth', 1.0, 'Displayname', ['\it \sigma_{N} \rm']);
set(H2, 'Color', 'r', 'LineWidth', 1.0, 'Displayname', ['\it u_{sfc} \rm']);
set(AX(1), 'FontSize', 12, 'xlim', [macro.dt*pks_sigma(1) 1e-5], 'ylim', [0 0.0000003], 'ytick', [0 1e-7 2e-7 3e-7], 'PlotBoxAspectRatio',[2 1 1]);
set(AX(2), 'FontSize', 12, 'xlim', [macro.dt*pks_sigma(1) 1e-5], 'ylim', [-40 0], 'ytick', [-40 -20 0], 'PlotBoxAspectRatio',[2 1 1]);

set(get(AX(1),'Ylabel'),'String','\it \sigma_{N} \rm', 'Color', [0 0 0], 'FontSize', 12);      set(AX(1), 'YColor', [0 0 0]);
set(get(AX(2),'Ylabel'),'String','\it u_{sfc} \rm', 'Color', 'r', 'FontSize', 12);      set(AX(2), 'YColor', 'r');
hold on;

xlabel('time', 'FontSize', 12);

hold off;
saveas(figure2, [path, 'stress_usfc.pdf']); 



%% Divide the surface velocity time series between the peaks in sigma
%  to look at the behavior of the post-seismic sfc velocities

%Find the last pks_sigma present in the element.top_left t-series
int = find(pks_sigma-length(element.top_left) < 1, 1, 'last');
pks_sigma(int:end) = [];

cc=gray(length(pks_sigma));
figure3 = figure;
axes3 = axes('Parent',figure3, 'XGrid', 'on', 'YGrid', 'on', 'PlotBoxAspectRatio',[1 1 1], ...
'XMinorTick', 'on', 'YMinorTick', 'on', 'XScale', 'Log', 'YScale', 'Log', 'FontSize', 14); 
box('on');
ylim([0.1 100]); 
ylabel('\it u_{sfc} \rm', 'FontSize', 14); 
xlabel('time', 'FontSize', 14);
hold('all');
    
for i = 1: length(pks_sigma)-1
    
    u_sfc = element.top_left(pks_sigma(i):pks_sigma(i+1)-1, 4);
    time = macro.dt*[1:1:length(u_sfc)];
    
    plot(time, -(u_sfc-U0/U0), 'Color', cc(i,:), 'LineWidth', 1.0, 'Parent', axes3);
    hold on;
    %instead of max(u_sfc), I should try to find a constant value as the
    %max surface velocity in the post-seismic period
    
    u_sfc = [];
    time = [];

end

saveas(figure3, [path, 'u_sfc.pdf']);
