%%This script plots the damage rate and brittle energy released PSDs
%for simulations using different Th values.
%It is used to produce figure 8: T_h

clear all;

%% Set the range of (dimensional) healing times
th = [10^8 10^9 10^10 10^11]; %(s)

%% Set the colormap
%cc=hsv(length(th));
%Custom map for 5 different data, to avoid having a yellow curve
cc = [1      0.800000011920929   0; ...  
      0      0.800000011920929   0.400000005960464; ...
      0      0.4                 1.0000; ...    
      0.8000 0                   1.0000];


%% Set the min and max cutoffs of the time series (remove first few cycles because outliers)
And the running mean window for the PSDs. 
min_cutoff = 5; %nb of cycles removed at the beginning of the t-series
max_cutoff = 500000; %max length of the t-series
window = 5; %frequency-averaging window (odd number) for running mean

%% Pre-allocate arrays for results
mean_rm_PDS_dam = NaN(length(th), max_cutoff/2);
mean_rm_PDS_Ebrit = NaN(length(th), max_cutoff/2);


%% Upload output files
%path_ = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_001/comp_th/dt_10_5/';   
%path_ = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_1/comp_th/dt_10_4/';
path_ = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_10/comp_th/dt_10_3/';


for m = 1: 1: length(th)    
%--------------------------------------------------------------------------

%Five simulations started with different noise on C for each mean PSD
filenb = [1 2 3 4 5];

if m == 1  
    path = [path_, 'th_10_8/alpha_4/ddam_10/'];
elseif m == 2  
    path = [path_, 'th_10_9/alpha_4/ddam_10/'];
elseif m == 3
    path = [path_, 'th_10_10/alpha_4/ddam_10/'];
else
    path = [path_, 'th_10_11/alpha_4/ddam_10/'];
end


for i = 1 : 1 : length(filenb)
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
macro.eps_y       = data(7,ib);         %y-displacement of the lateral boundaries
macro.sigma_x_top = data(8,ib);         %tangential stress integrated on the top boundary
macro.sigma_x_mid = data(9,ib);         %tangential stress integrated on the mid boundary
macro.sigma_x_bot = data(10,ib);        %tangential stress integrated on the bottom boundary
macro.E_brit      = data(11,ib);        %Energy dissipated in brittle failure
macro.area        = data(12,ib);        %Cumulated area of each avalanche
element.first_dam = find(data(2,ib), 1, 'first');
f                 = 1/(macro.dt);
data = [];

%% Remove a trend if there is one (not necessary, no impact on results)

dam_rate_fit = fit(macro.time', macro.dam_rate', 'poly1');
E_brit_fit = fit(macro.time', macro.E_brit', 'poly1');

%macro.dam_rate = macro.dam_rate - (dam_rate_fit.p1*macro.time + dam_rate_fit.p2);
%macro.E_brit = macro.E_brit - (E_brit_fit.p1*macro.time + E_brit_fit.p2);


%% Cut t-series to a max length
macro.dam_rate(max_cutoff:end) = [];
macro.E_brit(max_cutoff:end) = [];
macro.time(max_cutoff:end) = [];


%% Cut t-series after the first macro-rupture (first few cycles)
[pks, locs] = findpeaks(macro.sigma_x_top);
macro.dam_rate(1:locs(min_cutoff)) = [];
macro.E_brit(1:locs(min_cutoff)) = [];
macro.time(1:locs(min_cutoff)) = [];


%% CALCULATE PSD of the damage rate signal

%Fourier transform 

if (mod(length(macro.dam_rate),2) ~= 0)
macro.dam_rate(end) = [];
end
 
DFT_ = fft(macro.dam_rate(1:end));
DFT = DFT_(1:length(macro.dam_rate(1:end))/2+1);
rn2 = 1/(f*length(macro.dam_rate(1:end))) * abs(DFT).^2;
rn2(2:end-1) = 2*rn2(2:end-1);
s = 0:f/length(macro.dam_rate(1:end)):f/2;

PDS(1:length(rn2)) = rn2(1:length(rn2));


%% Calculating running means of the periodograms

rm_s(i,1:length(s)-2*floor(window/2)) = NaN;
rm_PDS_dam(i,1:length(s)-2*floor(window/2)) = NaN;


k = 0;
for j = 1+floor(window/2):1:length(s)-floor(window/2)
   k = k+1; 
    rm_s(i,k) = s(j); 
    rm_PDS_dam(i,k) = mean(PDS(j-floor(window/2):j+floor(window/2)));
   
end
s = []; PDS = [];

  
   

%% CALCULATE PSD of the Brittle energy release signal

%Fourier transform 

if (mod(length(macro.E_brit),2) ~= 0)
macro.E_brit(end) = [];
end
 
DFT_ = fft(macro.E_brit(1:end));
DFT = DFT_(1:length(macro.E_brit(1:end))/2+1);
rn2 = 1/(f*length(macro.E_brit(1:end))) * abs(DFT).^2;
rn2(2:end-1) = 2*rn2(2:end-1);
s = 0:f/length(macro.E_brit(1:end)):f/2;

PDS(1:length(rn2)) = rn2(1:length(rn2));


%% Calculating running means of the periodograms

rm_s(i,1:length(s)-2*floor(window/2)) = NaN;
rm_PDS_Ebrit(i,1:length(s)-2*floor(window/2)) = NaN;


k = 0;
for j = 1+floor(window/2):1:length(s)-floor(window/2)
   k = k+1; 
    rm_s(i,k) = s(j);
   
    rm_PDS_Ebrit(i,k) = mean(PDS(j-floor(window/2):j+floor(window/2)));
   
end
s = []; PDS = [];     

end


mean_rm_PDS_dam(m,1:length(rm_PDS_dam)) = sum(rm_PDS_dam, 1)/length(filenb);
mean_rm_PDS_Ebrit(m,1:length(rm_PDS_Ebrit)) = sum(rm_PDS_Ebrit, 1)/length(filenb);

end


%% --------------------------------------------------------------------------
%Plots

    
figure1 = figure('Name','Discrete Fourier Transform of the damage rate');
    axes1 = axes('Parent',figure1, 'XGrid', 'on', 'YGrid', 'on', 'YScale', 'log', 'XScale', 'log', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ... %'YTick',[1e2 1e3 1e4 1e5 1e6 1e7],... %'XTick',[1e-12 1e-10 1e-8 1e-6],...
    'FontSize', 14); 
    %xlim([1e0 1e5]);    
    %ylim([1e-15 1e-3]);
    box('on');
    hold('all');
    xlabel(' \it f \rm','FontSize', 12);
    ylabel(' PSD, damage rate \rm','FontSize', 12);
    
    
figure2 = figure('Name','Discrete Fourier Transform of E_{brit}');
    axes2 = axes('Parent',figure2, 'XGrid', 'on', 'YGrid', 'on', 'YScale', 'log', 'XScale', 'log', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ...%'YTick',[1e-45 1e-44 1e-43 1e-42 1e-41 1e-40 1e-39], 'XTick',[1e3 1e5 1e7 1e9],...
    'FontSize', 14); 
    %xlim([1e4 1e10]); %We = 0.001
    %xlim([1e3 1e9]); %We = 0.1
    xlim([1e2 1e8]); %We = 10
    %ylim([1e-48 1e-43]); %We = 0.001    
    %ylim([1e-45 1e-39]); %We = 0.1
    ylim([1e-42 1e-36]); %We = 10
    box('on');
    hold('all');    
    xlabel(' \it f \rm','FontSize', 12);
    ylabel(' PSD, \it E_{brit} \rm','FontSize', 12);
    

    
for m = 1:1:length(th)    
   
plot(rm_s(m, 2:end), mean_rm_PDS_dam(m, 2:length(rm_PDS_dam)), 'LineWidth', 1.0, 'Color', cc(m,:), 'DisplayName', ['\it T_h\rm = ' num2str(th(m))], 'Parent', axes1);
hold on;

plot(rm_s(m, 2:end), mean_rm_PDS_Ebrit(m, 2:length(rm_PDS_Ebrit)), 'LineWidth', 1.0, 'Color', cc(m,:), 'DisplayName', ['\it T_h\rm = ' num2str(th(m))], 'Parent', axes2);
hold on;

end 


%lgd = legend(axes1,'show');
%set(lgd, 'Location', 'southwest', 'FontSize', 10);
hold on;
saveas(figure1, [path_, 'fft_damrate.jpg']);
saveas(figure1, [path_, 'fft_damrate.pdf']);

%lgd = legend(axes2,'show');
%set(lgd, 'Location', 'southwest', 'FontSize', 10);
hold on;
saveas(figure2, [path_, 'fft_Ebrit.jpg']);
saveas(figure2, [path_, 'fft_Ebrit.pdf']);

