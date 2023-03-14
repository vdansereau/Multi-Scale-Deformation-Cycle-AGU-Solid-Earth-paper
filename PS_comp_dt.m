%%This script plots the damage rate and brittle energy released power PSD
%%from the uniaxial compression experiments


clear all;
path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_10/comp_dt/alpha_4/th_10_8/';

%dt = [10^4 10^5 10^6 10^7 10^8];   % We = 0.001: Delta t = 10^4 s, 10^5 s, 10^6 s, 10^7 s, 10^8 s
%dt = [10^3 10^4 10^5 10^6 10^7];    % We = 0.1  : Delta t = 10^3 s, 10^4 s, 10^5 s, 10^6 s, 10^7 s
dt = [10^2 10^3 10^4 10^5 10^6];   % We = 10   : Delta t = 10^2 s, 10^3 s, 10^4 s, 10^5 s, 10^6 s

filenb = [1 2 3 4 5];
cc=parula(length(filenb));


figure1 = figure('Name','Discrete Fourier Transform of the damage rate');
    axes1 = axes('Parent',figure1, 'XGrid', 'on', 'YGrid', 'on', 'YScale', 'log', 'XScale', 'log', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ... %'YTick',[1e2 1e3 1e4 1e5 1e6 1e7],... %'XTick',[1e-12 1e-10 1e-8 1e-6],...
    'FontSize', 14); 
    xlim([1e2 1e12]);    
    ylim([1e-35 1e-15]);
    box('on');
    hold('all');
    xlabel(' \it f \rm','FontSize', 14);
    ylabel(' PSD, \it \Sigma \rm','FontSize', 14);
    
    
figure2 = figure('Name','Discrete Fourier Transform of E_{brit}');
    axes2 = axes('Parent',figure2, 'XGrid', 'on', 'YGrid', 'on', 'YScale', 'log', 'XScale', 'log', ...
    'XMinorTick', 'on', 'YMinorTick', 'on', ... %'YTick',[1e-25 1e-20 1e-15 1e-10],... %'XTick',[1e-12 1e-10 1e-8 1e-6],...
    'FontSize', 12); 
    %xlim([1e0 1e5]);    
    %ylim([1e-25 1e-15]);
    box('on');
    hold('all');    
    xlabel(' \it f \rm','FontSize', 12);
    ylabel(' PSD, E_{brit} \rm','FontSize', 12);
    
    
for i = length(filenb) : -1 : 1
data = load([path, 'test_', int2str(filenb(i)), '.txt'])';

global macro element

[macro.time,ia,ic]= unique(data(1,:)); %model iteration
macro.time        = macro.time - data(1,1);
macro.dt          = macro.time(2)-macro.time(1); %model time step (s)
macro.dam_rate    = data(8,ia);        %avalanche count
macro.k           = data(3,ia);        %nb of sub-iterations
macro.eps_x_top   = data(4,ia);        %x-displacement of the top boundary
macro.eps_x_mid   = data(5,ia);        %x-displacement of the mid boundary
macro.eps_x_bot   = data(6,ia);        %x-displacement of the bottom boundary
macro.eps_y       = data(7,ia);        %y-displacement of the lateral boundaries
macro.sigma_x_top = data(8,ia);        %tangential stress integrated on the top boundary
macro.sigma_x_mid = data(9,ia);        %tangential stress integrated on the mid boundary
macro.sigma_x_bot = data(10,ia);       %tangential stress integrated on the bottom boundary
macro.E_brit      = data(11,ia);       %Energy dissipated in brittle failure
macro.friction    = data(12,ia);       %Nb of frictional elements
element.first_dam = find(data(2,ia), 1, 'first');
f                 = 1/(macro.dt);
data = [];



%% Cut t-series after the first macro-rupture
[pks, locs] = findpeaks(macro.sigma_x_top);
macro.dam_rate(1:locs(10)) = [];
macro.E_brit(1:locs(10)) = [];


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

window = 5; %frequency-averaging window (odd number)

rm_s(i,1:length(s)-2*floor(window/2)) = NaN;
rm_PDS_dam(i,1:length(s)-2*floor(window/2)) = NaN;


k = 0;
for j = 1+floor(window/2):1:length(s)-floor(window/2)
   k = k+1; 
    rm_s(i,k) = s(j); 
    rm_PDS_dam(i,k) = mean(PDS(j-floor(window/2):j+floor(window/2)));
   
end
%s = []; PDS = [];

   
    plot(rm_s(i,2:end), rm_PDS_dam(i,2:length(rm_s(i,:))), 'LineWidth', 1.0, 'Color', cc(i,:), 'DisplayName', ['\Delta\it t \rm = ', num2str(dt(i))], 'Parent', axes1);
    hold on;
    



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

window = 5; %frequency-averaging window (odd number)

rm_s(i,1:length(s)-2*floor(window/2)) = NaN;
rm_PDS_Ebrit(i,1:length(s)-2*floor(window/2)) = NaN;


k = 0;
for j = 1+floor(window/2):1:length(s)-floor(window/2)
   k = k+1; 
    rm_s(i,k) = s(j);
   
    rm_PDS_Ebrit(i,k) = mean(PDS(j-floor(window/2):j+floor(window/2)));
   
end
s = []; PDS = [];
    

    plot(rm_s(i,2:end), rm_PDS_Ebrit(i,2:length(rm_s(i,:))), 'LineWidth', 1.0, 'Color', cc(i,:), 'DisplayName', ['\Delta\it t \rm = ', num2str(dt(i))], 'Parent', axes2);
    hold on;

    

end


lgd = legend(axes1,'show');
set(lgd, 'Location', 'northeast', 'FontSize', 10);
hold on;
saveas(figure1, [path , '/fft_sigma.jpg']);
saveas(figure1, [path , '/fft_sigma.pdf']);

lgd = legend(axes2,'show');
set(lgd, 'Location', 'northeast', 'FontSize', 10);
hold on;
saveas(figure2, [path , '/fft_Ebrit.jpg']);
saveas(figure2, [path , '/fft_Ebrit.pdf']);

