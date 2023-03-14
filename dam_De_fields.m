%%This script plots snapshots of the damage field
%%and the associated distribution of the De number
$$for the case of De = 0.001
$$Used to produce figure 5: De_dam_fields


clear all;
path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_001/dt_10_5/C_10_4/th_10_10/alpha_4/ddam_10/';
filenb = [0];


%% Set the Deborah number of the simulation
%De = De0*(1-d)^(alpha-1);	%Effective Deborah number

De0 = 0.001			%Mean Deborah number of the domain: change for another set of simulations
De0_continental = 1/1.5*De0;  	%Deborah number of the continental crust
De0_oceanic = 1/0.5*De0;		%Deborah number of the oceanic crust


%% Set the frequency of field outputs, i.e., every output_frequency t-step after first dam event
output_frequency = 10; 		


%% Load data
data = load([path, 'test_', int2str(filenb), '.txt'])';

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
element.first_dam = find(data(3,ia), 1, 'first');
macro.dt          = macro.time(2)-macro.time(1); %model time step (s)
f                 = 1/(macro.dt);
data = [];



%% Plot a figure of the macro stress and macro damage increment vs time,
%  to identify good moments at which to plot the snapshots in damage and
%  distribution of De

figure1 = figure;
[AX, H1, H2] = plotyy(macro.time, macro.sigma_x_top, macro.time, macro.dam_rate);
set(H1, 'Color', 'k', 'LineWidth', 1.0, 'Displayname', ['disp_{top}']);
set(H2, 'Color', 'c', 'LineStyle', '-', 'LineWidth', 1.0, 'HandleVisibility', 'off');
set(AX(1), 'FontSize', 10, 'xlim', [0 0.5e-5], 'PlotBoxAspectRatio',[3 1 1]); 
set(AX(2), 'FontSize', 10, 'xlim', [0 0.5e-5], 'ylim', [0 1e-3], 'ytick', [0 1e-3], 'PlotBoxAspectRatio', [3 1 1]);

set(get(AX(1),'Ylabel'),'String','macro stress', 'Color', [0 0 0], 'FontSize', 10);  set(AX(1), 'YColor', [0 0 0]);
set(get(AX(2),'Ylabel'),'String','macro damage incr.', 'Color', 'c', 'FontSize', 10);       set(AX(2), 'YColor', 'c');
hold on;

xlabel('time', 'FontSize', 10);

saveas(figure1, [path, 'stress_dam.pdf']);
hold off;



%% Damage field files
%Check field file
nb_fields = 3;
field = ['output_' num2str(filenb) '-']; name_size = size(field,2);   %variable name
allfiles=dir([path field '*.vtk']);               %all vtk files for variable
 
%sorting the vtk files 
for i = 1 : size(allfiles,1)
     fname = allfiles(i).name;
     itrs(i) = str2num(fname(name_size+1:end-4)); %vtk file nb
end
itrs = sort(itrs);
 
%Load a field file to get the mesh grid information: store into element.
firstfile = mesh_vtk2_tensor([path field num2str(itrs(1)) '.vtk'], nb_fields);
 
 

%% Set colorbar
%n=255;
%gris=hot(n);
%rev_gris((n:-1:1),:)=gris(-(-n:-1),:);
%tmp=rev_gris(:,3); rev_gris(:,3)=rev_gris(:,1); rev_gris(:,1)=tmp;
%rev_gris(end,1)=0.99; rev_gris(end,2)=1.0; rev_gris(end,3)=1.0;
    
    
%% Plot damage and De field at at given (k) time for a stress minimum      

%Identify the time step corresponding to a stress minimum
time_step = 2.213e-6; 

%Identify the closest output file
k = floor(((time_step - element.first_dam*macro.dt)/macro.dt)/output_frequency);

data = mesh_vtk2mat([path field num2str(itrs(k)) '.vtk'], 3); 


%---------------------------------------------------------------------------
%% Figure of the damage field

    dam(1:element.Ne) = data(1,1:element.Ne);
 
    figure2 = figure('Color', [1 1 1]);
    axis off;
    pbaspect([2 1 1]);
    box on; 
    colormap([0.800000011920929 0.800000011920929 0.800000011920929;0.808695673942566 0.808695673942566 0.765217423439026;0.817391335964203 0.817391335964203 0.730434775352478;0.826086938381195 0.826086938381195 0.695652186870575;0.834782600402832 0.834782600402832 0.660869598388672;0.843478262424469 0.843478262424469 0.626086950302124;0.852173924446106 0.852173924446106 0.591304361820221;0.860869586467743 0.860869586467743 0.556521773338318;0.86956524848938 0.86956524848938 0.52173912525177;0.878260850906372 0.878260850906372 0.486956536769867;0.886956512928009 0.886956512928009 0.452173918485641;0.895652174949646 0.895652174949646 0.417391300201416;0.904347836971283 0.904347836971283 0.382608711719513;0.91304349899292 0.91304349899292 0.347826093435287;0.921739161014557 0.921739161014557 0.313043475151062;0.930434763431549 0.930434763431549 0.278260886669159;0.939130425453186 0.939130425453186 0.243478268384933;0.947826087474823 0.947826087474823 0.208695650100708;0.95652174949646 0.95652174949646 0.173913046717644;0.965217411518097 0.965217411518097 0.139130443334579;0.973913073539734 0.973913073539734 0.104347825050354;0.982608675956726 0.982608675956726 0.0695652216672897;0.991304337978363 0.991304337978363 0.0347826108336449;1 1 0;1 0.95652174949646 0;1 0.91304349899292 0;1 0.869565188884735 0;1 0.826086938381195 0;1 0.782608687877655 0;1 0.739130437374115 0;1 0.695652186870575 0;1 0.652173936367035 0;1 0.60869562625885 0;1 0.56521737575531 0;1 0.52173912525177 0;1 0.47826087474823 0;1 0.434782594442368 0;1 0.391304343938828 0;1 0.347826093435287 0;1 0.304347813129425 0;1 0.260869562625885 0;1 0.217391297221184 0;1 0.173913046717644 0;1 0.130434781312943 0;1 0.0869565233588219 0;1 0.0434782616794109 0;1 0 0;0.941176474094391 0 0;0.882352948188782 0 0;0.823529422283173 0 0;0.764705896377563 0 0;0.705882370471954 0 0;0.647058844566345 0 0;0.588235318660736 0 0;0.529411792755127 0 0;0.470588237047195 0 0;0.411764711141586 0 0;0.352941185235977 0 0;0.294117659330368 0 0;0.235294118523598 0 0;0.176470592617989 0 0;0.117647059261799 0 0;0.0588235296308994 0 0;0 0 0]);
    hold all;
    
    h2 = trisurf(element.num_node, element.node_x, element.node_y, zeros(size(element.node_x)), dam); view(2);
    set(h2, 'EdgeColor','none');
    set(gca, 'CLim', [0.0 1.0]);
    colorbar('EastOutside', 'Fontsize', 14);

    %title('Level of damage', 'FontSize', 12);
    saveas(figure2, [path, 'dam_field.pdf']);
    saveas(figure2, [path, 'dam_field.jpg']);
    hold off;
    
%---------------------------------------------------------------------------    
%% Figure of the De field

    De(1:element.Ne) = data(2,1:element.Ne);
    %Setting the De value of undammaged elements to NaN;
    De(find(De >= De0_continental)) = NaN;
    De(find(De == De0_oceanic)) = NaN;
 
    figure3 = figure();
    axis off;
    pbaspect([2 1 1]);
    box on; 
    colormap([0.800000011920929 0.800000011920929 0.800000011920929;0.808695673942566 0.808695673942566 0.765217423439026;0.817391335964203 0.817391335964203 0.730434775352478;0.826086938381195 0.826086938381195 0.695652186870575;0.834782600402832 0.834782600402832 0.660869598388672;0.843478262424469 0.843478262424469 0.626086950302124;0.852173924446106 0.852173924446106 0.591304361820221;0.860869586467743 0.860869586467743 0.556521773338318;0.86956524848938 0.86956524848938 0.52173912525177;0.878260850906372 0.878260850906372 0.486956536769867;0.886956512928009 0.886956512928009 0.452173918485641;0.895652174949646 0.895652174949646 0.417391300201416;0.904347836971283 0.904347836971283 0.382608711719513;0.91304349899292 0.91304349899292 0.347826093435287;0.921739161014557 0.921739161014557 0.313043475151062;0.930434763431549 0.930434763431549 0.278260886669159;0.939130425453186 0.939130425453186 0.243478268384933;0.947826087474823 0.947826087474823 0.208695650100708;0.95652174949646 0.95652174949646 0.173913046717644;0.965217411518097 0.965217411518097 0.139130443334579;0.973913073539734 0.973913073539734 0.104347825050354;0.982608675956726 0.982608675956726 0.0695652216672897;0.991304337978363 0.991304337978363 0.0347826108336449;1 1 0;1 0.95652174949646 0;1 0.91304349899292 0;1 0.869565188884735 0;1 0.826086938381195 0;1 0.782608687877655 0;1 0.739130437374115 0;1 0.695652186870575 0;1 0.652173936367035 0;1 0.60869562625885 0;1 0.56521737575531 0;1 0.52173912525177 0;1 0.47826087474823 0;1 0.434782594442368 0;1 0.391304343938828 0;1 0.347826093435287 0;1 0.304347813129425 0;1 0.260869562625885 0;1 0.217391297221184 0;1 0.173913046717644 0;1 0.130434781312943 0;1 0.0869565233588219 0;1 0.0434782616794109 0;1 0 0;0.941176474094391 0 0;0.882352948188782 0 0;0.823529422283173 0 0;0.764705896377563 0 0;0.705882370471954 0 0;0.647058844566345 0 0;0.588235318660736 0 0;0.529411792755127 0 0;0.470588237047195 0 0;0.411764711141586 0 0;0.352941185235977 0 0;0.294117659330368 0 0;0.235294118523598 0 0;0.176470592617989 0 0;0.117647059261799 0 0;0.0588235296308994 0 0;0 0 0]);
    hold all;
    
    h2 = trisurf(element.num_node, element.node_x, element.node_y, zeros(size(element.node_x)), log10(De)); view(2);
    set(h2, 'EdgeColor','none');
    set(gca, 'CLim', [-13 -3]);
    colorbar('EastOutside', 'Fontsize', 14);
    
    %title('De', 'FontSize', 12);
    saveas(figure3, [path, 'De_field.pdf']);
    hold off;
    

%--------------------------------------------------------------------------
%% Plot distribution of De over a window centred on a time k: around a stress minimum and a stress maximum
% I tried a window of 5 time steps but it does not make a difference between 1 and 5 time steps

%--------------------------------------------------------------------------
%Identify the time step corresponding to a stress minimum
time_step = 2.213e-6; 
%Identify the closest output file
k = floor(((time_step - element.first_dam*macro.dt)/macro.dt)/output_frequency);

window = 0; 
for i = 1:1:window*2+1;
    data = mesh_vtk2mat([path field num2str(itrs((k-window)+(i-1))) '.vtk'], 3); 
    De_min(element.Ne*(i-1)+1:element.Ne*i) = data(2,1:element.Ne);
end

    De_min(find(De_min >= De0_continental)) = NaN;
    De_min(find(De_min == De0_oceanic)) = NaN;
    

%--------------------------------------------------------------------------
%%Identify the time step corresponding to a stress maximum
%time_step = 3.77e-6; 
%Identify the closest output file
%k = floor(((time_step - element.first_dam*macro.dt)/macro.dt)/output_frequency);
    
%window = 0;
%for i = 1:1:window*2+1;
%    data = mesh_vtk2mat([path field num2str(itrs((k-window)+(i-1))) '.vtk'], 3); 
%    De_max(element.Ne*(i-1)+1:element.Ne*i) = data(2,1:element.Ne);
%end

%    De_max(find(De_max >= De0_continental)) = NaN;
%    De_max(find(De_max == De0_oceanic)) = NaN;

    
    
    
    figure4 = figure;
    axes4 = axes('Parent', figure4, 'FontSize', 8);
    pbaspect([1 2 1]);
    edges = [-13 : 1 : -3];
    set(axes4, 'XScale', 'log', 'Xlim', [10^(-13) 10^(-3)], 'xtick', [10^(-12) 10^(-10) 10^(-8) 10^(-6) 10^(-4)], 'XMinorTick','off', 'Ylim', [0 1]);
    ylabel('Probability', 'FontSize', 8); 
    xlabel('De', 'FontSize', 8);
    hold('all');
    
    
    
    histogram(De_min, 10.^edges, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'EdgeColor', 'k', 'FaceColor', 'k', 'LineWidth', 1.0, 'Displayname', [''], 'Parent', axes4);
    hold on;
    
%    %histogram(De_max, 10.^edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', 'c', 'LineWidth', 1.5, 'Displayname', [''], %'Parent', axes4);
%    %hold on;

    saveas(figure4, [path, 'De_hist.pdf']);