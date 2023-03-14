%%This script makes a movie of 
(i) the temporal evolution of the macro shear stress and macroscopic damage increments (lower panel)
(ii) the temporal evolution of the field of level of damage, d (log scale)
%It is used to make the animation presented in the Supplementary information. 


clear all;

%% Path to the output files
path = '/home/danserev/Documents/Rheolef/SEISMAZE/article1/faille/We_0_001/dt_10_5/C_10_4/th_10_10/alpha_4/ddam_10/';
filenb = [0];


%% Precise the temporal resolution of the damage field output files (set in the code)
dt_output = 10; %time-resolution of the field output


%% Load the output (text) files
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
element.first_dam = find(data(2,ia), 1, 'first');
macro.dt          = macro.time(2)-macro.time(1); %model time step (s)
f                 = 1/(macro.dt);
data = [];



%% Load the damage (field) files
%Check field file (distocrit)
nb_fields = 3;
field = 'output_0-'; name_size = size(field,2);   %variable name
allfiles=dir([path field '*.vtk']);               %all vtk files for variable
 
%sorting the vtk files (numbers get mixed in the transfer)
for i = 1 : size(allfiles,1);
     fname = allfiles(i).name;
     itrs(i) = str2num(fname(name_size+1:end-4)); %vtk file nb
end
itrs = sort(itrs);
 
%Load a field file to get the mesh grid information: store into element.
firstfile = mesh_vtk2_tensor([path field num2str(itrs(1)) '.vtk'], nb_fields);
 
 

%% Colorbar
    %    n=255;
    %    gris=hot(n);
    %    rev_gris((n:-1:1),:)=gris(-(-n:-1),:);
    %    tmp=rev_gris(:,3); rev_gris(:,3)=rev_gris(:,1); rev_gris(:,1)=tmp;
    %    rev_gris(end,1)=0.99; rev_gris(end,2)=1.0; rev_gris(end,3)=1.0;
    


%% Set the figure and figure axes
figure1 = figure('position',[100 100 850 600], 'Color', [1 1 1]);
%tstep = nn/2;     %t-step for movie frames
tmax = 5*10^-6;
tmin = 0;
K = floor((tmax-tmin)/(macro.dt*dt_output));        %nb of frames


for k = 1:1:K

    data = mesh_vtk2mat([path field num2str(itrs(k)) '.vtk'], 3); 
    dam(1:element.Ne) = data(1,1:element.Ne);
 
    subplot1 = subplot(2,2,[1 2], 'Parent', figure1);
    axis off;
    box(subplot1,'on');  
    %pbaspect([1 2 1]);
    hold(subplot1,'all');
    colormap([0.800000011920929 0.800000011920929 0.800000011920929;0.808695673942566 0.808695673942566 0.765217423439026;0.817391335964203 0.817391335964203 0.730434775352478;0.826086938381195 0.826086938381195 0.695652186870575;0.834782600402832 0.834782600402832 0.660869598388672;0.843478262424469 0.843478262424469 0.626086950302124;0.852173924446106 0.852173924446106 0.591304361820221;0.860869586467743 0.860869586467743 0.556521773338318;0.86956524848938 0.86956524848938 0.52173912525177;0.878260850906372 0.878260850906372 0.486956536769867;0.886956512928009 0.886956512928009 0.452173918485641;0.895652174949646 0.895652174949646 0.417391300201416;0.904347836971283 0.904347836971283 0.382608711719513;0.91304349899292 0.91304349899292 0.347826093435287;0.921739161014557 0.921739161014557 0.313043475151062;0.930434763431549 0.930434763431549 0.278260886669159;0.939130425453186 0.939130425453186 0.243478268384933;0.947826087474823 0.947826087474823 0.208695650100708;0.95652174949646 0.95652174949646 0.173913046717644;0.965217411518097 0.965217411518097 0.139130443334579;0.973913073539734 0.973913073539734 0.104347825050354;0.982608675956726 0.982608675956726 0.0695652216672897;0.991304337978363 0.991304337978363 0.0347826108336449;1 1 0;1 0.95652174949646 0;1 0.91304349899292 0;1 0.869565188884735 0;1 0.826086938381195 0;1 0.782608687877655 0;1 0.739130437374115 0;1 0.695652186870575 0;1 0.652173936367035 0;1 0.60869562625885 0;1 0.56521737575531 0;1 0.52173912525177 0;1 0.47826087474823 0;1 0.434782594442368 0;1 0.391304343938828 0;1 0.347826093435287 0;1 0.304347813129425 0;1 0.260869562625885 0;1 0.217391297221184 0;1 0.173913046717644 0;1 0.130434781312943 0;1 0.0869565233588219 0;1 0.0434782616794109 0;1 0 0;0.941176474094391 0 0;0.882352948188782 0 0;0.823529422283173 0 0;0.764705896377563 0 0;0.705882370471954 0 0;0.647058844566345 0 0;0.588235318660736 0 0;0.529411792755127 0 0;0.470588237047195 0 0;0.411764711141586 0 0;0.352941185235977 0 0;0.294117659330368 0 0;0.235294118523598 0 0;0.176470592617989 0 0;0.117647059261799 0 0;0.0588235296308994 0 0;0 0 0]);

        h2 = trisurf(element.num_node, element.node_x, element.node_y, zeros(size(element.node_x)), log(dam)); view(2);
        set(h2, 'EdgeColor','none', 'Parent', subplot1);
        set(gca, 'CLim', [-0.20 -0.01]);
        colorbar('peer',subplot1,'EastOutside', 'Fontsize', 14, 'Ticks', [-0.20 -0.15 -0.10 -0.05]);

        title('Level of damage, log_{10}(\it d \rm)', 'FontSize', 11);  

          
    subplot2 = subplot(2,2,[3 4], 'Parent', figure1);
    box(subplot2,'on');
    hold(subplot2,'all');
    colormap();
    

    [AX, H1, H2] = plotyy(  'Parent', subplot2, macro.time(1:k*dt_output), macro.sigma_x_top(1:k*dt_output), ...
                            'Parent', subplot2, macro.time(1:k*dt_output), macro.dam_rate(1:k*dt_output) );
    set(H1, 'Color', [0 0 0], 'DisplayName', '\it \sigma\ rm', 'LineWidth', 1.5);
    set(H2, 'Color', [0.6 0.6 0.6], 'LineStyle', '-', 'DisplayName', 'damage rate', 'LineWidth', 1.5);
    set(AX(1), 'FontSize', 14, 'xlim', [0 0.000005], 'ylim', [0 2e-7], 'YTick', [0 0.00000005 0.0000001 0.00000015 0.0000002]);
    set(AX(2), 'FontSize', 14, 'xlim', [0 0.000005], 'ylim', [0 1e-3], 'YTick', [0 0.001]);
    set(get(AX(1),'Ylabel'),'String','macro stress', 'Color', [0 0 0], 'FontSize', 14);      set(AX(1), 'YColor', [0 0 0]);
    set(get(AX(2),'Ylabel'),'String','macro damage incr.', 'Color', [0.6 0.6 0.6], 'FontSize', 14);     set(AX(2), 'YColor', [0.6 0.6 0.6]);

    xlabel('time', 'FontSize', 14);
  
    title('Macroscopic stress and damage increment evolution', 'FontSize', 11);

    M(k) = getframe(figure1);
    hold off;
    clf;
    
end


[h, w, p] = size(M(1).cdata);  % use 1st frame to get dimensions

v = VideoWriter([path, 'damage_movie2.avi']);
v.FrameRate = 60;    % Default 30
v.Quality = 100;    % Default 75
open(v);
writeVideo(v,M);
close(v);

