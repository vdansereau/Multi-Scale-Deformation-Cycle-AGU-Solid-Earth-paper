function data = mesh_vtk2mat(filename, nb_fields)  %Version reading tensors
% Important: 
%path = '/home/vdansereau/Documents/Rheolef/outputs/intermittency/output_4054100000/eta7/';
%field = 'output_6-';
%nb_fields=3;
%filename = [path field num2str(0) '.vtk'];
 %chargement d'un fichier VTK dans matlab, sous formes des structures utilisees pour Dadyn

global element

fid=fopen(filename);                             %Open file and skip first 4 lines
fgetl(fid); 
fgetl(fid); 
fgetl(fid); 
fgetl(fid);                                     

element.Nn=fscanf(fid, ['POINTS %d float']);    %Read nb of NODES

for i=1:element.Nn
	element.node_x(i)=fscanf(fid, '%f',1);      %x-coordinate of the nodes 
	element.node_y(i)=fscanf(fid, '%f',1);      %y-coordinate of the nodes
	element.node_z(i)=fscanf(fid, '%f',1);      %z-dimension : not used
end



fgetl(fid); 
temp=fscanf(fid, ['CELLS %d %d']); element.Ne=temp(1); %Read nb of ELEMENTS

for i=1:element.Ne
	if(fscanf(fid, '%d',1)~=3)
		disp('erreur lecture VTK: element non triangulaire?');
		return;
	else
		%on rajoute 1 a l'indice: vtk = regles du C : (0:n) alors que matlab: (1:n+1)
		element.num_node(i,1)=1+fscanf(fid, '%d',1);  %No corresponding to the three corners of a given element
		element.num_node(i,2)=1+fscanf(fid, '%d',1);
		element.num_node(i,3)=1+fscanf(fid, '%d',1);
	end
end

%convertion des numero de noeuds a la convention Dadyn
% element.num_node=sort(element.num_node(:,1:3),2);

%CELL TYPE (?)
fgetl(fid); 
temp=fscanf(fid, ['CELL_TYPES %d']); %Ne=temp(1);
for i=1:element.Ne
	    fscanf(fid, '%d',1);
end


%CELL DATA 

%dam
j = 1;  
fgetl(fid); 
temp=fscanf(fid, ['CELL_DATA %d']); data(1:nb_fields,1:element.Ne) = NaN;
fgetl(fid); 
fgetl(fid); 
fgetl(fid); 
    for i=1:element.Ne
        data(j,i) = fscanf(fid, '%f',1);
    end
%We
j = 2;            
fgetl(fid);  
fgetl(fid);
fgetl(fid);
fgetl(fid);
    for i=1:element.Ne
        data(j,i) = fscanf(fid, '%f',1);
    end
%sigma_h    
j = 3;
fgetl(fid);  
fgetl(fid);
fgetl(fid);  
fgetl(fid);
    for i=1:element.Ne
        data(j,i) = fscanf(fid, '%f',1);
    end    
    
        
%Calcul centre des elements
element.elem_x=mean(element.node_x(element.num_node(:,1:3)),2); %x-coordinate of element
element.elem_y=mean(element.node_y(element.num_node(:,1:3)),2); %y-coordinate of element



fclose(fid);
return


