%% analyze distance and orientation of electrodes in slice

clear all; close all; clc

%% file paths
fpathR = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Images\';
fpathP = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Images\';
fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\Processed Variables\'
% fpathR = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Slice Images\';
% fpathP = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Images\';

%% directory listings of all files
directR = dir(strcat(fpathR,'*.png*')); % raw matlab files
directP = dir(strcat(fpathP,'*.mat*')); % processed matlab files

%% cell arrays of files names
slicesR = cell(size(directR,1),1); %create empty cell array
for a = 1:size(directR,1);
    slicesR{a} = directR(a).name; % raw file names
end
slicesP = cell(size(directP,1),1); %create empty cell array
for a = 1:size(directP,1);
    
    slicesP{a} = directP(a).name; % processed file names
end

%% process new slices
if length(slicesR)~= length(slicesP);                                       % are there new slices?
    for i = 1:length(slicesR);                                              % loop over all slices
        if sum(strcmp(slicesP,strcat(slicesR{i}(1:end-4),'.mat')))==0                              % check if current slice is new
            im = imread(strcat(fpathR,slicesR{i}));
            [rows,cols,channels] = size(im);
            fig = imshow(im);hold on
            name = slicesR{i}(1:end-4);
            
            % coordinates of points along first field wire
            title('Select two points along first field wire')
            [x,y] = ginput(2);
            field_slope_1 = (y(2)-y(1))/(x(2)-x(1));
            field_int_1 = y(2)-field_slope_1*x(2);
            field_x_1 = 1:cols;
            field_y_1 = field_slope_1*field_x_1 + field_int_1;
            plot(field_x_1,field_y_1)
            
            % coordinates of points along second field wire
            title('Select two points along second field wire')
            [x,y] = ginput(2);
            field_slope_2 = (y(2)-y(1))/(x(2)-x(1));
            field_int_2 = y(2)-field_slope_2*x(2);
            field_x_2 = 1:cols;
            field_y_2 = field_slope_2*field_x_2 + field_int_2;
            plot(field_x_2,field_y_2)
            
            % soma coordinate
            title('Select somatic layer location')
            [soma_x,soma_y] = ginput(1);
            plot(soma_x,soma_y,'r.','MarkerSize',10)
            
            % stimulating electrode coordinate
            title('Select stim electrode location')
            [stim_x,stim_y] = ginput(1);
            plot(stim_x,stim_y,'r.','MarkerSize',10)
            
            % dendritic recording electrode
            title('Select dendritic electrode location')
            [rec_dend_x,rec_dend_y] = ginput(1);
            plot(rec_dend_x,rec_dend_y,'r.','MarkerSize',10)
            
             % somatic recording electrode
            title('Select somatic electrode location')
            [rec_soma_x,rec_soma_y] = ginput(1);
            plot(rec_soma_x,rec_soma_y,'r.','MarkerSize',10)
            
            % select points in consecutive mesh boxes (i.e. 400 um apart)
            title('Select mesh box edges')
            [mesh_x,mesh_y] = ginput(2);
            mesh = [mesh_x,mesh_y];
            dist = 2.5*norm(diff(mesh,1)); % measurement standard 1um
            
            % distance from stim to soma
            stim_to_soma = norm([stim_x,stim_y]-[soma_x,soma_y])/dist;
            
            % distance from stim to rec dend
            stim_to_dend = norm([stim_x,stim_y]-[rec_dend_x,rec_dend_y])/dist;
            
            % distance from rec to rec
            rec_to_rec = norm([rec_dend_x,rec_dend_y]-[rec_soma_x,rec_soma_y])/dist;
            
            save(strcat(fpathP,name,'.mat'))
        end
    end
end