%% variable: electrode locations
%==========================================================================
% notes
%==========================================================================

clear all
close all
clc

%% file paths
%==========================================================================
% desktop
% fpath_raw = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
% fpath_processed = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
% fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
% fpath_analysis = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
% fpath_filters = 'D:\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters
% fpath_processed_images = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Images\';

% laptop
fpath_raw = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
fpath_processed = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
fpath_variables = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
fpath_analysis = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
fpath_filters = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters
fpath_processed_images = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Images\'; % processed

%% load global slices structure
%==========================================================================
load(strcat(fpath_variables,'slices'));

%% load global electrode location structure
%==========================================================================
load(strcat(fpath_variables,'electrode_location'));

%% list processed images
%==========================================================================
dir_processed = dir(strcat(fpath_processed_images,'*.mat'));

%% process new slices
%==========================================================================
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    
%===================================== loop over experimental conditions
if isempty(slices{a,b,c,d,e})==0
    for f = 1:length(slices{a,b,c,d,e})
        % check for corresponding processed image
        if sum(strcmp({dir_processed(:).name}',slices{a,b,c,d,e}(f).name))==1
            % check if slice has already been stored in variable structure
            if sum(strcmp({electrode_location{a,b,c,d,e}(:).name}',slices{a,b,c,d,e}(f).name))==0
%====================================== loop over individual slices

% load processed data
%==========================================================
load(strcat(fpath_processed_images,slices{a,b,c,d,e}(f).name),...
    'dist','stim_to_soma','stim_to_dend','rec_to_rec',...
    'field_slope_1','field_int_1','field_slope_2',...
    'field_int_2')

% store electrode_location variable
%==================================================================
electrode_location{a,b,c,d,e}(f).dist = dist;
electrode_location{a,b,c,d,e}(f).stim_to_soma = stim_to_soma;
electrode_location{a,b,c,d,e}(f).stim_to_dend = stim_to_dend;
electrode_location{a,b,c,d,e}(f).rec_to_rec = rec_to_rec;
electrode_location{a,b,c,d,e}(f).field_wire1 = [field_slope_1 field_int_1];
electrode_location{a,b,c,d,e}(f).field_wire2 = [field_slope_2 field_int_2];
electrode_location{a,b,c,d,e}(f).slices = slices{a,b,c,d,e}(f);
electrode_location{a,b,c,d,e}(f).name = slices{a,b,c,d,e}(f).name;
                
            %===================================== end loop over individual slices
            end
        end
    end
end
%===================================== end loop over experimental conditions
                end
            end
        end
    end
end

save(strcat(fpath_variables,'electrode_location.mat'),'electrode_location')