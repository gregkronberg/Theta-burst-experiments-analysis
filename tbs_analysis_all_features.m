%% analysis: all features
%==========================================================================
% plot slopes over time for each condtion
%==========================================================================

clear all
close all
clc

%% file paths and global variables
%==========================================================================
current_path = pwd;
if current_path(1)=='D'
    % desktop paths
    fpath_raw = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
    fpath_processed = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
    fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
    fpath_analysis = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
    fpath_filters = 'D:\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters
else
    % laptop
    fpath_raw = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
    fpath_processed = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
    fpath_variables = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
    fpath_analysis = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
    fpath_filters = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters
    fpath_processed_images = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Images\'; % processed
end

% load global structures
load(strcat(fpath_variables,'slices.mat')); 
load(strcat(fpath_variables,'slopes.mat')); 
load(strcat(fpath_variables,'dcs_magnitude.mat')); 
load(strcat(fpath_variables,'drift.mat')); 
load(strcat(fpath_variables,'electrode_location.mat')); 
load(strcat(fpath_variables,'fiber_volley.mat')); 
load(strcat(fpath_variables,'soma_maxslope.mat'));
load(strcat(fpath_variables,'fepsp_area.mat'));

%% create list of all slices
%==========================================================================
cnt = 0
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        for f = 1:length(slices{a,b,c,d,e})
cnt = cnt+1;
ind_block = slices{a,b,c,d,e}(f).indBlock(1);
% slices
features(cnt).name = slices{a,b,c,d,e}(f).name;
features(cnt).date = slices{a,b,c,d,e}(f).date;
features(cnt).slice_num = slices{a,b,c,d,e}(f).slice_num;
features(cnt).ind_time = slices{a,b,c,d,e}(f).blocktimes(ind_block);
features(cnt).age = slices{a,b,c,d,e}(f).age;
features(cnt).hemi = slices{a,b,c,d,e}(f).hemi;
features(cnt).height = slices{a,b,c,d,e}(f).height;

% slopes
features(cnt).slopes_end = mean(slopes{a,b,c,d,e}(f).slopes_norm(end-9:end));
features(cnt).slopes_base_mean = slopes{a,b,c,d,e}(f).slopes_base_mean;
features(cnt).spikes_base_mean = slopes{a,b,c,d,e}(f).spikes_base_mean;
features(cnt).spikes_slopes_ratio_base = slopes{a,b,c,d,e}(f).spikes_slopes_ratio_base;

% drift
features(cnt).drift = drift{a,b,c,d,e}(f).slopes_drift;
features(cnt).drift_time = drift{a,b,c,d,e}(f).time;

% dcs magnitude
if length(dcs_magnitude{a,b,c,d,e})>1 && isempty(dcs_magnitude{a,b,c,d,e}(1).slices)==0
    features(cnt).dcs_magnitude = dcs_magnitude{a,b,c,d,e}(f).dcs_magnitude;
    % electrode loacation
    features(cnt).stim_to_soma = electrode_location{a,b,c,d,e}(f).stim_to_soma;
    features(cnt).stim_to_dend = electrode_location{a,b,c,d,e}(f).stim_to_dend;
else
    features(cnt).dcs_magnitude = [];
    features(cnt).stim_to_soma = [];
    features(cnt).stim_to_dend = [];
end

% fiber volley
features(cnt).fiber_volley = mean(fiber_volley{a,b,c,d,e}(f).fiber_volley_norm(1,:));

% burst area
features(cnt).burst_area_first = fepsp_area{a,b,c,d,e}(f).burst_area_norm(1);
features(cnt).burst_area_mean = fepsp_area{a,b,c,d,e}(f).burst_area_mean;
features(cnt).burst_area_max = fepsp_area{a,b,c,d,e}(f).burst_area_max;
features(cnt).burst_area_adapt = fepsp_area{a,b,c,d,e}(f).burst_area_adapt;

% soma max slope
features(cnt).soma_maxslope_first = soma_maxslope{a,b,c,d,e}(f).maxslope_norm(1);
features(cnt).soma_maxslope_mean = soma_maxslope{a,b,c,d,e}(f).maxslope_mean;
features(cnt).soma_maxslope_max = soma_maxslope{a,b,c,d,e}(f).maxslope_max;
features(cnt).soma_maxslope_adapt = soma_maxslope{a,b,c,d,e}(f).maxslope_adapt;
features(cnt).soma_maxslope_i_first = soma_maxslope{a,b,c,d,e}(f).maxslope_i_norm(1);
features(cnt).soma_maxslope_i_mean = soma_maxslope{a,b,c,d,e}(f).maxslope_i_mean;
features(cnt).soma_maxslope_i_max = soma_maxslope{a,b,c,d,e}(f).maxslope_i_max;
features(cnt).soma_maxslope_i_adapt = soma_maxslope{a,b,c,d,e}(f).maxslope_i_adapt;

                        end
                    end
                end
            end
        end
    end
end

% store all data as single matrix

% canonical correlation                           
                            
                            


                        
