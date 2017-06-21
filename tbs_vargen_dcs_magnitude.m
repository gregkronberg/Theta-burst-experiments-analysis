%% variable: DCS intensity
%==========================================================================
% Take slopes and population spikes of baseline respones and during
% induction.  Data are stored as a structure called slopes in the Matlab
% Variables folder.  slopes is organized according to the slices reference
% structure, where each entry contains the slopes and info for each slice
%==========================================================================

clear all
close all
clc

%% file paths and global variables
%==========================================================================
current_path = mfilename('fullpath');
if current_path(1)=='D'
    % desktop paths
    fpath_raw = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
    fpath_processed = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
    fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
    fpath_analysis = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
    fpath_filters = 'D:\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters
    fpath_processed_images = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Images\'; % processed
else
    % laptop
    fpath_raw = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
    fpath_processed = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
    fpath_variables = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
    fpath_analysis = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
    fpath_filters = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters
    fpath_processed_images = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Images\'; % processed
end

%% load global structures
%==========================================================================
load(strcat(fpath_variables,'slices'));
load(strcat(fpath_variables,'electrode_location'));
load(strcat(fpath_variables,'dcs_magnitude'));

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
        dcs_magnitude{a,b,c,d,e}(f).slices = [];
        if sum(strcmp({dir_processed(:).name}',slices{a,b,c,d,e}(f).name))==1
            if sum(strcmp({dcs_magnitude{a,b,c,d,e}(:).name}',slices{a,b,c,d,e}(f).name))==0
%====================================== loop over individual slices

% load processed data
load(strcat(fpath_processed,slices{a,b,c,d,e}(f).name),...
    'pulset','indS','indD','baseIndex','fs','tpre','apical','basal',...
    'perforant','soma')

% dcs onset time in samples for TBS protocol
dcs_on = 1*fs;

% distance between recording electrodes
dist = .001*electrode_location{a,b,c,d,e}(f).dist;%(mm)

% time window to average before DCS onset
t_avg = 0.1*fs;

% rise time of dcs step
t_rise = .001*fs;

% average voltage before dcs
pre_soma = mean(indS(dcs_on-t_avg:dcs_on));
pre_dend = mean(indD(dcs_on-t_avg:dcs_on));

% average voltage after dcs
post_soma = mean(indS(dcs_on+t_rise:dcs_on+t_rise+t_avg));
post_dend = mean(indD(dcs_on+t_rise:dcs_on+t_rise+t_avg));

% difference between before and after (mV)
diff_soma = 1000*(post_soma-pre_soma);
diff_dend = 1000*(post_dend-pre_dend);

% voltage difference between two channels (positive = anodal, negative = cathodal) 
if basal~=0
    voltage_diff = diff_soma - diff_dend;
else
    voltage_diff = diff_dend - diff_soma;
end

field = voltage_diff/dist;

% store electrode_location variable
%==================================================================
dcs_magnitude{a,b,c,d,e}(f).dcs_magnitude = field;
dcs_magnitude{a,b,c,d,e}(f).slices = slices{a,b,c,d,e}(f);
dcs_magnitude{a,b,c,d,e}(f).name = slices{a,b,c,d,e}(f).name;
                
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

save(strcat(fpath_variables,'dcs_magnitude.mat'),'dcs_magnitude')


