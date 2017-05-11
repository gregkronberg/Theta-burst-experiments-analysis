%% variable: add filtered time series
%==========================================================================
% notes
%==========================================================================

clear all
close all
clc

%% file paths
%==========================================================================
% desktop
fpath_raw = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
fpath_processed = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
fpath_analysis = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
fpath_filters = 'D:\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters

%% load global slices structure
%==========================================================================
load(strcat(fpath_variables,'slices'));

%% load filters
%==========================================================================
lowpass_1000 = load(strcat(fpath_filters,'iir_lowpass_1000Hz_fs10k.mat'));

%% process individual slices
%==========================================================================
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    
%===================================== loop over experimental conditions
if isempty(slices{a,b,c,d,e})==0
    for f = 1:length(slices{a,b,c,d,e})
%====================================== loop over individual slices

        % load processed data
        load(strcat(fpath_processed,slices{a,b,c,d,e}(f).name))
            
        %% apply filters
        %==================================================================
        % iir_bandpass_200Hz_600Hz_fs10k
        indD_filt_low_1000 = filtfilt(lowpass_1000.dfilt,indD);
        indS_filt_low_1000 = filtfilt(lowpass_1000.dfilt,indS);
        baseS_filt_low_1000 = filtfilt(lowpass_1000.dfilt,baseS);
        baseD_filt_low_1000 = filtfilt(lowpass_1000.dfilt,baseD);
        
        %% save to processed data folder
        %==================================================================
        % save all variables except slices
        save(strcat(fpath_processed,slices{a,b,c,d,e}(f).name), '-regexp', '^(?!(slices)$).') 
        
        %===================================== end loop over individual slices
    end
end
%===================================== end loop over experimental conditions
                end
            end
        end
    end
end