%% filter time series data and store with processed data files

clear all;
close all;
clc

%% file paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% desktop
fpath_processed = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
fpath_filters = 'D:\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters
% laptop
% fpathP = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';

%% list slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[temp,slices_filter] = folder_diff(fpath_processed,fpath_processed,'.mat');

%% load filters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bandpass_200_600 = load(strcat(fpath_filters,'iir_bandpass_200Hz_600Hz_fs10k.mat'));
highpass_1 = load(strcat(fpath_filters,'iir_highpass_1Hz_fs10k.mat'));

%%
for aa = 1:length(slices_filter)
    load(strcat(fpath_processed,slices_filter{aa}))
    
    %% apply filters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % iir_bandpass_200Hz_600Hz_fs10k
    indD_filt_band_200_600 = filtfilt(bandpass_200_600.dfilt,indD);
    indS_filt_band_200_600 = filtfilt(bandpass_200_600.dfilt,indS);
    baseS_filt_band_200_600 = filtfilt(bandpass_200_600.dfilt,baseS);
    baseD_filt_band_200_600 = filtfilt(bandpass_200_600.dfilt,baseD);
    
    % iir_highpass_1Hz_fs10k
    indD_filt_high_1 = filtfilt(highpass_1.dfilt,indD);
    indS_filt_high_1 = filtfilt(highpass_1.dfilt,indS);
    baseS_filt_high_1 = filtfilt(highpass_1.dfilt,baseS);
    baseD_filt_high_1 = filtfilt(highpass_1.dfilt,baseD);
    save(strcat(fpath_processed,slices_filter{aa}))
end
