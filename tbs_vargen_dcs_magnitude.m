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

%% file paths
%==========================================================================
% desktop
% fpath_raw = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
% fpath_processed = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
% fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
% fpath_analysis = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
% fpath_filters = 'D:\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters

% laptop
fpath_raw = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
fpath_processed = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
fpath_variables = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
fpath_analysis = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
fpath_filters = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters

%% load global slices structure
%==========================================================================
load(strcat(fpath_variables,'slices'));

%% load global slopes structure
%==========================================================================
load(strcat(fpath_variables,'electrode_location'));

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
        if sum(strcmp({dcs_magnitude{a,b,c,d,e}(:).name}',slices{a,b,c,d,e}(f).name))==0
%====================================== loop over individual slices

% load processed data
load(strcat(fpath_processed,slices{a,b,c,d,e}(f).name),...
    'pulset','indS','indD','baseIndex','fs','tpre')

% dcs onset time in samples for TBS protocol
dcs_on = 2*fs;

