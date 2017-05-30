%% create new variable 

clear all; close all; clc

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
else
    % laptop
    fpath_raw = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
    fpath_processed = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
    fpath_variables = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
    fpath_analysis = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
    fpath_filters = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters
    fpath_processed_images = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Images\'; % processed
end

%% load global slices structure
%==========================================================================
load(strcat(fpath_variables,'slices'));

for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})

                    slopes{a,b,c,d,e}.slices = [];
                    slopes{a,b,c,d,e}.name = 'temp';
                    
                end
            end
        end
    end
end

save(strcat(fpath_variables,'slopes.mat'),'slopes');
   
    