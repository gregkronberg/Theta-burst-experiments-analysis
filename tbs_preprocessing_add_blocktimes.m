clear all; close all; clc

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
load(strcat(fpath_variables,'slices.mat'));

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
        %====================================== loop over individual slices
        load(strcat(fpath_raw,slices{a,b,c,d,e}(f).name))
        
        date_string = datestr(blocktimes);
        time = str2num(date_string(:,end-7:end-6)) + str2num(date_string(:,end-4:end-3))./60;
        slices{a,b,c,d,e}(f).blocktimes = time;
    end
end
                end
            end
        end
    end
end

save(strcat(fpath_variables,'slices.mat'),'slices','conditions')