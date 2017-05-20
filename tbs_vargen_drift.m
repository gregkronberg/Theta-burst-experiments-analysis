%% variable: manual fEPSP slopes and population spikes
%==========================================================================
% Notes
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
load(strcat(fpath_variables,'slices.mat'));
load(strcat(fpath_variables,'slopes.mat'));
load(strcat(fpath_variables,'drift.mat'));

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
        if sum(strcmp({drift{a,b,c,d,e}(:).name}',slices{a,b,c,d,e}(f).name))==0
%====================================== loop over individual slices

indBlock = slices{a,b,c,d,e}(f).indBlock;

slopes_pre = slopes{a,b,c,d,e}(f).slopes(1:indBlock(1)-1);

smooth_span = 9;
slopes_smooth = smooth(slopes_pre,smooth_span);

slopes_drift = mean(slopes_smooth(end-9:end))/mean(slopes_smooth(20:25));

 %% store slopes variable
            %==================================================================
            drift{a,b,c,d,e}(f).slopes_smooth = slopes_smooth;
            drift{a,b,c,d,e}(f).slopes_drift = slopes_drift;
            drift{a,b,c,d,e}(f).time = slices{a,b,c,d,e}(f).blocktimes(indBlock(1));
            drift{a,b,c,d,e}(f).slices = slices{a,b,c,d,e}(f);
            drift{a,b,c,d,e}(f).name = slices{a,b,c,d,e}(f).name;

            %===================================== end loop over individual slices
        end
    end
end
%===================================== end loop over experimental conditions
                end
            end
        end
    end
end
save(strcat(fpath_variables,'drift.mat'),'drift')