%% analysis: electrode location
%==========================================================================
% plot slopes over time for each condtion
%==========================================================================

clear all
close all
clc

% file paths
%==========================================================================
% desktop
fpath_raw = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
fpath_processed = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
fpath_analysis = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
fpath_filters = 'D:\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters


% load global structures
%==========================================================================
load(strcat(fpath_variables,'slices')); % slices
load(strcat(fpath_variables,'slopes')); % slopes
load(strcat(fpath_variables,'electrode_location')); % slopes

% exclusion criteria
%==========================================================================
date_cut = [0 0 0];%[20170115 20170301 20170401]; % cutoff dates for [apical basal perforant]

% preallocate
%==========================================================================
tpre = 20;
tpost = 60;
t  = 1:tpre+tpost;
stim_color = {[0 0 0],[0 0 1],[1 0 0]};

%% store slopes for each condition in matrix
%==========================================================================
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})                    
%===================================== loop over experimental conditions
% check for slices
if isempty(slices{a,b,c,d,e})==0
    % apply exclusion criteria
    include = [slices{a,b,c,d,e}(:).date]'>date_cut(d);
    slices_temp{a,b,c,d,e} = slices{a,b,c,d,e}(include);
    slopes_temp{a,b,c,d,e} = slopes{a,b,c,d,e}(include);
    
    % preallocate
    for f = 1:length(slopes_temp{a,b,c,d,e})
%====================================== loop over individual slices

% index for baseline responses
base_index = slopes_temp{a,b,c,d,e}(f).baseIndex;

% store normalized slopes as a matrix
slopes_base_mean{a,b,c,d,e}(f) = mean(slopes_temp{a,b,c,d,e}(f).slopes(base_index(1:tpre)));
slopes_norm{a,b,c,d,e}(:,f) = slopes_temp{a,b,c,d,e}(f).slopes(base_index)'/slopes_base_mean{a,b,c,d,e}(f); % (blocks x slices)
ind_slopes_norm{a,b,c,d,e}(:,f)  = slopes_temp{a,b,c,d,e}(f).indSlopes/slopes_base_mean{a,b,c,d,e}(f);

% store normalized spikes as a matrix
spikes_base_mean{a,b,c,d,e}(f) = mean(slopes_temp{a,b,c,d,e}(f).spike((1:tpre)));
spikes_norm{a,b,c,d,e}(:,f) = slopes_temp{a,b,c,d,e}(f).spike'/spikes_base_mean{a,b,c,d,e}(f); % (blocks x slices)

    end
end
                end
            end
        end
    end
end

%% Stats
%==========================================================================
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0

% mean and standard error
slopes_mean{a,b,c,d,e} = mean(slopes_norm{a,b,c,d,e},2);
slopes_sem{a,b,c,d,e} = std(slopes_norm{a,b,c,d,e},0,2)/sqrt(size(slopes_norm{a,b,c,d,e},2));
slopes_end{a,b,c,d,e} = mean(slopes_norm{a,b,c,d,e}(end-10:end,:),1);
slopes_end_mean{a,b,c,d,e} = mean(slopes_end{a,b,c,d,e});
slopes_end_sem{a,b,c,d,e} = std(slopes_end{a,b,c,d,e},0,2)/sqrt(length(slopes_end{a,b,c,d,e}));

% pairwise comparisons to no DCS (1st parameter of condition b)
[h,slopes_p{a,b,c,d,e}] = ttest2(slopes_end{a,1,1,d,e},slopes_end{a,b,c,d,e});
                   
                    end
                end
            end
        end
    end
end

%% plot fiber volleys over time
%==========================================================================
for d = 1:length(conditions{4})
    figure;hold on
    for a = 1:length(conditions{1})
        for b = 1:length(conditions{2})
            for c = [1,3];%:length(conditions{3})
                    for e = 1:length(conditions{5})
                        if isempty(slices{a,b,c,d,e})==0
for f = 1:length(electrode_location{a,b,c,d,e})
    if isempty(electrode_location{a,b,c,d,e}(f).stim_to_soma)
        stim_soma{a,b,c,d,e}(f) = 0;
    else
        stim_soma{a,b,c,d,e}(f) = electrode_location{a,b,c,d,e}(f).stim_to_soma;
    end
end

%
plot((stim_soma{a,b,c,d,e}),slopes_end{a,b,c,d,e}(1:length(stim_soma{a,b,c,d,e})),'.','Color',...
stim_color{b},'MarkerSize',30)
ylim([0 4])
% 
% plot(abs(spikes_base_mean{a,b,c,d,e}(1:length(stim_soma{a,b,c,d,e}))./...
%     slopes_base_mean{a,b,c,d,e}(1:length(stim_soma{a,b,c,d,e}))),stim_soma{a,b,c,d,e},'.','Color',...
%     stim_color{b},'MarkerSize',30)
xlabel('Baseline pop spike amplitude')
ylabel('plasticity')
                        end
                    end
            end
        end
    end
end