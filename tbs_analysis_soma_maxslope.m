%% analysis: soma maxslope
%==========================================================================
% notes
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
load(strcat(fpath_variables,'slices.mat')); % slices
load(strcat(fpath_variables,'slopes.mat')); % slopes
load(strcat(fpath_variables,'fepsp_area.mat')); % slopes
load(strcat(fpath_variables,'soma_maxslope.mat')); % slopes

%% exclusion criteria
%==========================================================================
date_cut = [ 0 20170401 0];%[20170115 20170301 20170401]; % cutoff dates for [apical basal perforant]
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        % apply inclusion/exclusion criteria
                        include{a,b,c,d,e} = [slices{a,b,c,d,e}(:).date]'>date_cut(d);
                        slices_temp{a,b,c,d,e} =  slices{a,b,c,d,e}(include{a,b,c,d,e});
                    end
                end
            end
        end
    end
end

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
    for f = 1:length(slices{a,b,c,d,e})
%====================================== loop over individual slices

% max slope
soma_maxslope{a,b,c,d,e}(f).maxslope_norm = mean(soma_maxslope{a,b,c,d,e}(f).soma_maxslope_norm,1)';
soma_maxslope{a,b,c,d,e}(f).maxslope_mean = mean(soma_maxslope{a,b,c,d,e}(f).maxslope_norm);
soma_maxslope{a,b,c,d,e}(f).maxslope_max = max(soma_maxslope{a,b,c,d,e}(f).maxslope_norm);
soma_maxslope{a,b,c,d,e}(f).maxslope_adapt = mean(soma_maxslope{a,b,c,d,e}(f).maxslope_norm(end))...
    /mean(soma_maxslope{a,b,c,d,e}(f).maxslope_norm(1));

% timing
soma_maxslope{a,b,c,d,e}(f).maxslope_i = soma_maxslope{a,b,c,d,e}(f).soma_maxslope_index(:);
soma_maxslope{a,b,c,d,e}(f).maxslope_base_i = mean(soma_maxslope{a,b,c,d,e}(f).base_max_index);
soma_maxslope{a,b,c,d,e}(f).maxslope_i_norm = soma_maxslope{a,b,c,d,e}(f).maxslope_i - ...
    soma_maxslope{a,b,c,d,e}(f).maxslope_base_i;

soma_maxslope{a,b,c,d,e}(f).maxslope_i_mean = mean(soma_maxslope{a,b,c,d,e}(f).maxslope_i_norm);
soma_maxslope{a,b,c,d,e}(f).maxslope_i_max = max(soma_maxslope{a,b,c,d,e}(f).maxslope_i_norm);
soma_maxslope{a,b,c,d,e}(f).maxslope_i_adapt = mean(soma_maxslope{a,b,c,d,e}(f).maxslope_i_norm(end))/...
    mean(soma_maxslope{a,b,c,d,e}(f).maxslope_i_norm(1));
soma_trace{a,b,c,d,e}(:,:,f) = reshape(soma_maxslope{a,b,c,d,e}(f).ind_trace,[],60);
soma_maxslope{a,b,c,d,e}(f).reject = [slices{a,b,c,d,e}(f).date]'<date_cut(d);
%===================================== end loop over individual slices
    end
    reject = [soma_maxslope{a,b,c,d,e}(:).reject]';
    keep{a,b,c,d,e}  = ~reject ;%
    keep_i{a,b,c,d,e} = find(keep{a,b,c,d,e});
end
%===================================== end loop over experimental conditions
                end
            end
        end
    end
end

save(strcat(fpath_variables,'soma_maxslope.mat'),'soma_maxslope')

%% plot maxslope for each condition
figure;hold on
plot(1*ones(length(soma_maxslope{1,1,1,2,1}(keep_i{1,1,1,2,1})),1),[soma_maxslope{1,1,1,2,1}(keep_i{1,1,1,2,1}).maxslope_mean],'.k','MarkerSize',25)
plot(2*ones(length(soma_maxslope{1,3,3,2,1}(keep_i{1,3,3,2,1})),1),[soma_maxslope{1,3,3,2,1}(keep_i{1,3,3,2,1}).maxslope_mean],'.r','MarkerSize',25)
plot(0*ones(length(soma_maxslope{1,2,3,2,1}(keep_i{1,2,3,2,1})),1),[soma_maxslope{1,2,3,2,1}(keep_i{1,2,3,2,1}).maxslope_mean],'.b','MarkerSize',25)
title('mean maxslope')
xlim([-2 4])

figure;hold on
plot(1*ones(length(soma_maxslope{1,1,1,2,1}(keep_i{1,1,1,2,1})),1),[soma_maxslope{1,1,1,2,1}(keep_i{1,1,1,2,1}).maxslope_adapt],'.k','MarkerSize',25)
plot(2*ones(length(soma_maxslope{1,3,3,2,1}(keep_i{1,3,3,2,1})),1),[soma_maxslope{1,3,3,2,1}(keep_i{1,3,3,2,1}).maxslope_adapt],'.r','MarkerSize',25)
plot(0*ones(length(soma_maxslope{1,2,3,2,1}(keep_i{1,2,3,2,1})),1),[soma_maxslope{1,2,3,2,1}(keep_i{1,2,3,2,1}).maxslope_adapt],'.b','MarkerSize',25)
title('maxslope adaption (last/first)')
xlim([-2 4])

figure;hold on
plot(1*ones(length(soma_maxslope{1,1,1,2,1}(keep_i{1,1,1,2,1})),1),[soma_maxslope{1,1,1,2,1}(keep_i{1,1,1,2,1}).maxslope_max],'.k','MarkerSize',25)
plot(2*ones(length(soma_maxslope{1,3,3,2,1}(keep_i{1,3,3,2,1})),1),[soma_maxslope{1,3,3,2,1}(keep_i{1,3,3,2,1}).maxslope_max],'.r','MarkerSize',25)
plot(0*ones(length(soma_maxslope{1,2,3,2,1}(keep_i{1,2,3,2,1})),1),[soma_maxslope{1,2,3,2,1}(keep_i{1,2,3,2,1}).maxslope_max],'.b','MarkerSize',25)
title('max maxslope')
xlim([-2 4])

%%
[h,p] = ttest2([soma_maxslope{1,1,1,2,1}(keep_i{1,1,1,2,1}).maxslope_mean],[soma_maxslope{1,3,3,2,1}(keep_i{1,3,3,2,1}).maxslope_mean])