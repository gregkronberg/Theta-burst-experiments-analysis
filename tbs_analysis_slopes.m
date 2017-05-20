%% analysis: slopes
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
load(strcat(fpath_variables,'slices.mat')); % slices
load(strcat(fpath_variables,'slopes.mat')); % slopes

%% exclusion criteria
%==========================================================================
date_cut = [ 0 0 0];%[20170115 20170301 20170401]; % cutoff dates for [apical basal perforant]
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
    for f = 1:length(slopes{a,b,c,d,e})
%====================================== loop over individual slices

% index for baseline responses
base_index = slopes{a,b,c,d,e}(f).baseIndex;

% store normalized slopes as a matrix
slopes{a,b,c,d,e}(f).slopes_base_mean = mean(slopes{a,b,c,d,e}(f).slopes(base_index(1:tpre)));% (scalar)
slopes{a,b,c,d,e}(f).slopes_norm = slopes{a,b,c,d,e}(f).slopes(base_index)/slopes{a,b,c,d,e}(f).slopes_base_mean; % (blocks)
slopes{a,b,c,d,e}(f).ind_slopes_norm  = slopes{a,b,c,d,e}(f).indSlopes/slopes{a,b,c,d,e}(f).slopes_base_mean;% (induction pulses)

% store normalized spikes as a matrix
slopes{a,b,c,d,e}(f).spikes_base_mean = mean(slopes{a,b,c,d,e}(f).spike((1:tpre)));% (scalar)
slopes{a,b,c,d,e}(f).spikes_norm = slopes{a,b,c,d,e}(f).spike/slopes{a,b,c,d,e}(f).spikes_base_mean; % (blocks)

% baseline spike:epsp ratio
slopes{a,b,c,d,e}(f).spikes_slopes_ratio_base = slopes{a,b,c,d,e}(f).spikes_base_mean/slopes{a,b,c,d,e}(f).slopes_base_mean;
        
%===================================== end loop over individual slices
    end
end
%===================================== end loop over experimental conditions
                end
            end
        end
    end
end
save(strcat(fpath_variables,'slopes.mat'),'slopes')

%% Stats
%==========================================================================
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        % apply inclusion criteria
                        if isempty(slices_temp{a,b,c,d,e})==0

% mean and standard error
slopes_norm_temp = [slopes{a,b,c,d,e}(include{a,b,c,d,e}).slopes_norm];
slopes_mean{a,b,c,d,e} = mean(slopes_norm_temp,2);
slopes_sem{a,b,c,d,e} = std(slopes_norm_temp,0,2)/sqrt(size(slopes_norm_temp,2));
slopes_end{a,b,c,d,e} = mean(slopes_norm_temp(end-10:end,:),1);
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
end

%% plot slopes over time
%==========================================================================
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices_temp{a,b,c,d,e})==0
figure;hold on
errorbar(t,slopes_mean{a,b,c,d,e},slopes_sem{a,b,c,d,e},...
    '.','Color',stim_color{b},'MarkerSize',30);
errorbar(t,slopes_mean{a,1,1,d,e},slopes_sem{a,1,1,d,e},...
    '.','Color',stim_color{1},'MarkerSize',30);

% format figure
run('figure_format_slopes')
xlabel('Time (min)','FontSize',30,'FontWeight','bold')
ylabel('Normalize fEPSP slope','FontSize',30,'FontWeight','bold')
title(strcat('TBS with ',conditions{2}{b},', ',num2str(conditions{3}(c)),...
    'V/m, ',conditions{4}{d},', p = ',num2str(slopes_p{a,b,c,d,e})));
                    
                        end
                    end
                end
            end
        end
    end
end

%% plot final plasticity comparing within each day
%==========================================================================
plot_location = [0 -1 1];
for d = 1:length(conditions{4})
    figure;hold on
    for a = 1:length(conditions{1})
        for b = 1:length(conditions{2})
            for c = [1,3];%:length(conditions{3})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices_temp{a,b,c,d,e})==0

% plot all final data points for each condition                            
plot(plot_location(b)*ones(length(slopes_end{a,b,c,d,e})),slopes_end{a,b,c,d,e},'.',...
    'Color',stim_color{b},'MarkerSize',30)
title(strcat('TBS with',num2str(conditions{3}(c)),'V/m, ',conditions{4}{d}));
xlim([-2 2])

% draw line between data points for same day
date_list = [slices_temp{a,1,1,d,e}(:).date]'; 
for f = 1:length(slices_temp{a,b,c,d,e})
    date_match = date_list == slices_temp{a,b,c,d,e}(f).date;
    slopes_match = slopes_end{a,1,1,d,e}(date_match);
    if isempty(slopes_match)==0
        plot([plot_location(1),plot_location(b)],[slopes_match(1),slopes_end{a,b,c,d,e}(f)],...
            'Color',stim_color{b})
    end
end
ylabel('Plasticity')
                        end
                    end  
                end
            end
        end
    end
end

%% plasticity as a function of baseline slopes
%==========================================================================
for d = 1:length(conditions{4})
    figure;hold on
    for a = 1:length(conditions{1})
        for b = 1:length(conditions{2})
            for c = [1,3];%:length(conditions{3})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices_temp{a,b,c,d,e})==0
plot(abs([slopes{a,b,c,d,e}(include{a,b,c,d,e}).slopes_base_mean]),slopes_end{a,b,c,d,e},'.','Color',...
    stim_color{b},'MarkerSize',30)
xlabel('Baseline fEPSP slope')
ylabel('plasticity')
title(strcat('TBS with',num2str(conditions{3}(c)),'V/m, ',conditions{4}{d}));
                        end
                    end
                end
            end
        end
    end
end

%% plasticity as a function of baseline spikes
%==========================================================================
for d = 1:length(conditions{4})
    figure;hold on
    for a = 1:length(conditions{1})
        for b = 1:length(conditions{2})
            for c = [1,3];%:length(conditions{3})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices_temp{a,b,c,d,e})==0
plot(abs([slopes{a,b,c,d,e}(include{a,b,c,d,e}).spikes_base_mean]),slopes_end{a,b,c,d,e},'.','Color',...
    stim_color{b},'MarkerSize',30)
xlabel('Baseline pop spike amplitude')
ylabel('plasticity')
title(strcat('TBS with',num2str(conditions{3}(c)),'V/m, ',conditions{4}{d}));
                        end
                    end
                end
            end
        end
    end
end

%% plasticity as a function of baseline spike/slope ratio
%==========================================================================
for d = 1:length(conditions{4})
    figure;hold on
    for a = 1:length(conditions{1})
        for b = 1:length(conditions{2})
            for c = [1,3];%:length(conditions{3})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices_temp{a,b,c,d,e})==0


plot(abs([slopes{a,b,c,d,e}(include{a,b,c,d,e}).spikes_slopes_ratio_base]),...
    slopes_end{a,b,c,d,e},'.','Color',...
    stim_color{b},'MarkerSize',30);
xlabel('Baseline pop spikeamplitude:fEPSP slope ratio')
ylabel('plasticity')
title(strcat('TBS with',num2str(conditions{3}(c)),'V/m, ',conditions{4}{d}));


                        end
                    end
                end
            end
        end
    end
end