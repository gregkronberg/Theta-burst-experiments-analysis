%% analysis: slopes
%==========================================================================
% plot slopes over time for each condtion
%==========================================================================

clear all
close all

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
    % laptop paths
    fpath_raw = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
    fpath_processed = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
    fpath_variables = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
    fpath_analysis = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
    fpath_filters = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters
    fpath_processed_images = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Images\'; % processed
end

% load global structures
%---------------------------------
load(strcat(fpath_variables,'slices.mat')); % slices
load(strcat(fpath_variables,'slopes.mat')); % slopes

% exclusion criteria
%----------------------------------
date_cut = [ 0 20170401 0];%[20170115 20170301 20170401]; % cutoff dates for [apical basal perforant]

% time
%----------------------------------
tpre = 20;
tpost = 60;
t  = 1:tpre+tpost;

% figure parameters
%----------------------------------
stim_color = {[0 0 0],[0 .5 1],[1 0 .5]};

%% store slopes for each condition in matrix
%==========================================================================
% loop over experimental conditions
%--------------------------------
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})                    

% loop over slices
%--------------------------------
if isempty(slices{a,b,c,d,e})==0
    for f = 1:length(slopes{a,b,c,d,e})

% index for baseline responses
base_index = slopes{a,b,c,d,e}(f).baseIndex;

% remove noisy samples
%==========================================================================
% store cleaned traces in new vector
slopes{a,b,c,d,e}(f).slopes_clean = slopes{a,b,c,d,e}(f).slopes;
slopes{a,b,c,d,e}(f).spikes_clean = slopes{a,b,c,d,e}(f).spike;
% window size to compare against neighboring samples
win = 10; 
% tolerance for deviation from neighboring samples (in standard deviations)
std_tol = 2;

% remove pre-induction noise
%-------------------------------------
% slopes
replace = remove_noise(slopes{a,b,c,d,e}(f).slopes_clean(base_index(1:tpre)),win,std_tol);
slopes{a,b,c,d,e}(f).slopes_clean(base_index(find(replace)))  = slopes{a,b,c,d,e}(f).slopes_clean(base_index(find(replace)-1));
% spikes
replace = remove_noise(slopes{a,b,c,d,e}(f).spikes_clean(1:tpre),win,std_tol);
slopes{a,b,c,d,e}(f).spikes_clean(find(replace))  = slopes{a,b,c,d,e}(f).spikes_clean(find(replace)-1);

% remove post-induction noise
%--------------------------------------
% slopes
replace = remove_noise(slopes{a,b,c,d,e}(f).slopes_clean(base_index(tpre+1:tpre+tpost)),win,std_tol);
slopes{a,b,c,d,e}(f).slopes_clean(base_index(tpre+find(replace)))  = slopes{a,b,c,d,e}(f).slopes_clean(base_index(tpre+find(replace)-1));
% spikes
replace = remove_noise(slopes{a,b,c,d,e}(f).spikes_clean(tpre+1:tpre+tpost),win,std_tol);
slopes{a,b,c,d,e}(f).spikes_clean(tpre+find(replace))  = slopes{a,b,c,d,e}(f).spikes_clean(tpre+find(replace)-1);

% store normalized slopes
%---------------------------------------
slopes{a,b,c,d,e}(f).slopes_base_mean = mean(slopes{a,b,c,d,e}(f).slopes_clean(base_index(1:tpre)));% (scalar)
slopes{a,b,c,d,e}(f).slopes_norm = slopes{a,b,c,d,e}(f).slopes_clean(base_index)/slopes{a,b,c,d,e}(f).slopes_base_mean; % (blocks)
slopes{a,b,c,d,e}(f).ind_slopes_norm  = slopes{a,b,c,d,e}(f).indSlopes/slopes{a,b,c,d,e}(f).slopes_base_mean;% (induction pulses)

% store normalized spikes
%----------------------------------------
slopes{a,b,c,d,e}(f).spikes_base_mean = mean(slopes{a,b,c,d,e}(f).spikes_clean((1:tpre)));% (scalar)
slopes{a,b,c,d,e}(f).spikes_norm = slopes{a,b,c,d,e}(f).spikes_clean/slopes{a,b,c,d,e}(f).spikes_base_mean; % (blocks)

% normalized baseline spike:epsp ratio
%----------------------------------------
slopes{a,b,c,d,e}(f).spikes_slopes_ratio = slopes{a,b,c,d,e}(f).spikes_clean./slopes{a,b,c,d,e}(f).slopes_clean(base_index);
slopes{a,b,c,d,e}(f).spikes_slopes_ratio_base = mean(slopes{a,b,c,d,e}(f).spikes_slopes_ratio(1:tpre));
slopes{a,b,c,d,e}(f).spikes_slopes_ratio_norm = slopes{a,b,c,d,e}(f).spikes_slopes_ratio/slopes{a,b,c,d,e}(f).spikes_slopes_ratio_base;

% end loop over individual slices
%----------------------------------------
    end
end
% end loop over experimental conditions
%----------------------------------------
                end
            end
        end
    end
end

%% fit traces after induction and apply rejection criteria
%==========================================================================
% vector for fitting curves after induction
x = [0:59]'; 

% loop over experimental conditions
%---------------------------------
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
% loop over slices
%----------------------------------
if isempty(slices{a,b,c,d,e})==0
    cnt = 0;
        for f = 1:length(slices{a,b,c,d,e})

% data to fit
%----------------------------------
y = abs(slopes{a,b,c,d,e}(f).slopes_norm(tpre+1:tpre+60)); 

% linear fit
%----------------------------------
lin_fit = polyfit(x,y,1);
slopes{a,b,c,d,e}(f).lin_fit_coef = lin_fit';
slopes{a,b,c,d,e}(f).lin_fit = lin_fit(1)*x+lin_fit(2);
slopes{a,b,c,d,e}(f).lin_fit_error = y - slopes{a,b,c,d,e}(f).lin_fit;
slopes{a,b,c,d,e}(f).max = max(y);

% exponential fit
%----------------------------------
exp_fit = fit(x,y,'exp1');
slopes{a,b,c,d,e}(f).exp_fit = [exp_fit.a;exp_fit.b];
slopes{a,b,c,d,e}(f).exp_fit_post = exp_fit.a*exp(exp_fit.b*x);
slopes{a,b,c,d,e}(f).exp_fit_error = mean(y - slopes{a,b,c,d,e}(f).exp_fit_post);

% reject slices
%----------------------------------

% slopes{a,b,c,d,e}(f).reject = exp_fit.b > 0 | exp_fit.b < -5e-3;
% slopes{a,b,c,d,e}(f).reject = exp_fit.b> 0 |...
%     lin_fit(1)*90 + lin_fit(2) < 1 &...
%     [slices{a,b,c,d,e}(f).date]'>date_cut(d);

slopes{a,b,c,d,e}(f).reject = [slices{a,b,c,d,e}(f).date]'<date_cut(d)|...
    exp_fit.b > -0e-3 |...
    exp_fit.b < -8e-3  ...
    ;

% end loop over slices
%----------------------------------
    end
end
% end loop over conditions
%----------------------------------
                end
            end
        end
    end
end

%% save slopes variable
%==========================================================================
save(strcat(fpath_variables,'slopes.mat'),'slopes')

%% Stats
%==========================================================================
% loop over conditions
%--------------------------------
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0

% slices to keep/reject
%-------------------------------
reject = [slopes{a,b,c,d,e}(:).reject]';
keep{a,b,c,d,e}  = ~reject ;%
keep_i{a,b,c,d,e} = find(keep{a,b,c,d,e});

if isempty(slices{a,b,c,d,e}(keep{a,b,c,d,e}))==0

% slopes mean and standard error
%-------------------------------
slopes_norm_temp = [slopes{a,b,c,d,e}(keep{a,b,c,d,e}).slopes_norm];
slopes_mean{a,b,c,d,e} = mean(slopes_norm_temp,2);
slopes_std{a,b,c,d,e} = std(slopes_norm_temp,0,2);
slopes_sem{a,b,c,d,e} = std(slopes_norm_temp,0,2)/sqrt(size(slopes_norm_temp,2));
slopes_end{a,b,c,d,e} = mean(slopes_norm_temp(end-9:end,:),1);
slopes_end_mean{a,b,c,d,e} = mean(slopes_end{a,b,c,d,e});
slopes_end_std{a,b,c,d,e} = std(slopes_end{a,b,c,d,e},0,2);
slopes_end_sem{a,b,c,d,e} = std(slopes_end{a,b,c,d,e},0,2)/sqrt(length(slopes_end{a,b,c,d,e}));

% spikes ratio mean and standard error
%-------------------------------
spikes_ratio_norm_temp = [slopes{a,b,c,d,e}(keep{a,b,c,d,e}).spikes_slopes_ratio_norm];
spikes_ratio_norm_temp(:,isinf(spikes_ratio_norm_temp(end,:)))=[];
spikes_ratio_mean{a,b,c,d,e} = mean(spikes_ratio_norm_temp,2);
spikes_ratio_std{a,b,c,d,e} = std(spikes_ratio_norm_temp,0,2);
spikes_ratio_sem{a,b,c,d,e} = std(spikes_ratio_norm_temp,0,2)/sqrt(size(spikes_ratio_norm_temp,2));
spikes_ratio_end{a,b,c,d,e} = mean(spikes_ratio_norm_temp(end-9:end,:),1);
spikes_ratio_end_mean{a,b,c,d,e} = mean(spikes_ratio_end{a,b,c,d,e});
spikes_ratio_end_std{a,b,c,d,e} = std(spikes_ratio_end{a,b,c,d,e},0,2);
spikes_ratio_end_sem{a,b,c,d,e} = std(spikes_ratio_end{a,b,c,d,e},0,2)/sqrt(length(slopes_end{a,b,c,d,e}));

% pairwise comparisons to control (1st parameter of condition b)
%--------------------------------
[h,spikes_ratio_p{a,b,c,d,e}] = ttest2(spikes_ratio_end{a,1,1,d,e},spikes_ratio_end{a,b,c,d,e});
[h,slopes_p{a,b,c,d,e}] = ttest2(slopes_end{a,1,1,d,e},slopes_end{a,b,c,d,e});
end

% end loop over conditions
%--------------------------------
                    end
                end
            end
        end
    end
end

%% figure: slopes over time
%==========================================================================
% loop over conditions
%--------------------------------
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices{a,b,c,d,e}(keep{a,b,c,d,e}))==0
% figure
%--------------------------------
figure;hold on
% dcs
errorbar(t,slopes_mean{a,b,c,d,e},slopes_sem{a,b,c,d,e},...
    '.','Color',stim_color{b},'MarkerSize',30);

errorbar(t,slopes_mean{a,2,3,d,e},slopes_sem{a,2,3,d,e},...
    '.','Color',stim_color{2},'MarkerSize',30);
% control
errorbar(t,slopes_mean{a,1,1,d,e},slopes_sem{a,1,1,d,e},...
    '.','Color',stim_color{1},'MarkerSize',30);

% format figure
%---------------------------------
run('figure_format_slopes')
xlabel('Time (min)','FontSize',30,'FontWeight','bold')
ylabel('Normalize fEPSP slope','FontSize',30,'FontWeight','bold')
title(strcat('TBS with ',conditions{2}{b},', ',num2str(conditions{3}(c)),...
    'V/m, ',conditions{4}{d},', p = ',num2str(slopes_p{a,b,c,d,e})));

% figure
%--------------------------------
figure;hold on
% dcs
errorbar(t,spikes_ratio_mean{a,b,c,d,e},spikes_ratio_sem{a,b,c,d,e},...
    '.','Color',stim_color{b},'MarkerSize',30);
% control
errorbar(t,spikes_ratio_mean{a,1,1,d,e},spikes_ratio_sem{a,1,1,d,e},...
    '.','Color',stim_color{1},'MarkerSize',30);

% format figure
%---------------------------------
run('figure_format_slopes')
xlabel('Time (min)','FontSize',30,'FontWeight','bold')
ylabel('pop spike:epsp slope ratio','FontSize',30,'FontWeight','bold')
title(strcat('TBS with ',conditions{2}{b},', ',num2str(conditions{3}(c)),...
    'V/m, ',conditions{4}{d},', p = ',num2str(spikes_ratio_p{a,b,c,d,e})));

                        end
                    end
                end
            end
        end
    end
end

%% figure: final plasticity comparing within each day
%==========================================================================
plot_location = [0 -1 1];

% loop over conditions
%--------------------------------
for d = 1:length(conditions{4})
    figure;hold on
    for a = 1:length(conditions{1})
        for b = 1:length(conditions{2})
            for c = [1,3];%:length(conditions{3})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices{a,b,c,d,e}(keep{a,b,c,d,e}))==0

% plot final plasticity for each slice
%----------------------------------
plot(plot_location(b)*ones(length(slopes_end{a,b,c,d,e})),slopes_end{a,b,c,d,e},'.',...
    'Color',stim_color{b},'MarkerSize',30)
title(strcat('TBS with',num2str(conditions{3}(c)),'V/m, ',conditions{4}{d}));
xlim([-2 2])

% draw line between data points for same day
%-----------------------------------
% kept slices
slices_keep{a,b,c,d,e} = slices{a,b,c,d,e}(keep{a,b,c,d,e});
% dates of all kept CONTROL slices
date_list = [slices_keep{a,1,1,d,e}.date]'; 
% loop over kept slices for all conditions
for f = 1:length(slices_keep{a,b,c,d,e})
    % check if there is a control slice with same date
    date_match = date_list == slices_keep{a,b,c,d,e}(f).date;
    % get plasticity of matching control
    slopes_match = slopes_end{a,1,1,d,e}(date_match);
    % plot a line between the matches
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
% loop over conditions
%--------------------------------
for d = 1:length(conditions{4})
    figure;hold on
    for a = 1:length(conditions{1})
        for b = 1:length(conditions{2})
            for c = [1,3];%:length(conditions{3})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices{a,b,c,d,e}(keep{a,b,c,d,e}))==0
                        
plot(abs([slopes{a,b,c,d,e}(keep{a,b,c,d,e}).slopes_base_mean]),slopes_end{a,b,c,d,e},'.','Color',...
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
% loop over conditions
%--------------------------------
for d = 1:length(conditions{4})
    figure;hold on
    for a = 1:length(conditions{1})
        for b = 1:length(conditions{2})
            for c = [1,3];%:length(conditions{3})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices{a,b,c,d,e}(keep{a,b,c,d,e}))==0

plot(abs([slopes{a,b,c,d,e}(keep{a,b,c,d,e}).spikes_base_mean]),slopes_end{a,b,c,d,e},'.','Color',...
    stim_color{b},'MarkerSize',30)
[R,P] = corrcoef([[slopes{a,b,c,d,e}(keep{a,b,c,d,e}).spikes_base_mean]',slopes_end{a,b,c,d,e}']);
xlabel('Baseline pop spike amplitude')
ylabel('plasticity')
title(strcat('TBS with',num2str(conditions{3}(c)),'V/m, ',conditions{4}{d},...
    ', R = ',num2str(R(1,2)),'P = ',num2str(P(1,2))));
                        end
                    end
                end
            end
        end
    end
end

%% plasticity as a function of baseline spike/slope ratio
%==========================================================================
% loop over conditions
%--------------------------------
for d = 1:length(conditions{4})
    figure;hold on
    for a = 1:length(conditions{1})
        for b = 1:length(conditions{2})
            for c = [1,3];%:length(conditions{3})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        if isempty(slices{a,b,c,d,e}(keep{a,b,c,d,e}))==0

plot(abs([slopes{a,b,c,d,e}(keep{a,b,c,d,e}).spikes_slopes_ratio_base]),...
    slopes_end{a,b,c,d,e},'.','Color',...
    stim_color{b},'MarkerSize',30);
[R,P] = corrcoef([[slopes{a,b,c,d,e}(keep{a,b,c,d,e}).spikes_slopes_ratio_base]',slopes_end{a,b,c,d,e}']);
xlabel('Baseline pop spikeamplitude:fEPSP slope ratio')
ylabel('plasticity')
title(strcat('TBS with',num2str(conditions{3}(c)),'V/m, ',conditions{4}{d},...
    ', R = ',num2str(R(1,2)),'P = ',num2str(P(1,2))));

                        end
                    end
                end
            end
        end
    end
end

%% plot height vs plasticity
figure;hold on
for b = 1:length(conditions{2})
    for c = [1,3];
        if isempty(slices{1,b,c,2,1}(keep{1,b,c,2,1}))==0
            plot([slices{1,b,c,2,1}(keep{1,b,c,2,1}).height],...
                slopes_end{1,b,c,2,1},'.','Color',stim_color{b},'MarkerSize',20)
        end
    end
end