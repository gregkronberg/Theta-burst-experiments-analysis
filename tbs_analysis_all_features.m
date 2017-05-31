%% analysis: all features
%==========================================================================
% 
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
load(strcat(fpath_variables,'slices.mat')); 
load(strcat(fpath_variables,'slopes.mat')); 
load(strcat(fpath_variables,'dcs_magnitude.mat')); 
load(strcat(fpath_variables,'drift.mat')); 
load(strcat(fpath_variables,'electrode_location.mat')); 
load(strcat(fpath_variables,'fiber_volley.mat')); 
load(strcat(fpath_variables,'soma_maxslope.mat'));
load(strcat(fpath_variables,'fepsp_area.mat'));

%% create feature list of all slices
%==========================================================================
features = []; 
cnt = 0;
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        for f = 1:length(slices{a,b,c,d,e})
                            if slopes{a,b,c,d,e}(f).reject ==0
cnt = cnt+1;
ind_block = slices{a,b,c,d,e}(f).indBlock(1);
% slices
% features(cnt).name = slices{a,b,c,d,e}(f).name;
features(cnt).date = slices{a,b,c,d,e}(f).date;
features(cnt).slice_num = slices{a,b,c,d,e}(f).slice_num;
features(cnt).location = d;
features(cnt).dcs_polarity = b;
features(cnt).dcs_intensity = c;
features(cnt).ind_time = slices{a,b,c,d,e}(f).blocktimes(ind_block);
features(cnt).age = slices{a,b,c,d,e}(f).age;
features(cnt).height = slices{a,b,c,d,e}(f).height;
if strcmp(slices{a,b,c,d,e}(f).hemi,'l')
    features(cnt).hemi = 1;
elseif strcmp(slices{a,b,c,d,e}(f).hemi,'r')
    features(cnt).hemi = 2;
else
    features(cnt).hemi = NaN;
end

% dcs magnitude
if length(dcs_magnitude{a,b,c,d,e})>1 
    if isempty(dcs_magnitude{a,b,c,d,e}(f).slices)==0
        features(cnt).dcs_measured = dcs_magnitude{a,b,c,d,e}(f).dcs_magnitude;
%         % electrode loacation
        features(cnt).stim_to_soma = electrode_location{a,b,c,d,e}(f).stim_to_soma;
        features(cnt).stim_to_dend = electrode_location{a,b,c,d,e}(f).stim_to_dend;
    else
        features(cnt).dcs_measured = NaN;
        features(cnt).stim_to_soma = NaN;
        features(cnt).stim_to_dend = NaN;
    end
else
    features(cnt).dcs_measured = NaN;
    features(cnt).stim_to_soma = NaN;
    features(cnt).stim_to_dend = NaN;
end

% Physiological features
% slopes
features(cnt).slopes_base_mean = slopes{a,b,c,d,e}(f).slopes_base_mean;
features(cnt).spikes_base_mean = slopes{a,b,c,d,e}(f).spikes_base_mean;
features(cnt).spikes_slopes_ratio_base = slopes{a,b,c,d,e}(f).spikes_slopes_ratio_base;

% % drift
features(cnt).drift = drift{a,b,c,d,e}(f).slopes_drift;
features(cnt).drift_time = drift{a,b,c,d,e}(f).time;

% fiber volley
features(cnt).fiber_volley = mean(fiber_volley{a,b,c,d,e}(f).fiber_volley_norm(1,:));

% burst area
features(cnt).burst_area_first = fepsp_area{a,b,c,d,e}(f).burst_area_norm(1);
features(cnt).burst_area_mean = fepsp_area{a,b,c,d,e}(f).burst_area_mean;
features(cnt).burst_area_max = fepsp_area{a,b,c,d,e}(f).burst_area_max;
features(cnt).burst_area_adapt = fepsp_area{a,b,c,d,e}(f).burst_area_adapt;

% soma max slope
features(cnt).soma_maxslope_first = soma_maxslope{a,b,c,d,e}(f).maxslope_norm(1);
features(cnt).soma_maxslope_mean = soma_maxslope{a,b,c,d,e}(f).maxslope_mean;
features(cnt).soma_maxslope_max = soma_maxslope{a,b,c,d,e}(f).maxslope_max;
features(cnt).soma_maxslope_adapt = soma_maxslope{a,b,c,d,e}(f).maxslope_adapt;
features(cnt).soma_maxslope_i_first = soma_maxslope{a,b,c,d,e}(f).maxslope_i_norm(1);
features(cnt).soma_maxslope_i_mean = soma_maxslope{a,b,c,d,e}(f).maxslope_i_mean;
features(cnt).soma_maxslope_i_max = soma_maxslope{a,b,c,d,e}(f).maxslope_i_max;
features(cnt).soma_maxslope_i_adapt = soma_maxslope{a,b,c,d,e}(f).maxslope_i_adapt;

% slopes
features(cnt).slopes_end = mean(slopes{a,b,c,d,e}(f).slopes_norm(end-9:end));
                        
                            end
                        end
                    end
                end
            end
        end
    end
end

%% store all features as matrix
%==========================================================================
% list features
feature_names = fieldnames(features);
feature_cell = cell(length(features),length(feature_names));
% store features in a cell
for a = 1:length(features)
    for b = 1:length(feature_names)
        feature_cell{a,b} = features(a).(feature_names{b}); % (examples x features)
    end
end

% convert to matrix
for a = 1:length(feature_names)
    % raw feature matrix
    feature_mat(:,a) = double(cell2mat(feature_cell(:,a))); % (examples x features)
end

save(strcat(fpath_variables,'features.mat'),'features','feature_mat')

%% canonical correlation
%==========================================================================
% gui for selecting features
%------------------------------------
% fig = figure;hold on
% fig.Position = [100 100 400 length(feature_names)*21]
% set(gca,'Visible','off')
% title('select features')
% for a = 1:length(feature_names)   
%     gui(a) = uicontrol('Style','checkbox','String',feature_names{a},'Position',...
%         [20 a*20 400 20],'Value',1);
% end
% 
% h = uicontrol('Position',[20 0 400 20],'String','Continue',...
%               'Callback','uiresume(gcbf)');
% uiwait(gcf);    
    
% % indices of features to include
% feature_include = logical([gui(:).Value])';

% code for selecting features
feature_list = feature_names;%...
%     {'date';...
%     'location';...
%     'dcs_measured';...
%     'spikes_base_mean';...
%     'fiber_volley';...
%     'burst_area_first';...
%     'burst_area_adapt';...
%     'soma_maxslope_first';...
%     'soma_maxslope_adapt';...
%     'slopes_end'};

feature_include = logical(zeros(length(feature_names),1));
for a = 1:length(feature_names)
    if sum(strcmp(feature_list,feature_names{a}))==1
        feature_include(a) = logical(1);
    end
end

% feature matrix
feature_mat_temp = feature_mat(:,feature_include);

% list kept features
feature_names_temp = {feature_names{feature_include}}';

% index of  dcs polarity feature
dcs_i = strcmp(feature_names_temp,'dcs_polarity'); % (features)

loc_feat_i = strcmp(feature_names_temp,'location');
% index of recording location feature(1 = apical, 2 = basal, 3 = perforant)
loc_i = feature_mat_temp(:,strcmp(feature_names_temp,'location'))==2; %(features)

% new target variable index
target_i = strcmp(feature_names_temp,'slopes_end');

% index of missing data 
nans = sum(isnan(feature_mat_temp) + isinf(feature_mat_temp) + (feature_mat_temp==0),2) > 0; % (examples)

% index of examples to keep for analysis
keep = ~nans & loc_i; % examples

% normalize (zero mean,unit variance)
feature_mat_norm = (feature_mat_temp(keep,:) - ...
    ones(size(feature_mat_temp(keep,:),1),1)*mean(feature_mat_temp(keep,:),1))./...
    (ones(size(feature_mat_temp(keep,:),1),1)*std(feature_mat_temp(keep,:),1));

% canonical correlation                           
[A,B,r,U,V,stats] = canoncorr(feature_mat_norm(:,~target_i&~loc_feat_i),...
    feature_mat_norm(:,target_i));

% figures
dcs_color = {[0 0 0],[0 0 1],[1 0 0]};% plot color for dcs conditions ([control,cathodal,anodal])
feature_mat_nonans = feature_mat_temp(keep,:);
figure;hold on
for a = 1:size(U,1)
    dcs_polarity  = feature_mat_nonans(a,dcs_i);
    plot(U(a),V(a),'.','Color',dcs_color{dcs_polarity},'MarkerSize',30)
end

% sort coefficient in ascending order
[A_sorted,A_sorted_i] = sort(abs(A));
% shift index due to target variable
A_sorted_i(A_sorted_i>=find(loc_feat_i)) = A_sorted_i(A_sorted_i>=find(loc_feat_i))+1;
% list features in ascending order
features_sorted = feature_names_temp(A_sorted_i);


%% plot single feature
feat = strcmp(feature_names_temp,'dcs_measured');
targ = strcmp(feature_names_temp,'slopes_end');
x = feature_mat_norm(:,feat);
y = feature_mat_norm(:,targ);
figure;hold on
plot(x,y,'.','MarkerSize',30);


%% regress out "noise" features
%% canonical correlation for multiple feature combinations
%=========================================================================
% number of features
% num = 5; % scalar
% for a = 1:num
%     :length(feature_names)
%     feature1 = strcmp(feature_names,feature_names{a});
%     for b = 1:length(feature_names())
                        
