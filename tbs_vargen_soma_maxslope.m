%% variable: max slope in somatic recordings during induction
%==========================================================================
% notes
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
load(strcat(fpath_variables,'slices'));

%% load global slopes structure
%==========================================================================
load(strcat(fpath_variables,'soma_maxslope'));

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
        if sum(strcmp({soma_maxslope{a,b,c,d,e}(:).name}',slices{a,b,c,d,e}(f).name))==0
%====================================== loop over individual slices
            % load processed data
            load(strcat(fpath_processed,slices{a,b,c,d,e}(f).name),...
                'indS_filt_band_200_600','baseS_filt_band_200_600','indBlock','indD4','baseIndex','fs','pulset')
            
            % tbs induction
            %==============================================================
            fs = 10000;
            indTime = 7*fs;                                                             % total time of induction block(samples)
            tbsOn = 2*fs;                                                               % tbs start time in s
            tbsOff = 5*fs;                                                              % tbs end time in s
            tBurst = .04*fs;                                                            % duration of each burst
            nBurst = 15;                                                                % number of burst
            nPulse = 4;                                                                 % number of pulses in a burst
            rBurst = 100;                                                               % burst frequency
            bipolarWidth = 2.5e-3*fs;                                                       % width of bipolar stimulus
            slope_window = bipolarWidth:(fs/rBurst-bipolarWidth);
            
            % cropped soma baseline trace
            baseS_filt_crop = baseS_filt_band_200_600((pulset+1):pulset+(fs/rBurst),indBlock(1)-20:indBlock(1)-1); % (time x blocks)
            
            % Make each column the time series of a single burst
            indS_filt_crop = indS_filt_band_200_600(tbsOn+1:tbsOff); % traces during TBS (time)
            indS_filt_crop = reshape(indS_filt_crop,[],nBurst); % reshape (time x bursts)

            % remove time between bursts
            indS_filt_crop = indS_filt_crop(1:tBurst,:); % (time x bursts)
            
            % reshape
            indS_filt_crop = reshape(indS_filt_crop,[],nPulse,nBurst); % (time x pulse x burst)
            
            % take derivatives
            baseS_filt_diff = diff(baseS_filt_crop,1); %(time x blocks) 
            indS_filt_diff = diff(indS_filt_crop,1);% (time x pulse x burst)
            
            % find max slope and index of max slope
            %===============================================================
            % baseline
            [base_max,base_max_i] = max(-baseS_filt_diff(slope_window,:),[],1);
            base_max = squeeze(base_max);% (blocks)
            base_max_i = squeeze(base_max_i)+bipolarWidth;% (blocks)
            base_max_mean = mean(base_max);% (scalar)
            
            % induction
            [ind_max,ind_max_i] = max(-indS_filt_diff(slope_window,:,:),[],1);
            ind_max = squeeze(ind_max);% (pulses x bursts)
            ind_max_i = squeeze(ind_max_i)+bipolarWidth;% (pulses x bursts)
            ind_max_norm = ind_max/base_max_mean;% (pulses x bursts)
            
            %% store fepsp_area variable
            %==================================================================
            soma_maxslope{a,b,c,d,e}(f).soma_maxslope = ind_max;
            soma_maxslope{a,b,c,d,e}(f).soma_maxslope_norm = ind_max_norm;
            soma_maxslope{a,b,c,d,e}(f).soma_maxslope_index = ind_max_i;
            soma_maxslope{a,b,c,d,e}(f).base_maxslope = base_max;
            soma_maxslope{a,b,c,d,e}(f).base_max_index = base_max_i;
            soma_maxslope{a,b,c,d,e}(f).base_trace = baseS_filt_crop;
            soma_maxslope{a,b,c,d,e}(f).ind_trace = indS_filt_crop;
            soma_maxslope{a,b,c,d,e}(f).slices = slices{a,b,c,d,e}(f);
            soma_maxslope{a,b,c,d,e}(f).name = slices{a,b,c,d,e}(f).name;
             
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

save(strcat(fpath_variables,'soma_maxslope.mat'),'soma_maxslope')