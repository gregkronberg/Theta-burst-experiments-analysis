%% variable: fEPSP area during induction
%==========================================================================
% notes
%==========================================================================

clear all
close all
clc

%% file paths
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

%% load global slices structure
%==========================================================================
load(strcat(fpath_variables,'slices'));

%% load global slopes structure
%==========================================================================
load(strcat(fpath_variables,'fEPSP_area'));

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
        if sum(strcmp({fepsp_area{a,b,c,d,e}(:).name}',slices{a,b,c,d,e}(f).name))==0
%====================================== loop over individual slices
            % load processed data
            load(strcat(fpath_processed,slices{a,b,c,d,e}(f).name),...
                'indD_filt_high_1','baseD_filt_high_1','indBlock','indD4','baseIndex','fs','pulset')
            
            % tbs induction
            fs = 10000;
            tbsOn = 2*fs;                                                               % tbs start time in s
            tbsOff = 5*fs;                                                              % tbs end time in s
            tBurst = .04*fs;                                                            % duration of each burst
            nBurst = 15;                                                                % number of burst
            nPulse = 4;                                                                 % number of pulses in a burst
            rBurst = 100;                                                               % burst frequency
            bipolarWidth = 2e-3*fs;                                                       % width of bipolar stimulus
            bipolarT  = ([0:fs/rBurst:(nPulse-1)*fs/rBurst]' + ...
                [1:bipolarWidth])';
            bipolarT = bipolarT(:);
            
            % cropped trace to region for taking area
            baseD_filt_crop = baseD_filt_high_1((pulset+bipolarWidth):pulset+(fs/rBurst),indBlock(1)-20:indBlock(1)-1); % (time x blocks)
            
            % normalize to first sample
            baseD_filt_crop_norm = baseD_filt_crop - ones(size(baseD_filt_crop,1),1)*baseD_filt_crop(1,:);
            
            % Area under fEPSP
            base_area = mean(sum(abs(baseD_filt_crop),1)); % average area under basleine fEPSP
            
            % Make each column the time series of a single burst
            % 3rd dimension is slice number
            indD_filt_crop = indD_filt_high_1(tbsOn+1:tbsOff); % traces during TBS (time)
            indD_filt_crop = reshape(indD_filt_crop,[],nBurst); % reshape (time x bursts)

            % remove time between bursts
            indD_filt_crop = indD_filt_crop(1:tBurst,:); % (time x burst)
            
            % reshape to (time x pulse x burst)
            indD_filt_crop = reshape(indD_filt_crop,[],nPulse,nBurst);
            
            % remove bipolar pulse
            indD_filt_crop = indD_filt_crop(bipolarWidth:end,:,:); % (time x burst x slice)
            
            % reshape to time x bursts
            indD_filt_crop = reshape(indD_filt_crop,[],nBurst);
            
            % Normalize first sample to zero
            indD_filt_crop_norm = indD_filt_crop-...
                ones(size(indD_filt_crop,1),1)*indD_filt_crop(1,:); % time x bursts
            
            % reshape to time x pulses x bursts
            indD_filt_crop_norm = reshape(indD_filt_crop_norm,[],nPulse,nBurst);
            
            % area under each pulse
            burst_area = squeeze(sum(abs(indD_filt_crop_norm),1)); % pulses x bursts
            
            % normalize to baseline
            burst_area_norm = burst_area/base_area; % pulse x bursts
            
            %% store fepsp_area variable
            %==================================================================
            fepsp_area{a,b,c,d,e}(f).fepsp_area = burst_area;
            fepsp_area{a,b,c,d,e}(f).base_area = base_area;
            fepsp_area{a,b,c,d,e}(f).fepsp_area_norm = burst_area_norm;
            fepsp_area{a,b,c,d,e}(f).base_trace = baseD_filt_crop_norm;
            fepsp_area{a,b,c,d,e}(f).ind_trace = indD_filt_crop_norm;
            fepsp_area{a,b,c,d,e}(f).slices = slices{a,b,c,d,e}(f);
            fepsp_area{a,b,c,d,e}(f).name = slices{a,b,c,d,e}(f).name;
            
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

save(strcat(fpath_variables,'fepsp_area.mat'),'fepsp_area')