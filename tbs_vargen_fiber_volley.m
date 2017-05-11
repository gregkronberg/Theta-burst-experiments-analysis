%% variable: fiber volleys
%==========================================================================
% store fiber volleys during baseline and induction
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
load(strcat(fpath_variables,'fiber_volley'));

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
        if sum(strcmp({fiber_volley{a,b,c,d,e}(:).name}',slices{a,b,c,d,e}(f).name))==0
%====================================== loop over individual slices
            % load processed data
            load(strcat(fpath_processed,slices{a,b,c,d,e}(f).name),...
                'indD4','baseD','indBlock',...
                'baseIndex','fs','pulset','tpre')
            
            %% time parameters
            %==============================================================
            % tbs induction
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
            t_fv = .0013*fs:.0025*fs;                                                    % time window for taking fiber volley
            
            % normalize traces
            baseD = baseD - ones(size(baseD,1),1)*baseD(t_fv(1),:);
            baseD_mean =  mean(baseD(:,indBlock(1)-20:indBlock(1)-1),2);% 
            for j = 1:size(indD4,3)
                indD4(:,:,j) = indD4(:,:,j) -  ones(size(indD4(:,:,j),1),1)*indD4(t_fv(1),:,j);
            end
            
            % take derivatives
            indD_diff1 = diff(indD4,1,1);% first derivative (time x pulses x bursts)
            baseD_diff1 = diff(baseD_mean,1,1);% first derivative (time x blocks)
            
            % find zero crossings
            baseD_sign = sign(baseD_diff1);% (time x blocks)
            baseD_sign(baseD_diff1==0)=1;% (time x blocks)
            baseD_zerocross = diff(baseD_sign)<0; 
            
            % choose peak time to measure fiber volley amplitude
            t2 =  min(find(baseD_zerocross));
            t2 = t2 +t_fv(1);
            t1 = .001*fs;
            t3 = t1 +(t2-t1)/2;
            
            % measure baseline fiber volley
            base_fv = zeros(size(baseD,2),1);
            for j = 1:size(baseD,2)
                fitx = [ones(2,1),[t1;t2]];                 % x values to fit a line between the positive peaks
                fity = baseD([round(pulset+t1);round(pulset+t2)],j);                                      % peak values to fit line
                fit = fitx\fity;                                        % fit line between the two peaks
                top = fit(1) + fit(2)*round(t3);                         % extrapolate from negative peak to fitted line
                bottom = baseD(pulset+round(t3),j);           % value at negative peak
                base_fv(j) = top-bottom;                             % fiber volley amplitude
            end
            base_fv_mean = mean(base_fv(indBlock(1)-tpre:indBlock(1)-1));
            
            % induction fiber volleys
            ind_fv = zeros(size(indD4,2),size(indD4,3));
            for j = 1:size(indD4,2)
                for k =1:size(indD4,3)
                    fitx = [ones(2,1),[t1;t2]];                 % x values to fit a line between the positive peaks
                    fity = indD4([round(t1);round(t2)],j,k);                                      % peak values to fit line
                    fit = fitx\fity;                                        % fit line between the two peaks
                    top = fit(1) + fit(2)*round(t3);                         % extrapolate from negative peak to fitted line
                    bottom = indD4(round(t3),j,k);           % value at negative peak
                    ind_fv(j,k) = top-bottom;                             % fiber volley amplitude
                end
            end
            ind_fv_norm = ind_fv/base_fv_mean;
            
            %% store fepsp_area variable
            %==================================================================
            fiber_volley{a,b,c,d,e}(f).fiber_volley = ind_fv;
            fiber_volley{a,b,c,d,e}(f).fiber_volley_norm = ind_fv_norm;
            fiber_volley{a,b,c,d,e}(f).fiber_volley_base = base_fv;
            fiber_volley{a,b,c,d,e}(f).fiber_volley_base_mean = base_fv_mean;
            fiber_volley{a,b,c,d,e}(f).fiber_volley_t = [t1 t3 t2];
            fiber_volley{a,b,c,d,e}(f).baseD = baseD;
            fiber_volley{a,b,c,d,e}(f).indD4 = indD4;
            fiber_volley{a,b,c,d,e}(f).slices = slices{a,b,c,d,e}(f);
            fiber_volley{a,b,c,d,e}(f).name = slices{a,b,c,d,e}(f).name;

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

save(strcat(fpath_variables,'fiber_volley.mat'),'fiber_volley')