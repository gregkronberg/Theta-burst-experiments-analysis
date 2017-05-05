%% variable: manual fEPSP slopes and population spikes
%==========================================================================
% Take slopes and population spikes of baseline respones and during
% induction.  Data are stored as a structure called slopes in the Matlab
% Variables folder.  slopes is organized according to the slices reference
% structure, where each entry contains the slopes and info for each slice
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
load(strcat(fpath_variables,'slopes'));

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
        if sum(strcmp({slopes{a,b,c,d,e}(:).name}',slices{a,b,c,d,e}(f).name))==0
%====================================== loop over individual slices
        
            % load processed data
            load(strcat(fpath_processed,slices{a,b,c,d,e}(f).name),...
                'pulset','baseS','baseD','indD4','baseIndex','fs','tpre')

            %% population spikes of baseline traces
            %==================================================================
            % manually find spikes
            t1 = pulset + .001*fs;                                      % start time for plot pop spikes (samples)
            t2 = pulset + .01*fs;                                       % end time for plotting pop spikes (samples)
            figure(2)

            % plot pre-induction spikes
            subplot(1,2,1)
            plot(t1:t2,baseS(t1:t2,baseIndex(1:tpre))...
                -ones(length(t1:t2),1)*baseS(t1,baseIndex(1:tpre)));      % plot all pop spikes pre-induction
            title({'Select times of 3 peak times', 'for population spike'})
            ylim([-.005 .005])
            [peakt,peakv] = ginput(3);                                  % select 3 peaks
            spike = zeros(length(baseIndex),1);                 
            for j = 1:20                                               % loop over pre-induction blocks
                fitx = [ones(2,1) round(peakt([1 3]))];                 % x values to fit a line between the positive peaks
                fity = baseS([round(peakt(1)) round(peakt(3))],...          
                    baseIndex(j));                                      % peak values to fit line
                fit = fitx\fity;                                        % fit line between the two peaks
                top = fit(1) + fit(2)*peakt(2);                         % extrapolate from negative peak to fitted line
                bottom = baseS(round(peakt(2)),baseIndex(j));           % value at negative peak
                spike(j) = abs(top-bottom);                             % pop spike amplitude
            end

            % repeat for pop spikes after induction
            % plot post-induction spikes
            subplot(1,2,2)
            plot(t1:t2,baseS(t1:t2,baseIndex(tpre+1:end))...
                -ones(length(t1:t2),1)*baseS(t1,baseIndex(tpre+1:end)));
            title('Select times of 3 consecutive peaks for population spike')
            [peakt,peakv] = ginput(3);
            for j = 21:length(baseIndex);                               
                fitx = [ones(2,1) round(peakt([1 3]))];
                fity = baseS([round(peakt(1)) round(peakt(3))],baseIndex(j));
                fit = fitx\fity;
                top = fit(1) + fit(2)*peakt(2);
                bottom = baseS(round(peakt(2)),baseIndex(j));
                spike(j) = abs(top-bottom);
            end
            % normalized population spike amplitude
            spikeN = spike/mean(spike(1:tpre));

            %% take slopes of baseline traces
            %==================================================================
            % manually take slopes
            figure(1)
            t1 = pulset + .001*fs;                                      % time window to plot responses (pulset is the time of bipolar pulse in samples)
            t2 = pulset + .01*fs;
            plot(t1:t2,baseD(t1:t2,baseIndex)-ones(length(t1:t2)...
                ,1)*mean(baseD(t1:t2,baseIndex),1));                % plot all baseline traces
            [tslope,yslope] = ginput(2);                            % select slope start and end time points
            slopesD = zeros(size(baseD,2),1);
            for j = 1:size(baseD,2)                             % store all slopes 
                slopesD(j) = mean( diff( baseD(...
                    round(tslope(1)):round(tslope(2)),j) ,1) )*fs;
            end
            
            tslope = tslope-pulset;
            
            %% slopes during tetanus
            %==============================================================
            indSlopes = zeros(size(indD4,2),size(indD4,3));
            for j = 1:size(indD4,3)
                for k = 1:size(indD4,2)
                    indSlopes(k,j) = mean( diff( indD4(...
                        round(tslope(1)):round(tslope(2)),k,j) ,1) )*fs; 
                end
            end
            indSlopes = indSlopes(:);

            %% store slopes variable
            %==================================================================
            slopes{a,b,c,d,e}(f).slopes = slopesD;
            slopes{a,b,c,d,e}(f).spike = spike;
            slopes{a,b,c,d,e}(f).indSlopes = indSlopes;
            slopes{a,b,c,d,e}(f).baseIndex = baseIndex;
            slopes{a,b,c,d,e}(f).tslope = tslope;
            slopes{a,b,c,d,e}(f).slices = slices{a,b,c,d,e}(f);
            slopes{a,b,c,d,e}(f).name = slices{a,b,c,d,e}(f).name;

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
save(strcat(fpath_variables,'slopes.mat'),'slopes')

