% analyze fEPSP area during induction
clear all; close all;clc

%% file paths and directories
fpath = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';
fpath_processed = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\Processed Variables\';
% fpath = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';
direct = dir(strcat(fpath,'*.mat*')); % processed matlab files

%% define conditions
induction = {'_TBS'};                                                       % plasticity induction protocol
stim = {'_control';'_cathodal';'_anodal'};                                  % DCS stimulation conditions
intensity = [0,5,20]; % stimulation intensity in V/m
position = {'_apical';'_basal';'_perforant'}; % recording location in slice
locs = {'A','B','P','S'}; % loc 
drug = {'_none';'_mk801'};

control = find(strcmp(stim,'_control'));
cathodal = find(strcmp(stim,'_cathodal'));
anodal = find(strcmp(stim,'_anodal'));
soma = find(strcmp(position,'_soma'));
apical = find(strcmp(position,'_apical'));
basal = find(strcmp(position,'_basal'));
perforant = find(strcmp(position,'_perforant'));

stimcolor = {[0,0,0],[0,0,1],[1,0,0]};% figure colors for each stim condition
stimcolor_shade = {[.5 .5 .5],[.7 .7 1],[1 .7 .7]};

%% time conditions
% baseline responses
tbase = 20;
tpost = 60;
t = [1:(tbase+tpost)]';

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

%% design filters
% high pass
Fstop = 0.1;
Fpass = 1;
Astop = 30;
Apass = 0.5;
dfiltHigh1 = designfilt('highpassiir','StopbandFrequency',Fstop ,...
  'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
  'PassbandRipple',Apass,'SampleRate',fs);
dfiltHigh2 = designfilt('bandpassiir','FilterOrder',20,...
         'HalfPowerFrequency1',200,'HalfPowerFrequency2',1000, ...
         'SampleRate',fs);

%% arrange file names by condition
slices = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slice = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
cutoff(apical) = 0;
cutoff(basal) = 0;
cutoff(perforant) = 0;
for a = 1:length(induction)
    for b = 1:length(stim)
        for c = 1:length(intensity)
            for d = 1:length(position)
                for e = 1:length(drug)
                        slice{a,b,c,d,e} = dir(strcat(fpath,'*',...
                           induction{a},stim{b},strcat('_',num2str(intensity(c)),'Vm'),...
                           position{d},drug{e},'*.mat'));
                       % choose slices after a certain date
                       g = 0;
                       for f = 1:length(slice{a,b,c,d,e})
                           if cutoff(d)<str2double(slice{a,b,c,d,e}(f).name(1:8))
                               g = g+1;
                               slices{a,b,c,d,e}(g) = slice{a,b,c,d,e}(f);
                           end
                       end
                end
            end
        end
    end
end

%% load each data file and store processed data
baseArea = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indFilt1 = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indFilt2 = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indAll = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesN = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopes = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesEnd = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
heights = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
dates = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
hemis = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
burstArea = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
burstAreaN = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
burstAreaNmean = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
baseFilt1 = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
baseFilt2 = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
X = [];
for d = 1:length(position)
    figure(d);hold on
    for a = 1:length(induction)
        for b = 1:length(stim)
            for c = 1:length(intensity)
                for e = 1:length(drug)
                    if isempty(slices{a,b,c,d,e}) == 0
                        baseArea{a,b,c,d,e} = zeros(1,length(slices{a,b,c,d,e}));
                        for f = 1:length(slices{a,b,c,d,e});
                            
                            %% load raw data
                            load(strcat(fpath,slices{a,b,c,d,e}(f).name),...
                                'indBlock','spike','spikeN',...
                                'height','hemi','age','current','date','indS4',...
                                'indS','fs','baseS','pulset','indD','indD4','baseD','slopesD',...
                                'baseD_filt_high_1','indD_filt_high_1');
                            
                            % check that there are enough recording blocks
                            if length(slopesD)>=indBlock(1)+tpost-1
                                slopes{a,b,c,d,e}(:,f) = slopesD([indBlock(1)-tbase:indBlock(1)-1,indBlock(1):indBlock(1)+tpost-1]); % slopes x slices
                            end
                            indAll{a,b,c,d,e}(:,f) = indD;  % induction traces (time x slices)
                            base = baseD((pulset+bipolarWidth):pulset+(fs/rBurst),:); % cropped baseline traces (time x blocks)
                            baseAll = baseD; % full baseline traces (time x blocks)
                            
                            % slopes and heights
                            % normalizes slopes
                            slopesN{a,b,c,d,e}(:,f) = slopes{a,b,c,d,e}(:,f)/...
                                mean(slopes{a,b,c,d,e}(1:tbase,f)); % (blocks x slices)
                        
                            % slope after 60 min for each slice
                            slopesEnd{a,b,c,d,e}(f) = mean(slopesN{a,b,c,d,e}(end-9:end,f)); % (slices x 1)
                            
                            % height, date, hemisphere
                            heights{a,b,c,d,e}(f) = height;% (slices x 1)
                            dates{a,b,c,d,e}(f) = date;% (slices x 1)
                            hemis{a,b,c,d,e}(f) = hemi;% (slices x 1)
                            
                            % area under baseline fEPSPs
                            % Normalize first sample to zero
                            base1 = ones(size(base,1),1)*base(1,:); % each column is the value of the first sample for that block
                            base = base-base1; % normalized and cropped baseline trace (time x blocks)
                            
                            % filter time series during baseline to remove
                            % drift
                            baseFilt1{a,b,c,d,e}{f} = baseD_filt_high_1;
                            
                            % cropped trace to region for taking area
                            baseFilt2{a,b,c,d,e}{f} = baseFilt1{a,b,c,d,e}{f}((pulset+bipolarWidth):pulset+(fs/rBurst),:);
                            
                            % Area under fEPSP
                            baseArea{a,b,c,d,e}(f) = mean(sum(abs(baseFilt2{a,b,c,d,e}{f}(:,...
                                indBlock(1)-20:indBlock(1)-1)),1)); % average area under basleine fEPSP
                            
                            % filter time series recordings during induction to remove drift
                            indFilt1{a,b,c,d,e}(:,f) = indD_filt_high_1; 
                        end
                        
                        % Make each column the time series of a single burst
                        % 3rd dimension is slice number
                        indFilt1{a,b,c,d,e} = indFilt1{a,b,c,d,e}(tbsOn+1:tbsOff,:); % traces during TBS (time x slices)
                        indFilt1{a,b,c,d,e} = reshape(indFilt1{a,b,c,d,e},[],...
                            nBurst,length(slices{a,b,c,d,e})); % reshape (time x burst x slice)
                        
                        % remove time between bursts
                        indFilt1{a,b,c,d,e} = indFilt1{a,b,c,d,e}(1:tBurst,:,:); % (time x burst x slice)
                        
                        % remove bipolar pulse
                        indFilt1{a,b,c,d,e}(bipolarT,:,:) = []; % (time x burst x slice)
                        
                        % Normalize first sample to zero
                        indFilt1{a,b,c,d,e} = indFilt1{a,b,c,d,e}-...
                            indFilt1{a,b,c,d,e}(1,:,:); %(time x burst x slice)
                        
                        indFilt2{a,b,c,d,e} = reshape(indFilt1{a,b,c,d,e},[],size(indFilt1{a,b,c,d,e},3));
                        
                        % integrate under each burst envelope
                        burstArea{a,b,c,d,e} = permute(sum(abs(indFilt1{a,b,c,d,e}),1),[2,3,1]); % (burst x slices)
                        
                        % Normalize to area during baseline
                        burstAreaN{a,b,c,d,e} = burstArea{a,b,c,d,e}./...
                            (ones(nBurst,1)*baseArea{a,b,c,d,e});
                        
                        % Average burst area
                        burstAreaNmean{a,b,c,d,e} = mean(burstAreaN{a,b,c,d,e},1);
                        
                        % load data into a single matrix for
                        % regression/correlation
                        X = [X;d*ones(length(slopesEnd{a,b,c,d,e}),1),...
                                b*ones(length(slopesEnd{a,b,c,d,e}),1),...
                                burstAreaNmean{a,b,c,d,e}',heights{a,b,c,d,e}',...
                                slopesEnd{a,b,c,d,e}'];
                        

                        
                        %% figures
                        figure(d)
                        plot(burstAreaNmean{a,b,c,d,e},slopesEnd{a,b,c,d,e},'.','Color',stimcolor{b},'MarkerSize',15)
                        xlabel('Normalized area under theta burst dendritic recording')
                        ylabel('Synaptic strength')
                        title(position{d}(2:end))
                    end
                end
            end
        end
    end
    % correlation for each recording location
    loc_ind = X(:,1) == d;
    [R{d},P{d}] = corrcoef(X(loc_ind,[3,5]));
end

%% figures
burst_plot = 3; % which burst to plot
for d = 1:length(position)
    figure;hold on
    for a = 1:length(induction)
        for b = 1:length(stim)
            for c = 1:length(intensity)
                for e = 1:length(drug)
                    if isempty(slices{a,b,c,d,e}) == 0
                        plot(squeeze(mean(indFilt1{a,b,c,d,e}(:,burst_plot,:),3)),'Color',stimcolor{b},'LineWidth',2)
%                         plot(squeeze(mean(indFilt1{a,b,c,d,e}(:,burst_plot,:),3)) + squeeze(std(indFilt1{a,b,c,d,e}(:,burst_plot,:),0,3)/sqrt(size(indFilt1{a,b,c,d,e},3))),'Color',stimcolor_shade{b})
%                         plot(squeeze(mean(indFilt1{a,b,c,d,e}(:,burst_plot,:),3)) - squeeze(std(indFilt1{a,b,c,d,e}(:,burst_plot,:),0,3)/sqrt(size(indFilt1{a,b,c,d,e},3))),'Color',stimcolor_shade{b})
                    end
                end
            end
        end
    end
end
save(strcat(fpath_processed,'fepsp_area.mat'),...
    'indFilt1','indFilt2','baseFilt1','baseFilt2','baseArea','burstAreaN','slices')