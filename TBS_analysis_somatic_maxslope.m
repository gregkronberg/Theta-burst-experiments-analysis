% analyze spiking during induction
clear all; close all; clc

%% file paths and directories
% fpath = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';
fpath = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';
direct = dir(strcat(fpath,'*.mat*')); % processed matlab files

slice_remove = ['20170117_1_TBS_control_0Vm_apical_none_rig1.mat',...
    '20161019_1_TBS_anodal_20Vm_apical_none_rig1.mat'];
%% define conditions
induction = {'_TBS'}; % plasticity induction protocol
stim = {'_control';'_cathodal';'_anodal'};% DCS stimulation conditions
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
         'HalfPowerFrequency1',100,'HalfPowerFrequency2',800, ...
         'SampleRate',fs);
     
%% arrange file names by condition
slices = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
for a = 1:length(induction)
    for b = 1:length(stim)
        for c = 1:length(intensity)
            for d = 1:length(position)
                for e = 1:length(drug)
                        slices{a,b,c,d,e} = dir(strcat(fpath,'*',...
                           induction{a},stim{b},strcat('_',num2str(intensity(c)),'Vm'),...
                           position{d},drug{e},'*.mat'));
                end
            end
        end
    end
end
slice_remove_n = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
for d = 1:length(position)
    figure(d);hold on
    for a = 1:length(induction)
        for b = 1:length(stim)
            for c = 1:length(intensity)
                for e = 1:length(drug)
                    if isempty(slices{a,b,c,d,e}) == 0
                        baseVar{a,b,c,d,e} = zeros(1,length(slices{a,b,c,d,e}));
                        for f = 1:length(slices{a,b,c,d,e});
                            if isempty(strfind(slice_remove,slices{a,b,c,d,e}(f).name))==0
                                slice_remove_n{a,b,c,d,e} = [slice_remove_n{a,b,c,d,e};f];
                                continue
                            end
                            %% load raw data
                            load(strcat(fpath,slices{a,b,c,d,e}(f).name),...
                                'indBlock','spike','spikeN',...
                                'height','hemi','age','current','date','indS4',...
                                'indS','fs','baseS','pulset','indD','indD4','baseD','slopesD');
                            if length(slopesD)>=indBlock(1)+tpost-1
                                slopes{a,b,c,d,e}(:,f) = slopesD([indBlock(1)-tbase:indBlock(1)-1,indBlock(1):indBlock(1)+tpost-1]);
                            end
                            indAll{a,b,c,d,e}(:,f) = indD;
                            base = baseD((pulset+bipolarWidth):pulset+(fs/rBurst),:);
                            baseAll = baseD;
                            
                            % slopes and heights
                            % normalizes slopes
                            slopesN{a,b,c,d,e}(:,f) = slopes{a,b,c,d,e}(:,f)/...
                            mean(slopes{a,b,c,d,e}(1:tbase,f));
                        
                            % slope after 60 min for each slice
                            slopesEnd{a,b,c,d,e}(f) = mean(slopesN{a,b,c,d,e}(end-9:end,f));
                            
                            % height, date, hemisphere
                            heights{a,b,c,d,e}(f) = height;
                            dates{a,b,c,d,e}(f) = date;
                            hemis{a,b,c,d,e}(f) = hemi;
                            
                             % load somatic recording during baseline
                            baseSall = baseS;
                            baseS = baseS((pulset+bipolarWidth):pulset+(fs/rBurst),:);
                            
                             % load somatic recording during induction
                            indSall{a,b,c,d,e}(:,f) = indS;
                            
                            % high pass filter somatic recordings during baseline 
                            baseSfilt21 = filtfilt(dfiltHigh2,baseSall);
                            baseSfilt2{a,b,c,d,e}{f}(:,:) = baseSfilt21((pulset+bipolarWidth):pulset+(fs/rBurst),:);
                            baseSfilt{a,b,c,d,e}(:,:,f) = baseSfilt2{a,b,c,d,e}{f}(:,indBlock(1)-20:indBlock(1)-1);
                            
                            % high pass somatic recordings during induction
                            indSfilt{a,b,c,d,e}(:,f) = filtfilt(dfiltHigh2,...
                                indSall{a,b,c,d,e}(:,f));
                            
                        end
                        if isempty(slice_remove_n{a,b,c,d,e})~=1;
                                slopesEnd{a,b,c,d,e}(slice_remove_n{a,b,c,d,e}) = [];
                        end
                        % organize so each column is the time series of
                        % a single burst, 3rd dimension is slice number
                        indSfilt{a,b,c,d,e} = indSfilt{a,b,c,d,e}(tbsOn+1:tbsOff,:);
                        indSfilt{a,b,c,d,e} = reshape(indSfilt{a,b,c,d,e},[],...
                            nBurst,length(slices{a,b,c,d,e}));
                        
                        % remove time between bursts
                        indSfilt{a,b,c,d,e} = indSfilt{a,b,c,d,e}(1:tBurst,:,:);
                        
                        % take derivatives
                        baseSfilt_diff{a,b,c,d,e} = diff(baseSfilt{a,b,c,d,e},1,1);
                        indSfilt_diff{a,b,c,d,e} = diff(indSfilt{a,b,c,d,e},1,1);
                        
                        % find max slope during each pulse
                        baseSfilt_max{a,b,c,d,e} = squeeze(max(-baseSfilt_diff{a,b,c,d,e},[],1));
                        baseSfilt_max_mean{a,b,c,d,e} = mean(baseSfilt_max{a,b,c,d,e},1);
                        baseSfilt_maxN{a,b,c,d,e} = baseSfilt_max{a,b,c,d,e}./(ones(size(baseSfilt_max{a,b,c,d,e},1),1)*baseSfilt_max_mean{a,b,c,d,e});
                        
                        % reshape induction indSfilt to fidn max during
                        % each pulse
                        indSfilt_diff2{a,b,c,d,e} = reshape(indSfilt_diff{a,b,c,d,e}(1:396,:,:),[],nPulse,nBurst,length(slices{a,b,c,d,e}));
                        indSfilt_max{a,b,c,d,e} = squeeze(max(-indSfilt_diff2{a,b,c,d,e}(bipolarWidth:end-bipolarWidth,:,:,:),[],1)); % pulses x bursts x slices
                        
                        % normalize
                        for g = 1:length(slices{a,b,c,d,e})
                            indSfilt_maxN{a,b,c,d,e}(:,:,g) = indSfilt_max{a,b,c,d,e}(:,:,g)/baseSfilt_max_mean{a,b,c,d,e}(g);
                        end
                        indSfilt_maxN{a,b,c,d,e} = reshape(indSfilt_maxN{a,b,c,d,e},nPulse*nBurst,[]);
                        indSfilt_max{a,b,c,d,e} = reshape(indSfilt_max{a,b,c,d,e},nPulse*nBurst,[]);
                        if isempty(slice_remove_n{a,b,c,d,e})==0;
                            indSfilt_maxN{a,b,c,d,e}(:,slice_remove_n{a,b,c,d,e}) = [];
                        end
                    end
                end
            end
        end
    end
end

%% figures
figure;plot(mean(indSfilt_maxN{1,3,3,1,1},2),'r')
hold on;plot(mean(indSfilt_maxN{1,1,1,1,1},2),'k')
plot(mean(indSfilt_maxN{1,2,3,1,1},2),'b')
title('apical')
xlabel('TBS pulse number')
ylabel('Normalize peak slope in somatic recording')

figure;plot(mean(indSfilt_maxN{1,3,3,2,1},2),'r')
hold on;plot(mean(indSfilt_maxN{1,1,1,2,1},2),'k')
plot(mean(indSfilt_maxN{1,2,3,2,1},2),'b')
title('basal')
xlabel('TBS pulse number')
ylabel('Normalize peak slope in somatic recording')

figure;plot(mean(indSfilt_maxN{1,3,3,3,1},2),'r')
hold on;plot(mean(indSfilt_maxN{1,1,1,3,1},2),'k')
plot(mean(indSfilt_maxN{1,2,3,3,1},2),'b')
title('perforant')
xlabel('TBS pulse number')
ylabel('Normalize peak slope in somatic recording')

figure;plot(indSfilt_maxN{1,3,3,1,1},'r')
hold on;plot((indSfilt_maxN{1,1,1,1,1}),'k')
plot((indSfilt_maxN{1,2,3,1,1}),'b')
title('apical')
xlabel('TBS pulse number')
ylabel('Normalize peak slope in somatic recording')

figure;plot((indSfilt_maxN{1,3,3,2,1}),'.r')
hold on;plot((indSfilt_maxN{1,1,1,2,1}),'.k')
plot((indSfilt_maxN{1,2,3,2,1}),'.b')
title('basal')
xlabel('TBS pulse number')
ylabel('Normalize peak slope in somatic recording')

figure;plot((indSfilt_maxN{1,3,3,3,1}),'r')
hold on;plot((indSfilt_maxN{1,1,1,3,1}),'k')
plot((indSfilt_maxN{1,2,3,3,1}),'b')
title('perforant')
xlabel('TBS pulse number')
ylabel('Normalize peak slope in somatic recording')

%% 
figure;plot(mean(indSfilt_maxN{1,3,3,1,1},1),slopesEnd{1,3,3,1,1},'.')
xlabel('Normalize peak slope in somatic recording')
ylabel('LTP')