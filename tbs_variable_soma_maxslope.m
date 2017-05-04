% analyze spiking during induction
clear all; close all; clc

%% file paths and directories
fpath = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';
% fpath = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';

fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\Processed Variables\';

direct = dir(strcat(fpath,'*.mat*')); % processed matlab files

slice_remove = ['20170117_1_TBS_control_0Vm_apical_none_rig1.mat',...
    '20161019_1_TBS_anodal_20Vm_apical_none_rig1.mat'];

%% define conditions
induction = {'_TBS'}; % plasticity induction protocol
stim = {'_control';'_cathodal';'_anodal'};% DCS stimulation conditions
intensity = [0,5,20]; % stimulation intensity in V/m
position = {'_apical';'_basal';'_perforant'}; % recording location in slice
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
indTime = 7*fs;                                                             % total time of induction block(samples)
tbsOn = 2*fs;                                                               % tbs start time in s
tbsOff = 5*fs;                                                              % tbs end time in s
tBurst = .04*fs;                                                            % duration of each burst
nBurst = 15;                                                                % number of burst
nPulse = 4;                                                                 % number of pulses in a burst
rBurst = 100;                                                               % burst frequency
bipolarWidth = 2.5e-3*fs;                                                       % width of bipolar stimulus
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
         'HalfPowerFrequency1',200,'HalfPowerFrequency2',600, ...
         'SampleRate',fs);
     
%% arrange file names by condition
slices = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slice = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
cutoff(apical) = 20170115;
cutoff(basal) = 20170301;
cutoff(perforant) = 0;
for a = 1:length(induction)
    for b = 1:length(stim)
        for c = 1:length(intensity)
            for d = 1:length(position)
                for e = 1:length(drug)
                    slice{a,b,c,d,e} = dir(strcat(fpath,'*',...
                       induction{a},stim{b},strcat('_',num2str(intensity(c)),'Vm'),...
                       position{d},drug{e},'*.mat'));
                    
                    % choose slices
                    g = 0;
                    for f = 1:length(slice{a,b,c,d,e})
                        if isempty(strfind(slice_remove,slice{a,b,c,d,e}(f).name))==0 % check for bad slices (listed at beginning)
                            continue
                        elseif cutoff(d)<str2double(slice{a,b,c,d,e}(f).name(1:8)) % check for cutoff date
                            g = g+1;
                            slices{a,b,c,d,e}(g) = slice{a,b,c,d,e}(f);
                        end
                    end
                end
            end
        end
    end
end

%% preallocate
slopes = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesN = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesEnd = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indAll = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indSall = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indSfilt = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indSfilt_diff2 = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indSfilt_max = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indSfilt_maxN = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
baseSfilt = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
baseSfilt2 = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
baseSfilt_diff = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
baseSfilt_max = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
baseSfilt_maxN = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
baseSfilt_max_mean = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
indSfilt_diff = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
heights = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
dates = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
hemis = cell(length(induction),length(stim),length(intensity),length(position),length(drug));

% loop over conditions
for d = 1:length(position)
    figure(d);hold on
    for a = 1:length(induction)
        for b = 1:length(stim)
            for c = 1:length(intensity)
                for e = 1:length(drug)
                    % check for slices in current condition
                    if isempty(slices{a,b,c,d,e}) == 0
% preallocate
indAll{a,b,c,d,e} = zeros(indTime,length(slices{a,b,c,d,e}));
indSall{a,b,c,d,e} = zeros(indTime,length(slices{a,b,c,d,e}));
slopesN{a,b,c,d,e} = zeros(tbase+tpost,length(slices{a,b,c,d,e}));
slopesEnd{a,b,c,d,e} = zeros(length(slices{a,b,c,d,e}),1);
baseSfilt2{a,b,c,d,e} = cell(length(slices{a,b,c,d,e}),1);

% loop over slices
for f = 1:length(slices{a,b,c,d,e})
    %% load raw data
    load(strcat(fpath,slices{a,b,c,d,e}(f).name),...
        'indBlock','spike','spikeN',...
        'height','hemi','age','current','date','indS4',...
        'indS','fs','baseS','pulset','indD','indD4','baseD','slopesD');
    
    % analyze slopes
    if length(slopesD)>=indBlock(1)+tpost-1 % check that there are enough blocks
        % store slopes (blocks x slices)
        slopes{a,b,c,d,e}(:,f) = slopesD([indBlock(1)-tbase:indBlock(1)-1,indBlock(1):indBlock(1)+tpost-1]);
    end
    % full dendritic induction traces (time x slices)
    indAll{a,b,c,d,e}(:,f) = indD;  
    % cropped baseline dendritic traces, remove bipolar artifact (time x blocks)
    base = baseD((pulset+bipolarWidth):pulset+(fs/rBurst),:);   
    
    % normalize slopes
    slopesN{a,b,c,d,e}(:,f) = slopes{a,b,c,d,e}(:,f)/...
        mean(slopes{a,b,c,d,e}(1:tbase,f)); % (blocks x slices)

    % slope after 60 min for each slice
    slopesEnd{a,b,c,d,e}(f) = mean(slopesN{a,b,c,d,e}(end-9:end,f)); % (slices x 1)
    
    % height, date, hemisphere
    heights{a,b,c,d,e}(f) = height;
    dates{a,b,c,d,e}(f) = date;
    hemis{a,b,c,d,e}(f) = hemi;
    name{a,b,c,d,e}{f} = slices{a,b,c,d,e}(f).name;
    
    % somatic recording during baseline
    % store full baseline traces before cropping (time x blocks)
    baseSall = baseS; 
     % cropped baseline traces, remove bipolar artifact (time x slices)
    baseS = baseS((pulset+bipolarWidth):pulset+(fs/rBurst),:);

    % somatic recording during induction
    % full somatic recording during induction (time x slices)
    indSall{a,b,c,d,e}(:,f) = indS; 

    % high pass filter somatic recordings during baseline 
    baseSfilt21 = filtfilt(dfiltHigh2,baseSall); % full baseline recording filtered (time x blocks)
    baseSfilt2{a,b,c,d,e}{f}(:,:) = baseSfilt21((pulset+bipolarWidth):pulset+(fs/rBurst),:); % cropped filtered baseline recordings {slices}(time x blocks), all blocks
    baseSfilt{a,b,c,d,e}(:,:,f) = baseSfilt2{a,b,c,d,e}{f}(:,indBlock(1)-20:indBlock(1)-1); % keep 20 pre-induction blocks and store in matrix (time x blocks x slices)

    % high pass somatic recordings during induction
    indSfilt{a,b,c,d,e}(:,f) = filtfilt(dfiltHigh2,...
        indSall{a,b,c,d,e}(:,f)); % full somatic induction recordings filtered (time x slices)

end
% organize so each column is the time series of
% a single burst, 3rd dimension is slice number
indSfilt{a,b,c,d,e} = indSfilt{a,b,c,d,e}(tbsOn+1:tbsOff,:); % filtered signal during TBS (time x slices)
indSfilt{a,b,c,d,e} = reshape(indSfilt{a,b,c,d,e},[],...
    nBurst,length(slices{a,b,c,d,e})); % (time x burst x slices)

% remove time between bursts
indSfilt{a,b,c,d,e} = indSfilt{a,b,c,d,e}(1:tBurst,:,:); % (time x burst x slice)

% take derivatives
baseSfilt_diff{a,b,c,d,e} = diff(baseSfilt{a,b,c,d,e},1,1); % first derivate of baseline filtered signal (time x blocks x slices)
indSfilt_diff{a,b,c,d,e} = diff(indSfilt{a,b,c,d,e},1,1); % first derivative of cropped filtered somatic signal during induction (time x burst x slice) 

% find max slope during each pulse
[baseSfilt_max{a,b,c,d,e},baseSfilt_max_i{a,b,c,d,e}] = max(-baseSfilt_diff{a,b,c,d,e},[],1);
baseSfilt_max{a,b,c,d,e} = squeeze(baseSfilt_max{a,b,c,d,e});
baseSfilt_max_i{a,b,c,d,e} = squeeze(baseSfilt_max_i{a,b,c,d,e});
baseSfilt_max_mean{a,b,c,d,e} = mean(baseSfilt_max{a,b,c,d,e},1);
baseSfilt_maxN{a,b,c,d,e} = baseSfilt_max{a,b,c,d,e}./(ones(size(baseSfilt_max{a,b,c,d,e},1),1)*baseSfilt_max_mean{a,b,c,d,e});

% reshape induction indSfilt to fidn max during
% each pulse
indSfilt_diff2{a,b,c,d,e} = reshape(indSfilt_diff{a,b,c,d,e}(1:396,:,:),[],nPulse,nBurst,length(slices{a,b,c,d,e}));
[indSfilt_max{a,b,c,d,e},indSfilt_max_i{a,b,c,d,e}] = (max(-indSfilt_diff2{a,b,c,d,e}(bipolarWidth:end-bipolarWidth,:,:,:),[],1)); % pulses x bursts x slices
indSfilt_max{a,b,c,d,e} = squeeze(indSfilt_max{a,b,c,d,e});
indSfilt_max_i{a,b,c,d,e} = squeeze(indSfilt_max_i{a,b,c,d,e})+bipolarWidth;

% normalize
for g = 1:length(slices{a,b,c,d,e})
    indSfilt_maxN{a,b,c,d,e}(:,:,g) = indSfilt_max{a,b,c,d,e}(:,:,g)/baseSfilt_max_mean{a,b,c,d,e}(g);
end
indSfilt_maxN{a,b,c,d,e} = reshape(indSfilt_maxN{a,b,c,d,e},nPulse*nBurst,[]);
indSfilt_max2{a,b,c,d,e} = reshape(indSfilt_max{a,b,c,d,e},nPulse*nBurst,[]);
                    end
                end
            end
        end
    end
end

save(strcat(fpath_variables,'somatic_maxslope.mat'),'indSfilt','baseSfilt',...
    'baseSfilt_max','baseSfilt_max_i','indSfilt_maxN','indSfilt_max_i','slices')
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
figure;plot(mean(indSfilt_maxN{1,3,3,2,1},1),slopesEnd{1,3,3,2,1},'.r','MarkerSize',10)
hold on;plot(mean(indSfilt_maxN{1,1,1,2,1},1),slopesEnd{1,1,1,2,1},'.k','MarkerSize',10)
hold on;plot(mean(indSfilt_maxN{1,2,3,2,1},1),slopesEnd{1,2,3,2,1},'.b','MarkerSize',10)
xlabel('Normalize peak slope in somatic recording')
ylabel('LTP')