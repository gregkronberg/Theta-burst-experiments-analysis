%% Synaptic plasticity experiments preprocessing

% inputs are matlab files exported from labchart

% outputs are features of each slice (e.g. animal age, brain hemisphere,
% stimulation conditions), time series data organized into matrices (time x
% recording block), and fEPSP slopes over time and population spikes over
% time

clear all; close all; clc

%% file paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% desktop
fpathR = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
fpathP = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
% laptop
% fpathR = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\';
% fpathP = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';

%% list new slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[slicesN,slicesR,slicesP] = folder_diff(fpathR,fpathP,'.mat');

%% process new slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(slicesN)
            
    load(strcat(fpathR,slicesN{i}));    % load new slice
    fs = mode(samplerate(:));           % sampling rate
    l = dataend-datastart+1;            % length of each block in samples

    %% extract info from comments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % induction blocks and type
    numInd = 4;                                                     % number of possible induction periods
    indBlock = zeros(numInd,1);                                     % recording block number for plasticity induction
    indType = cell(numInd,1);                                      % type of plasticity induction, e.g. TBS 
    indBlock(1) = comment_find('induction_','block',com,comtext);   % first induction
    indType{1} = comment_find('induction','value',com,comtext);
    indBlock(2) = comment_find('induction2','block',com,comtext);   % second induction
    indType{2} = comment_find('induction2','value',com,comtext);
    indBlock(3) = comment_find('induction3','block',com,comtext);   % third induction
    indType{3} = comment_find('induction3','value',com,comtext);
    indBlock(4) = comment_find('induction4','block',com,comtext);   % fourth induction
    indType{4} = comment_find('induction4','value',com,comtext);
    baseIndex = [indBlock(1)-20:indBlock(1)-1,indBlock(1)+1:indBlock(1)+60];    % index of baseline responses (pre/post induction)
    
    % channel number for each recording location
    soma = comment_find('Soma','channel',com,comtext);              % 
    apical = comment_find('Apical','channel',com,comtext);
    basal = comment_find('Basal','channel',com,comtext);
    perforant = comment_find('Perforant','channel',com,comtext);
    
    % physiological parameters
    age = comment_find('age','value',com,comtext);
    hemi = comment_find('hemi','value',com,comtext);
    height = comment_find('height','value',com,comtext);
    
    % stimulation parameters
    dcs = 
    current = comment_find('current','value',com,comtext);
    isolator = comment_find('isolator','value',com,comtext);
    path2 = comment_find('2path','value',com,comtext);

    % date and slice number
    date = str2double(slicesN{i}(1:8));                             % date of recording
    slice_num = str2double(slicesN{i}(10));                        % slice number used in chamber
    
    % timing of bipolar pulse
    if date>=20161013                                              % changed the timing of probe pulses, slices that have a rig number have the new pulse time
        pulset = .5*fs;                                             % start time of bipolar pulse (samples)
    else pulset = 0*fs ;                                             
    end
    pulset_path2 = 30.5*fs;                                         % start time of bipolar pulse 

    %% organize baseline time series traces into matrices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % base_ is a matrix where each column is trace from baseline
    baseS = zeros(mode(l(:)),size(datastart,2));                    % baseline somatic traces (time x blocks)
    baseD = zeros(mode(l(:)),size(datastart,2));                    % baseline dendritic traces, first path (time x blocks)
    baseD_path2 = zeros(mode(l(:)),size(datastart,2));              % baseline dendritic traces, second path (time x blocks)
    
    % loop over blocks
    for j = 1:size(datastart,2)                                    
        % skip induction blocks
        if l(:,j)==mode(l,2)                                       
            % store probe traces
            % soma
            if soma ~= 0
                baseS(:,j) =  data(datastart(soma,j):dataend(soma,j))';
            end
            % apical
            if apical ~= 0
                baseD(:,j) =  data(datastart(apical,j):dataend(apical,j))';
            % basal
            elseif basal ~= 0
                baseD(:,j) =  data(datastart(basal,j):dataend(basal,j))';
            % perforant    
            elseif perforant ~= 0
                baseD(:,j) =  data(datastart(perforant,j):dataend(perforant,j))';    
            end
        end
    end

    %% take population spikes of baseline traces
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % manually find spikes
    t1 = pulset + .001*fs;                                      % start time for plot pop spikes (samples)
    t2 = pulset + .01*fs;                                       % end time for plotting pop spikes (samples)
    figure(2)

    % plot pre-induction spikes
    subplot(1,2,1)
    plot(t1:t2,baseS(t1:t2,baseIndex(1:20))...
        -ones(length(t1:t2),1)*baseS(t1,baseIndex(1:20)));      % plot all pop spikes pre-induction
    title('Select times of 3 consecutive peaks for population spike')
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
    plot(t1:t2,baseS(t1:t2,baseIndex(21:end))...
        -ones(length(t1:t2),1)*baseS(t1,baseIndex(21:end)));
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
    spikeN = spike/mean(spike(1:20));

    %% take slopes of baseline traces
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    %% store induction block
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tbsOn = 2;                                                      % tbs start time in s
    tbsOff = 5;                                                     % tbs end time in s
    tBurst = .04;                                                   % duration of each burst
    nBurst = 15;                                                    % number of burst
    nPulse = 4;                                                     % number of pulses in a burst
    rBurst = 100;                                                   % burst frequency

    % soma traces during induction
    if soma ~= 0                                                   % if there is a somatic recording
        indS =  data(datastart(soma,indBlock(1)):...
            dataend(soma,indBlock(1)))';                               % extract data from induction block
        indS1 = indS(tbsOn*fs+1:tbsOff*fs);                         % take data only during bursts
        indS2 = reshape(indS1,[],nBurst);                           % organize into 2D matrix (time x bursts)
        indS3 = indS2(1:round(tBurst*fs),:);                        % remove time in between bursts so that length of each pulse is the same
        indS4 = reshape(indS3,[],nPulse,nBurst);                    % reshape to 3D matrix (time x pulses x bursts)
    end

    % repeat for dendritic traces during induction 
    % check for recording locations
    if apical ~= 0
        indD =  data(datastart(apical,indBlock(1)):dataend(apical,indBlock(1)))';
    elseif basal~= 0 
        indD =  data(datastart(basal,indBlock(1)):dataend(basal,indBlock(1)))';
    elseif perforant ~= 0 
        indD =  data(datastart(perforant,indBlock(1)):dataend(perforant,indBlock(1)))';
    end
    indD1 = indD(tbsOn*fs+1:tbsOff*fs);
    indD2 = reshape(indD1,[],nBurst);
    indD3 = indD2(1:round(tBurst*fs),:);
    indD4 = reshape(indD3,[],nPulse,nBurst);

    %% save workspace variables
    save(strcat(fpathP,slicesN{i}));
end