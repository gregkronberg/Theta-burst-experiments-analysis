clear all; close all; clc

%% file paths
fpathR = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\';
fpathP = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';
% fpathR = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\';
% fpathP = 'C:\Users\Greg Kronberg\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\';

%% directory listings of all files
directR = dir(strcat(fpathR,'*.mat*')); % raw matlab files
directP = dir(strcat(fpathP,'*.mat*')); % processed matlab files

%% cell arrays of files names
slicesR = cell(size(directR,1),1); %create empty cell array
for a = 1:size(directR,1)
    slicesR{a} = directR(a).name; % raw file names
end
slicesP = cell(size(directP,1),1); %create empty cell array
for a = 1:size(directP,1)
    slicesP{a} = directP(a).name; % processed file names
end

%% process new slices
if length(slicesR)~= length(slicesP)                                       % are there new slices?
    for i = 1:length(slicesR)                                              % loop over all slices
        if sum(strcmp(slicesP,slicesR{i}))==0                              % check if current slice is new
            
            load(strcat(fpathR,slicesR{i}));                                % load new slice
            fs = mode(samplerate(:));                                       % sampling rate
            l = dataend-datastart+1;                                        % length of each block in samples
            
            
            %% extract info from comments
            numInd = 4;
            indBlock = zeros(numInd,1);
            indBlock_shift = zeros(numInd,1);
            soma=0;apical=0;basal=0;perforant=0;age=0;hemi=0;height=0;current=0;        % preallocate
            z = 0;
            path2 = {};
            for j = 1:size(comtext,1)                                      % loop through all comments
                
                if isempty(strfind(comtext(j,:),'induction_')) == 0         % find plasticity induction block
                   indBlock(1) = com(com(:,5)==j,2);
                end
                if isempty(strfind(comtext(j,:),'induction2_')) == 0         % find plasticity induction block
                   indBlock(2) = com(com(:,5)==j,2);
                end
                if isempty(strfind(comtext(j,:),'induction3_')) == 0         % find plasticity induction block
                   indBlock(3) = com(com(:,5)==j,2);
                end
                if isempty(strfind(comtext(j,:),'induction4_')) == 0         % find plasticity induction block
                   indBlock(4) = com(com(:,5)==j,2);
                end
                
                if isempty(strfind(comtext(j,:),'Soma')) == 0          % find channel with somatic recording
                    soma = com(com(:,5)==j,1);
                end
                
                if isempty(strfind(comtext(j,:),'Apical')) == 0        % channel with apical recording
                    apical = com(com(:,5)==j,1);
                end
                
                if isempty(strfind(comtext(j,:),'Basal')) == 0         % channel with basal recording
                    basal = com(com(:,5)==j,1);
                end
                
                if isempty(strfind(comtext(j,:),'Perforant')) == 0         % channel with perforant recording
                    perforant = com(com(:,5)==j,1);
                end
                
                if isempty(strfind(comtext(j,:),'age')) == 0           % age of animal
                    ageI = strfind(comtext(j,:),'age');
                    ageL = length('age_');
                    age = str2double(comtext(j,ageI+ageL:ageI+ageL+1));
                end
                
                if isempty(strfind(comtext(j,:),'hemi')) == 0          % brain hemisphere
                    hemiI = strfind(comtext(j,:),'hemi');
                    hemiL = length('hemi_');
                    hemi = comtext(j,hemiI+hemiL);
                end
                if isempty(strfind(comtext(j,:),'height')) == 0        % position along dorsal ventral axis of hippocampus (x0.1 mm from IA from paxinos atlas)
                    heightI = strfind(comtext(j,:),'height');
                    heightL = length('height_');
                    height = str2double(comtext(j,heightI+heightL:heightI+heightL+1));
                end
                
                if isempty(strfind(comtext(j,:),'current')) == 0       % current injected for DCS (in uA)
                    currentI = strfind(comtext(j,:),'current');
                    currentL = length('current_');
                    current = str2double(comtext(j,currentI+currentL:end));
                end
                
                if isempty(strfind(comtext(j,:),'isolator')) == 0       % DCS isolator on throughout the experiment?
                    isolatorI = strfind(comtext(j,:),'isolator');
                    isolatorL = length('isolator_');
                    isolator = str2double(comtext(j,currentI+currentL:end));
                end
                
                if isempty(strfind(comtext(j,:),'2path_apical')) == 0         % find if there was a second pathway stimulated and which pathway it was
                    path2 = 'apical';
                elseif isempty(strfind(comtext(j,:),'2path_basal')) == 0
                    path2 = 'basal'
                elseif isempty(strfind(comtext(j,:),'2path_basal')) == 0
                    path2 = 'perforant'
                    indBlock(1) = com(com(:,5)==j,2);
                end
            end
            
            date = str2double(slicesR{i}(1:8));                             % date of recording
            slice_num = str2double(slicesR{i}(10));                        % slice number used in chamber
            if date>=20161013                                              % changed the timing of probe pulses, slices that have a rig number have the new pulse time
                pulset = .5*fs;                                             % start time of bipolar pulse (samples)
            else pulset = 0*fs ;                                             
            end
            pulset_path2 = 30.5*fs;                                         % start time of bipolar pulse 
            
            
            
            baseIndex = [indBlock(1)-20:indBlock(1)-1,indBlock(1)+1:indBlock(1)+60];    % index of baseline responses (pre/post induction)

            %% organize baseline traces into column vectors
            % base_ is a matrix where each column is trace from baseline
            baseS = zeros(mode(l(:)),size(datastart,2));                    % baseline somatic traces
            baseD = zeros(mode(l(:)),size(datastart,2));                    % baseline dendritic traces (first path)
            baseD_path2 = zeros(mode(l(:)),size(datastart,2));              % baseline dendritic traces (second path)
            
            for j = 1:size(datastart,2);                                    % loop over blocks
                if l(:,j)==mode(l,2);                                       % skip induction blocks
                    
                    % store probe traces
                    if soma ~= 0;
                        baseS(:,j) =  data(datastart(soma,j):dataend(soma,j))';
                    end
                    
                    if apical ~= 0;
                        baseD(:,j) =  data(datastart(apical,j):dataend(apical,j))';
                    elseif basal ~= 0
                        baseD(:,j) =  data(datastart(basal,j):dataend(basal,j))';
                    elseif perforant ~= 0
                        baseD(:,j) =  data(datastart(perforant,j):dataend(perforant,j))';    
                    end
                    
                    
                end
            end
            if size(baseD,2)<max(baseIndex)
                baseD = [baseD,baseD(:,end)];
            end
            if size(baseS,2)<max(baseIndex)
                baseS = [baseS,baseS(:,end)];
            end
            
            
            %% take population spikes of baseline traces
            % manually find spikes
            t1 = pulset + .001*fs;                                      % start time for plot pop spikes (samples)
            t2 = pulset + .01*fs;                                       % end time for plotting pop spikes (samples)
            figure(2)

            % plot pre-induction spikes
            subplot(1,2,1)
            plot(t1:t2,baseS(t1:t2,baseIndex(1:20))...
                -ones(length(t1:t2),1)*baseS(t1,baseIndex(1:20)));      % plot all pop spikes pre-induction
            [peakt,peakv] = ginput(3);                                  % select 3 peaks
            spike = zeros(length(baseIndex),1);                 
            for j = 1:20;                                               % loop over pre-induction blocks
                fitx = [ones(2,1) round(peakt([1 3]))];                 % x values to fit a line between the positive peaks
                fity = baseS([round(peakt(1)) round(peakt(3))],...          
                    baseIndex(j));                                      % peak values to fit line
                fit = fitx\fity;                                        % fit line between the two peaks
                top = fit(1) + fit(2)*peakt(2);                         % extrapolate from negative peak to fitted line
                bottom = baseS(round(peakt(2)),baseIndex(j));           % value at negative peak
                spike(j) = abs(top-bottom);                             % pop spike amplitude
            end

            % plot post-induction spikes
            subplot(1,2,2)
            plot(t1:t2,baseS(t1:t2,baseIndex(21:end))...
                -ones(length(t1:t2),1)*baseS(t1,baseIndex(21:end)));
            [peakt,peakv] = ginput(3);
            for j = 21:length(baseIndex);                               % repeat for pop spikes after induction
                fitx = [ones(2,1) round(peakt([1 3]))];
                fity = baseS([round(peakt(1)) round(peakt(3))],baseIndex(j));
                fit = fitx\fity;
                top = fit(1) + fit(2)*peakt(2);
                bottom = baseS(round(peakt(2)),baseIndex(j));
                spike(j) = abs(top-bottom);
            end
            spikeN = spike/mean(spike(1:20));
            
            %% take slopes of baseline traces
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
            tbsOn = 2;                                                      % tbs start time in s
            tbsOff = 5;                                                     % tbs end time in s
            tBurst = .04;                                                   % duration of each burst
            nBurst = 15;                                                    % number of burst
            nPulse = 4;                                                     % number of pulses in a burst
            rBurst = 100;                                                   % burst frequency
            
            % indices to remove bipolar artifact
            stimI = tbsOn*fs+1:round(.01*fs):tbsOn*fs+ 3*round(.01*fs);
            
            % soma traces during induction
            if soma ~= 0;                                                   % if there is a somatic recording
                indS =  data(datastart(soma,indBlock(1)):...
                    dataend(soma,indBlock(1)))';                               % extract data from induction block
                indS1 = indS(tbsOn*fs+1:tbsOff*fs);                         % take data only during bursts
                indS2 = reshape(indS1,[],nBurst);                           % organize into 3D matrix, 1: time, 2: pulse number, 3: burst number
                indS3 = indS2(1:round(tBurst*fs),:);                        % remove time in between bursts so that length of each pulse is the same
                indS4 = reshape(indS3,[],nPulse,nBurst);                    % reshape back to 3D matrix
                
            end
            
            % repeat for dendritic traces during induction 
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
            save(strcat(fpathP,slicesR{i}));
        end
    end
end