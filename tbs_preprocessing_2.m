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
for a = 1:size(directR,1);
    slicesR{a} = directR(a).name; % raw file names
end
slicesP = cell(size(directP,1),1); %create empty cell array
for a = 1:size(directP,1);
    slicesP{a} = directP(a).name; % processed file names
end

%% process new slices
if length(slicesR)~= length(slicesP);                                       % are there new slices?
    for i = 1:length(slicesR);                                              % loop over all slices
        if sum(strcmp(slicesP,slicesR{i}))==0;                              % check if current slice is new
            
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
            for j = 1:size(comtext,1);                                      % loop through all comments
                
                if isempty(strfind(comtext(j,:),'induction_')) == 0;         % find plasticity induction block
                   indBlock(1) = com(com(:,5)==j,2);
                end
                if isempty(strfind(comtext(j,:),'induction2_')) == 0;         % find plasticity induction block
                   indBlock(2) = com(com(:,5)==j,2);
                end
                if isempty(strfind(comtext(j,:),'induction3_')) == 0;         % find plasticity induction block
                   indBlock(3) = com(com(:,5)==j,2);
                end
                if isempty(strfind(comtext(j,:),'induction4_')) == 0;         % find plasticity induction block
                   indBlock(4) = com(com(:,5)==j,2);
                end
                
                if isempty(strfind(comtext(j,:),'Soma')) == 0;          % find channel with somatic recording
                    soma = com(com(:,5)==j,1);
                end
                
                if isempty(strfind(comtext(j,:),'Apical')) == 0;        % channel with apical recording
                    apical = com(com(:,5)==j,1);
                end
                
                if isempty(strfind(comtext(j,:),'Basal')) == 0;         % channel with basal recording
                    basal = com(com(:,5)==j,1);
                end
                
                if isempty(strfind(comtext(j,:),'Perforant')) == 0;         % channel with perforant recording
                    perforant = com(com(:,5)==j,1);
                end
                
                if isempty(strfind(comtext(j,:),'age')) == 0;           % age of animal
                    ageI = strfind(comtext(j,:),'age');
                    ageL = length('age_');
                    age = str2double(comtext(j,ageI+ageL:ageI+ageL+1));
                end
                
                if isempty(strfind(comtext(j,:),'hemi')) == 0;          % brain hemisphere
                    hemiI = strfind(comtext(j,:),'hemi');
                    hemiL = length('hemi_');
                    hemi = comtext(j,hemiI+hemiL);
                end
                if isempty(strfind(comtext(j,:),'height')) == 0;        % position along dorsal ventral axis of hippocampus (x0.1 mm from IA from paxinos atlas)
                    heightI = strfind(comtext(j,:),'height');
                    heightL = length('height_');
                    height = str2double(comtext(j,heightI+heightL:heightI+heightL+1));
                end
                
                if isempty(strfind(comtext(j,:),'current')) == 0;       % current injected for DCS (in uA)
                    currentI = strfind(comtext(j,:),'current');
                    currentL = length('current_')
                    current = str2double(comtext(j,currentI+currentL:end));
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
            if date>=20161013;                                              % changed the timing of probe pulses, slices that have a rig number have the new pulse time
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
            autospikes = 0;                                                 % takes spikes automatically? (1=yes, 0=no)
            
            % manually find spikes
            if autospikes == 0;
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
            end
            
            % automatically find spikes
            if autospikes == 1;
                % filter raw traces
                L1 = 5; % length of MA filter
                co1 = ones(L1,1)/L1; % MA coefficients
                co2 = 1; % Ar coefficients
                baseSsm = filtfilt(co1,co2,baseS); % zero phase low pass filter

                % derivative of smoothed trace
                dbaseS = [diff(baseS,1);zeros(1,size(baseS,2))]; % derivative of smoothed original signal
                dbaseS = filtfilt(co1,co2,dbaseS);% filter derivative

                % find where derivative crosses zero, low pass filter
                L2 = 10;                                                    % length of MA filter
                co3 = ones(L2,1)/L2;                                        % MA coefficients
                co4 = 1; % Ar coefficients
                dbaseS_sign = sign(filtfilt(co3,co4,sign(dbaseS)));         % sign of the derivative (smoothed)
                dbaseS_sign(dbaseS==0)=1;                                   % set zeros equal to one
                ddbaseS_sign = diff(dbaseS_sign,1,1);
                dbaseS_zero = [logical(ddbaseS_sign~=0);...
                    false(1,size(dbaseS_sign,2))];                          % find zero crossings (sign change)
                baseSt = repmat(((1:size(baseSsm,1))/fs)',1,size(baseSsm,2));% matrix of time vectors
                
                % find peak values of pop spike trace
                t1 = pulset + .001*fs;
                t2 = pulset + .01*fs;
                figure(1)
                plot(t1:t2,baseS(t1:t2,...
                    baseIndex)-mean(baseS(t1:t2,baseIndex),1));
                [tspike,yspike] = ginput(3);
                delay1 = .0025*fs;
                delay2 = .02*fs;
                t1 = pulset + delay1;
                t2 = pulset + delay2;
                baseST = baseSt(t1:t2,baseIndex);
                baseSSM = baseSsm(t1:t2,baseIndex);
                dbaseS_ZERO = dbaseS_zero(t1:t2,baseIndex);
                peakt = zeros(3,length(baseIndex));
                peakv = zeros(3,length(baseIndex));
                spike = zeros(length(baseIndex),1);
                for j = 1:length(baseIndex);
                    peakt1 = baseST(dbaseS_ZERO(:,j),j);
                    peakv1 = baseSSM(dbaseS_ZERO(:,j),j);
                    if length(peakt1)>=3 & mean(diff(peakt1(1:3))) > .001;
%                         peakt2 = [peakt1(1);peakt1((peakt1-peakt1(1))>.0012)];
%                         peakv1 = [peakv1(1);peakv1((peakt1-peakt1(1))>.0012)];
                        peakt(:,j) = peakt1(1:3);
                        peakv(:,j) = peakv1(1:3);        
                    else
                        plot(baseSSM(:,j));
                        title('click peaks in order max>min>max')
                        [peakt(:,j),peakv(:,j)] = ginput(3);
                    end
                    fitx = [ones(2,1) peakt([1 3],j)];
                    fity = peakv([1 3],j);
                    fit = fitx\fity;
                    top = fit(1) + fit(2)*peakt(2,j);
                    bottom = peakv(2,j);
                    spike(j) = abs(top-bottom);
                end
                spikeN = spike/mean(spike(1:20));          
            end
            
            %% take slopes of baseline traces
            autoslopes = 0;                                                 % takes slopes automatically? (1=yes, 0=no)
            
            % manually take slopes   
            if autoslopes == 0                                             % manually take slopes by selecting time window with cursors
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
                % remove slopes that are zero (noise/amplifier saturation)
%                 slopes_before = slopesD(1:indBlock-1);
%                 slopes_after = slopesD(indBlock+1:end);
%                 slopes_base = slopesD(baseIndex);
%                 [removeA,removeAI] = find(abs(slopes_after) <...
%                         0.2*abs(mean(slopes_base)));
%                 [removeB,removeBI] = find(abs(slopes_before) <...
%                         0.2*abs(mean(slopes_base)));
%                 if length(slopes_after)-length(removeA)>=60
%                     slopes_after(removeAI) = []; 
%                 end
%                 if length(slopes_before)-length(removeB)>=20
%                     slopes_before(removeBI) = []; 
%                     indBlock_shift(1) = indBlock(1)-length(removeBI);
%                 end
%                 slopesD = [slopes_before;slopes_after];
%                 baseIndex_shift = [indBlock_shift(1)-20:indBlock_shift(1)-1,...
%                     indBlock_shift(1)+1:indBlock_shift(1)+60];    % index of baseline responses (pre/post induction)
            end            
            % automatically find slopes
            if autoslopes == 1;                                               % take slopes automatically
                slope1 = pulset + .0018*fs;
                slope2 = pulset + .0024*fs;
                ds = .0001*fs;
                slopeWin = .0006*fs; 
                slopesA = zeros(size(datastart,2),round((slope2-slope1)./ds));
                slopesB = zeros(size(datastart,2),round((slope2-slope1)./ds));
                if apical ~= 0;
                    for j = 1:size(datastart,2);
                        for k = 1:round((slope2-slope1)./ds)+1;
                            slopesA(j,k) = mean(diff(baseA(slope1+(k-1)*ds:...
                                slope1+(k-1)*ds+slopeWin,j)))./(1/fs);
                        end      
                    end
                    temp = mean(slopesA(indBlock+51:indBlock+60,:),1);%/mean(slopesA(indBlock-20:indBlock-1,:),1);
                    [a,b] = max(abs(temp));
                    slopesA = slopesA(:,b);
                end
                              
                
                
                if basal ~= 0;
                    for j = 1:size(datastart,2);
                        for k = 1:round((slope2-slope1)./ds);
                            slopesB(j,k) = mean(diff(baseB(slope1+(k-1)*ds:...
                                slope1+(k-1)*ds+slopeWin,j)))./(1/fs);
                        end
                    end
                    temp = mean(slopesB(indBlock+51:indBlock+60,:),1);
                    [a,b] = max(abs(temp));
                    slopesB = slopesB(:,b);
                end
                
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
            
            %% take slopes during induction

            %% high pass filter somatic recordings from induction block
            
            % design high pass IIR filter
%             Fstop = 200;
%             Fpass = 300;
%             Astop = 65;
%             Apass = 0.5;
%             dfilt = designfilt('highpassiir','StopbandFrequency',Fstop ,...
%               'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
%               'PassbandRipple',Apass,'SampleRate',fs,'DesignMethod','butter');
            
            % design high pass FIR filter
%             Fstop = 200;
%             Fpass = 300;
%             Astop = 65;
%             Apass = 0.5;
%             dfilt = designfilt('highpassfir','StopbandFrequency',Fstop ,...
%               'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
%               'PassbandRipple',Apass,'SampleRate',fs,'DesignMethod','kaiserwin');
% 
%             % filter somatic recording 
%             % remove bipolar artifact
%             bipolarW = .000*fs;
%             jumps = (0:nPulse-1)'*(1/rBurst)*fs;
%             bipolarI = (1:bipolarW)+jumps;
%             
%             S2I = 1:size(indS2,1);
%             indSfilt1 = indS2(setdiff(S2I,bipolarI),:);                   % traces with bipolar stimulus artifact removed
%             indSfilt1 = indSfilt1(:);
%             indSfilt = [indS(1:tbsOn*fs);indSfilt1;indS(tbsOff*fs:end)];
%             
%             % apply filter to somatic recordings
%             indShigh = filtfilt(dfilt,indSfilt);
%             indS1high = indShigh(tbsOn*fs+1:(tbsOn*fs+length(indSfilt1)));                          % take data only during bursts
%             indS2high = reshape(indS1high,[],nBurst);                             % organize into 3D matrix, 1: time, 2: pulse number, 3: burst number
%             indS3high = indS2high(1:round(tBurst*fs-length(bipolarI(:))),:);                         % remove time in between bursts so that length of each pulse is the same
%             indS4high = reshape(indS3high,[],nPulse,nBurst);                    % reshape back to 3D matrix
%             figure
%             plot(indS4high(:,:,1)-indS4high(1,:,1))
%             
%             % dendritic recordings
%             
%             % apical 
%             if apical ~= 0;
%                 A2I = 1:size(indA2,1);
%                 indAfilt1 = indA2(setdiff(A2I,bipolarI),:);                   % traces with bipolar stimulus artifact removed
%                 indAfilt1 = indAfilt1(:);
%                 indAfilt = [indA(1:tbsOn*fs);indAfilt1;indA(tbsOff*fs:end)];
% 
%                 % apply filter to somatic recordings
%                 indAhigh = filtfilt(dfilt,indAfilt);
%                 indA1high = indAhigh(tbsOn*fs+1:(tbsOn*fs+length(indAfilt1)));                          % take data only during bursts
%                 indA2high = reshape(indA1high,[],nBurst);                             % organize into 3D matrix, 1: time, 2: pulse number, 3: burst number
%                 indA3high = indA2high(1:round(tBurst*fs-length(bipolarI(:))),:);                         % remove time in between bursts so that length of each pulse is the same
%                 indA4high = reshape(indA3high,[],nPulse,nBurst);                    % reshape back to 3D matrix
%                 figure
%                 plot(indA4high(:,:,1)-indA4high(1,:,1))
%             end
%             
%             % basal
%             if basal ~= 0;
%                 B2I = 1:size(indB2,1);
%                 indBfilt1 = indB2(setdiff(B2I,bipolarI),:);                   % traces with bipolar stimulus artifact removed
%                 indBfilt1 = indBfilt1(:);
%                 indBfilt = [indB(1:tbsOn*fs);indBfilt1;indB(tbsOff*fs:end)];
% 
%                 % apply filter to somatic recordings
%                 indBhigh = filtfilt(dfilt,indBfilt);
%                 indB1high = indBhigh(tbsOn*fs+1:(tbsOn*fs+length(indBfilt1)));                          % take data only during bursts
%                 indB2high = reshape(indB1high,[],nBurst);                             % organize into 3D matrix, 1: time, 2: pulse number, 3: burst number
%                 indB3high = indB2high(1:round(tBurst*fs-length(bipolarI(:))),:);                         % remove time in between bursts so that length of each pulse is the same
%                 indB4high = reshape(indB3high,[],nPulse,nBurst);                    % reshape back to 3D matrix
%                 figure
%                 plot(indB4high(:,:,1)-indB4high(1,:,1))
%             end
%             
            %% linear fit between soma and dendrite recording
            % covariance
%             covyx = (indS3high-mean(indS3high,1))'*((indA3high-mean(indA3high,1)));
%             covxx = (indA3high-mean(indA3high,1))'*((indA3high-mean(indA3high,1)))';
%             Ahat = covyx/covxx; 
%             nsub = (indS3high-mean(indS3high,1))' - Ahat*(indA3high-mean(indA3high,2));
%             figure;plot(nsub-nsub(1,:))
            
            
            %% save workspace variables
            save(strcat(fpathP,slicesR{i}));
        end
    end
end