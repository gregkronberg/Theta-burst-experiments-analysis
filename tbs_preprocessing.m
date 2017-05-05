%% preprocessing
%==========================================================================
% Time series processing
    % Raw time series data from LabChart are loaded and organized into
    % matrices.  Data during probing blocks (i.e. testing baseline synaptic
    % efficacy are organized as (time x blocks).  Data during induction has 
    % multiple formats including (time, 1D),(time x bursts),(time x pulses x bursts).  
    % Data are also filtered for later processing.
    
    % Processed time series data for each slice are saved as individual
    % files in the folder "Processed Matlab Data"

% Info about experimental conditions
    % Experimental parameters are extracted from comments for each slice and
    % stored in one structure called "slices", which is saved in the "Matlab
    % Variables" folder.  slices is used as a reference to extract info about
    % each slice in downstream processing steps.  Slices is organized as
    % slices{conditions}(slice number).parameter. "conditions" is a cell that
    % keeps track of the various experimental conditions and is stored in the
    % same file as slices to be used in downstream processing
%==========================================================================

clear all; close all; clc

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

%% define experimental conditions
%==========================================================================
induction = {'_TBS_'};                      % plasticity induction protocol
stim = {'control';'cathodal';'anodal'};     % DCS stimulation conditions
stim_color = {[0 0 0],[0 0 1],[1 0 0]};     % plot colors for each DCS condition
intensity = [0,5,20];                       % stimulation intensity in V/m
position = {'apical';'basal';'perforant'};  % recording location in slice
position_mark = {'.','x','*'};              % plotting symbol for each position
drug = {'_none';'_mk801'};                  % drugs used
% store conditions in single cell
conditions = {induction,stim,intensity,position,drug};

%% list slices for each condition
%==========================================================================
% preallocate
% list of new slices
slices_new = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
% complete list of raw slices
slices_raw = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
% complete list of processed slices
slices_processed = cell(length(induction),length(stim),length(intensity),length(position),length(drug));


for a = 1:length(induction)
    for b = 1:length(stim)
        for c = 1:length(intensity)
            for d = 1:length(position)
                for e = 1:length(drug)
%===================================== loop over experimental conditions                  

% list slices for each condition
[slices_new{a,b,c,d,e},slices_raw{a,b,c,d,e},slices_processed{a,b,c,d,e}] = ... 
    folder_diff(fpath_raw,fpath_processed,strcat(...
    induction{a},stim{b},strcat('_',num2str(intensity(c)),'Vm_'),...
    position{d},drug{e},'*.mat'));

%===================================== end loop over conditions
                end
            end
        end
    end
end

%% load filters
%==========================================================================
bandpass_200_600 = load(strcat(fpath_filters,'iir_bandpass_200Hz_600Hz_fs10k.mat'));
highpass_1 = load(strcat(fpath_filters,'iir_highpass_1Hz_fs10k.mat'));


%% process new slices
%==========================================================================
for a = 1:length(induction)
    for b = 1:length(stim)
        for c = 1:length(intensity)
            for d = 1:length(position)
                for e = 1:length(drug)
%===================================== loop over experimental conditions
% check for slices in this condition
if isempty(slices_new{a,b,c,d,e})==0
    % loop over new slices in current condition
    for f = 1:length(slices_new{a,b,c,d,e})
%====================================== loop over individual slices
       
        % number of slices already stored in this condition
        num = length(slices{a,b,c,d,e});
        
        % load new slices
        load(strcat(fpath_raw,slices_new{a,b,c,d,e}{f}));

        % sampling rate
        fs = mode(samplerate(:));           

        % length of each block in samples
        l = dataend-datastart+1;            

        %% extract info from comments
        %==================================================================
        % during of sampling before/after induction
        if induction{a} == '_TBS_'
            tpre = 20;% recording time before induction (minutes)
            tpost = 60;% recording time after induction (minutes)
        end
        % find induction blocks and type of induction
        numInd = 4;                                                     % number of possible induction periods
        indBlock = zeros(numInd,1);                                     % recording block number for plasticity induction
        indType = cell(numInd,1);                                       % type of plasticity induction, e.g. TBS 
        indBlock(1) = comment_find('induction_','block',com,comtext);   % first induction
        indType{1} = comment_find('induction','value',com,comtext);
        indBlock(2) = comment_find('induction2','block',com,comtext);   % second induction
        indType{2} = comment_find('induction2','value',com,comtext);
        indBlock(3) = comment_find('induction3','block',com,comtext);   % third induction
        indType{3} = comment_find('induction3','value',com,comtext);
        indBlock(4) = comment_find('induction4','block',com,comtext);   % fourth induction
        indType{4} = comment_find('induction4','value',com,comtext);
        baseIndex = [indBlock(1)-tpre:indBlock(1)-1,indBlock(1)+1:indBlock(1)+tpost];    % index of baseline responses (pre/post induction)

        % channel number for each recording location
        soma = comment_find('Soma','channel',com,comtext);               
        apical = comment_find('Apical','channel',com,comtext);
        basal = comment_find('Basal','channel',com,comtext);
        perforant = comment_find('Perforant','channel',com,comtext);

        % physiological parameters
        age = comment_find('age','value',com,comtext);
        hemi = comment_find('hemi','value',com,comtext);
        height = comment_find('height','value',com,comtext);

        % stimulation parameters
        current = comment_find('current','value',com,comtext);
        isolator = comment_find('isolator','value',com,comtext);
        path2 = comment_find('2path','value',com,comtext);

        % date and slice number
        date = str2double(slices_new{a,b,c,d,e}{f}(1:8));                             % date of recording
        slice_num = str2double(slices_new{a,b,c,d,e}{f}(10));                        % slice number used in chamber

        % timing of bipolar pulse
        if date>=20161013                                              % changed the timing of probe pulses, slices that have a rig number have the new pulse time
            pulset = .5*fs;                                             % start time of bipolar pulse (samples)
        else
            pulset = 0*fs ;
        end
        pulset_path2 = 30.5*fs;                                         % start time of bipolar pulse 

        %% store info from comments
        %==================================================================
        % add new info to global slice info structure
        slices{a,b,c,d,e}(num+1).name = slices_new{a,b,c,d,e}{f};
        slices{a,b,c,d,e}(num+1).fs = fs;
        slices{a,b,c,d,e}(num+1).indBlock = indBlock;
        slices{a,b,c,d,e}(num+1).indType = indType;
        slices{a,b,c,d,e}(num+1).baseIndex = baseIndex;
        slices{a,b,c,d,e}(num+1).soma = soma;
        slices{a,b,c,d,e}(num+1).apical = apical;
        slices{a,b,c,d,e}(num+1).basal = basal;
        slices{a,b,c,d,e}(num+1).perforant = perforant;
        slices{a,b,c,d,e}(num+1).age = age;
        slices{a,b,c,d,e}(num+1).hemi = hemi;
        slices{a,b,c,d,e}(num+1).height = height;
        slices{a,b,c,d,e}(num+1).current = current;
        slices{a,b,c,d,e}(num+1).isolator = isolator;
        slices{a,b,c,d,e}(num+1).path2 = path2;
        slices{a,b,c,d,e}(num+1).date = date;
        slices{a,b,c,d,e}(num+1).slice_num = slice_num;
        slices{a,b,c,d,e}(num+1).pulset = pulset;
        slices{a,b,c,d,e}(num+1).pulset_path2 = pulset_path2;

        %% organize baseline time series traces into matrices
        %======================================================================
        % baseS and baseD are matrices of baseline response time series (time x block) 
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

        %% store induction block as matrices
        %==================================================================
        % theta burst parameters
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

        % dendritic traces during induction 
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
        
        
        %% apply filters
        %==================================================================
        % iir_bandpass_200Hz_600Hz_fs10k
        indD_filt_band_200_600 = filtfilt(bandpass_200_600.dfilt,indD);
        indS_filt_band_200_600 = filtfilt(bandpass_200_600.dfilt,indS);
        baseS_filt_band_200_600 = filtfilt(bandpass_200_600.dfilt,baseS);
        baseD_filt_band_200_600 = filtfilt(bandpass_200_600.dfilt,baseD);

        % iir_highpass_1Hz_fs10k
        indD_filt_high_1 = filtfilt(highpass_1.dfilt,indD);
        indS_filt_high_1 = filtfilt(highpass_1.dfilt,indS);
        baseS_filt_high_1 = filtfilt(highpass_1.dfilt,baseS);
        baseD_filt_high_1 = filtfilt(highpass_1.dfilt,baseD);
    
        %% save to processed data folder
        %==================================================================
        % save all variables except slices
        save(strcat(fpath_processed,slices_new{a,b,c,d,e}{f}), '-regexp', '^(?!(slices)$).') 


        %===================================== end loop over individual slices
    end
end
%===================================== end loop over experimental conditions
                end
            end
        end
    end
end

%% save updated global slices variable
%==========================================================================
save(strcat(fpath_variables,'slices.mat'),'slices','conditions');

    
    
    
                       