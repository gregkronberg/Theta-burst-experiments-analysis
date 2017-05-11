%% Design filters for synaptic plasticity experiments

fpath_filt = 'D:\Google Drive\Work\Research Projects\Theta LTP\Filters\';

filt_labels = {'iir_highpass_1Hz_fs10k','iir_bandpass_200Hz_1000Hz_fs10k',...
    'iir_bandpass_200Hz_600Hz_fs10k','iir_bandpass_100Hz_1000Hz_fs10k','iir_lowpass_1000Hz_fs10k'};

[temp,filt_saved] = folder_diff(fpath_filt,fpath_filt,'.mat');

% loop over filters
for a = 1:length(filt_labels)
    % check if filter is already saved
    if sum(strcmp(filt_saved,filt_labels{a}))==0
        % create and save filters
        
        % iir_highpass_1Hz_fs10k
        %==================================================================
        if strcmp(filt_labels{a},'iir_highpass_1Hz_fs10k')
            fs =10000;
            Fstop = 0.1;
            Fpass = 1;
            Astop = 30;
            Apass = 0.5;
            dfilt = designfilt('highpassiir','StopbandFrequency',Fstop ,...
              'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
              'PassbandRipple',Apass,'SampleRate',fs);
            save(strcat(fpath_filt,filt_labels{a},'.mat'),'dfilt')
            
        % iir_bandpass_200Hz_1000Hz_fs10k   
        %==================================================================
        elseif strcmp(filt_labels{a},'iir_bandpass_200Hz_1000Hz_fs10k')
            fs =10000;
            dfilt = designfilt('bandpassiir','FilterOrder',20,...
                'HalfPowerFrequency1',200,'HalfPowerFrequency2',1000, ...
                'SampleRate',fs);
            save(strcat(fpath_filt,filt_labels{a},'.mat'),'dfilt')
            
        % iir_bandpass_200Hz_600Hz_fs10k 
        %==================================================================
        elseif strcmp(filt_labels{a},'iir_bandpass_200Hz_600Hz_fs10k')
            fs = 10000;
            dfilt = designfilt('bandpassiir','FilterOrder',20,...
                'HalfPowerFrequency1',200,'HalfPowerFrequency2',600, ...
                'SampleRate',fs);
            save(strcat(fpath_filt,filt_labels{a},'.mat'),'dfilt')
        
        % iir_bandpass_100Hz_1000Hz_fs10k
        %==================================================================
        elseif strcmp(filt_labels{a},'iir_bandpass_100Hz_1000Hz_fs10k')
            fs = 10000;
            dfilt = designfilt('bandpassiir','FilterOrder',20,...
                'HalfPowerFrequency1',200,'HalfPowerFrequency2',600, ...
                'SampleRate',fs);
            save(strcat(fpath_filt,filt_labels{a},'.mat'),'dfilt')
        
        % iir_lowpass_1000Hz_fs10k
        %==================================================================
        elseif strcmp(filt_labels{a},'iir_lowpass_1000Hz_fs10k')
            fs =10000;
            designfilt('lowpassiir','FilterOrder',20, ...
                'PassbandFrequency',1000,'PassbandRipple',0.2, ...
                'SampleRate',fs);
            save(strcat(fpath_filt,filt_labels{a},'.mat'),'dfilt')
        
        end 
    end
end