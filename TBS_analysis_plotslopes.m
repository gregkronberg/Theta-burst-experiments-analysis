%% analysis: fEPSP slopes
%==========================================================================
% notes
%==========================================================================

clear all
close all
clc

%==========================================================================
%% file paths
%==========================================================================
% desktop
fpath_raw = 'D:\Google Drive\Work\Research Projects\Theta LTP\Raw Matlab Data\'; % raw
fpath_processed = 'D:\Google Drive\Work\Research Projects\Theta LTP\Processed Matlab Data\'; % processed
fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
fpath_analysis = 'D:\Google Drive\Work\Research Projects\Theta LTP\Analysis\';% analysis
fpath_filters = 'D:\Google Drive\Work\Research Projects\Theta LTP\Filters\'; % filters

%==========================================================================
%% load global structures
%==========================================================================
% slices
load(strcat(fpath_variables,'slices'));
% slopes
load(strcat(fpath_variables,'slopes'));

%% exclusion criteria
%==========================================================================
date_cut = [20170115 20170301 20170401];

%==========================================================================
%% preallocate
%==========================================================================
tpre = 20;
tpost = 60;
t  = 1:tpre+tpost;
slopes_norm = cell(length(conditions{1}),length(conditions{2}),length(conditions{3}),length(conditions{4}),length(conditions{5}));
stim_color = {[0 0 0],[0 0 1],[1 0 0]};
%% store slopes for each condition in matrix
%==========================================================================
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    
%===================================== loop over experimental conditions
if isempty(slopes{a,b,c,d,e})==0
    % apply exclusion criteria
    include = [slices{1,1,1,1,1}(:).date]'>date_cut(d);
    slices_temp = slices{a,b,c,d,e}(include);
    slopes_temp = slopes{a,b,c,d,e}(include);
    % preallocate
    slopes_norm{a,b,c,d,e} = zeros(tpre+tpost,length(slopes_temp{a,b,c,d,e}));
    for f = 1:length(slopes_temp{a,b,c,d,e})
%====================================== loop over individual slices
        
base_index = slopes_temp{a,b,c,d,e}(f).baseIndex;

% store normalized slopes as a matrix
slopes_base_mean = mean(slopes_temp{a,b,c,d,e}(f).slopes(base_index(1:tpre)));
slopes_norm{a,b,c,d,e}(:,f) = slopes_temp{a,b,c,d,e}(f).slopes(base_index)'/slopes_base_mean; % (blocks x slices)
        
        %===================================== end loop over individual slices
    end
end
%===================================== end loop over experimental conditions
                end
            end
        end
    end
end

%==========================================================================
%% Stats and plots
%==========================================================================
for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slopes_temp{a,b,c,d,e})==0
% stats
%===================================
% mean and standard error
slopes_mean = mean(slopes_norm{a,b,c,d,e},2);
slopes_sem = std(slopes_norm{a,b,c,d,e},0,2)/sqrt(size(slopes_norm{a,b,c,d,e},2));
slopes_end{a,b,c,d,e} = mean(slopes_norm{a,b,c,d,e}(end-10:end,:),1);
slopes_end_mean{a,b,c,d,e} = mean(slopes_end{a,b,c,d,e});
slopes_end_sem{a,b,c,d,e} = std(slopes_end{a,b,c,d,e},0,2)/sqrt(length(slopes_end{a,b,c,d,e}));
% pairwise comparisons to no DCS (1st parameter of condition b)
[h,p{a,b,c,d,e}] = ttest2(slopes_end{a,1,c,d,e},slopes_end{a,b,c,d,e});

% plot slopes over time
%====================================
figure;hold on
errorbar(t,slopes_mean{a,b,c,d,e},slopes_sem{a,b,c,d,e},...
    '.','Color',stim_color{b},'MarkerSize',30);
errorbar(t,slopesMean{a,1,1,d,e},slopesStd{a,1,1,d,e}/sqrt(size(slopesN{a,1,1,d,e},2)),...
    '.','Color',stim_color{1},'MarkerSize',30);


% format figure
xSize = 15; ySize = 10;
set(gcf,'PaperUnits','centimeters')
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
%                             set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0.5 0.5 xSize*50 ySize*50])
%set(gca,'TickDir','out','TickLength',[0.03 0.03]); % make the tick sizes bigger and outside
set(gcf, 'Color', 'w'); 
xlim([t(1) t(end)])
ylim([.6 1.8])
set(gca,'TickDir','out')
set(gca,'TickLength',[.03 .03],'FontSize',10,'FontWeight','bold','LineWidth',2)
xlabel('Time (min)','FontSize',30,'FontWeight','bold')
ylabel('Normalize fEPSP slope','FontSize',30,'FontWeight','bold')
title(strcat('TBS with ',conditions{2}{b},', ',num2str(conditions{3}(c)),'V/m, ',conditions{4}{d}));
                    end
                end
            end
        end
    end
end



                    
                    
                    
                    
        %% arrange file names by condition
slices = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
for a = 1:length(induction);
    for b = 1:length(stim);
        for c = 1:length(intensity);
            for d = 1:length(position);
                for e = 1:length(drug);
                        slices{a,b,c,d,e} = dir(strcat(fpath,'*',...
                           induction{a},stim{b},strcat('_',num2str(intensity(c)),'Vm_'),...
                           position{d},drug{e},'*.mat'));
                end
            end
        end
    end
end

%% store slopes and pop spikes by condition
slopes = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesN = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
spikes = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
spikesN = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesInd = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesIndN1 = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesIndN2= cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesMean = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesStd = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesEnd = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
spikesMean = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
spikesStd = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
spikesEnd = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
heights = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
dates = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
hemis = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
rig = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
ages = cell(length(induction),length(stim),length(intensity),length(position),length(drug));
currents= cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopes_raw= cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopes_extra= cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopesN_extra= cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopes_extra_mean= cell(length(induction),length(stim),length(intensity),length(position),length(drug));
slopes_extra_std= cell(length(induction),length(stim),length(intensity),length(position),length(drug));
for a = 1:length(induction);
    for b = 1:length(stim);
        for c = 1:length(intensity);
            for d = 1:length(position);
                for e = 1:length(drug);
                    if isempty(slices{a,b,c,d,e}) == 0
                        
                        slopes{a,b,c,d,e} = zeros(tbase+tpost,length(slices{a,b,c,d,e}));
                        slopesN{a,b,c,d,e} = zeros(tbase+tpost,length(slices{a,b,c,d,e}));
                        spikes{a,b,c,d,e} = zeros(tbase+tpost,length(slices{a,b,c,d,e}));
                        spikesN{a,b,c,d,e} = zeros(tbase+tpost,length(slices{a,b,c,d,e}));
                        rig{a,b,c,d,e} = zeros(length(slices{a,b,c,d,e}),1);
                        num_ind{a,b,c,d,e} = zeros(length(slices{a,b,c,d,e}),1);
                        slopesN_extra{a,b,c,d,e} = cell(4,1);
                        slopes_extra{a,b,c,d,e} = cell(4,1);
                        slopes_extra_mean{a,b,c,d,e} = cell(4,1);
                        slopes_extra_std{a,b,c,d,e} = cell(4,1);
                        for f = 1:length(slices{a,b,c,d,e});
                            load(strcat(fpath,slices{a,b,c,d,e}(f).name),...
                                'slopesD','indBlock',...
                                'spike','spikeN','height','hemi','age',...
                                'current','date');
                            
                            % which rig?
                            if isempty(strfind(slices{a,b,c,d,e}(f).name,'rig'))==1
                                rig{a,b,c,d,e}(f) = 1;
                            elseif isempty(strfind(slices{a,b,c,d,e}(f).name,'rig'))==0
                                rig_loc = strfind(slices{a,b,c,d,e}(f).name,'rig');
                                rig{a,b,c,d,e}(f) = str2double(slices{a,b,c,d,e}(f).name(rig_loc+3));
                            end
                            
                            % how many induction periods
%                             num_ind{a,b,c,d,e}(f) = max(find(indBlock~=0));
                            if length(slopesD)>=indBlock(1)+tpost
                                slopes{a,b,c,d,e}(:,f) = slopesD([indBlock(1)-tbase:indBlock(1)-1,indBlock(1)+1:indBlock(1)+tpost]);
                            end
                            for g  = 2:length(indBlock)
                                if indBlock(g)~=0
                                    if indBlock(g)+tpost2<length(slopesD)
                                        num_ind{a,b,c,d,e}(f)=g;
                                        slopes_extra{a,b,c,d,e}{g}(:,f)  = slopesD(indBlock(g)+1:indBlock(g)+tpost2);
                                        slopesN_extra{a,b,c,d,e}{g}(:,f) =  slopes_extra{a,b,c,d,e}{g}(:,f)/mean(slopes{a,b,c,d,e}(1:tbase,f));
                                    end
                                end
                                
                            end
                                    
                            slopes_raw{a,b,c,d}{f} = slopesD;
                            spikes{a,b,c,d,e}(:,f) = spike;%popSpike([indBlock-tbase:indBlock-1,indBlock+1:indBlock+tpost]);
                            slopesN{a,b,c,d,e}(:,f) = slopes{a,b,c,d,e}(:,f)/mean(slopes{a,b,c,d,e}(1:tbase,f));
                            spikesN{a,b,c,d,e}(:,f) = spikeN;%spikes{a,b,c,d,e}(:,f)/mean(spikes{a,b,c,d,e}(1:tbase,f));
%                             slopesIndN1{a,b,c,d,e}(:,f) = slopesInd{a,b,c,d,e}(:,f)/mean(slopes{a,b,c,d,e}(1:tbase,f));
%                             slopesIndN2{a,b,c,d,e}(:,f) = slopesInd{a,b,c,d,e}(:,f)/slopesInd{a,b,c,d,e}(1,f);
                            heights{a,b,c,d,e}(f) = height;
                            hemis{a,b,c,d,e}(f)  =hemi;
                            dates{a,b,c,d,e}(f) = date;
                            ages{a,b,c,d,e}(f) = age;
                            currents{a,b,c,d,e}(f) = current;
                        end
                        spikesMean{a,b,c,d,e} = mean(spikesN{a,b,c,d,e},2);
                        spikesStd{a,b,c,d,e} = std(spikesN{a,b,c,d,e},0,2);
                        spikesEnd{a,b,c,d,e} = mean(spikesN{a,b,c,d,e}(end-9:end,:),1);
                        
                    end
                end
            end
        end
    end
end
save(strcat(fpath_variables,'fepsp_slopes'),'slopes','slopes_raw','slopesN',...
    'spikes','spikesN','slices')

height_range = 22:28;
cutoff{apical} = 20170115;
cutoff{basal} = 20170301;
cutoff{perforant} = 20170401;
%% stats
for a = 1:length(induction);
    for b = 1:length(stim);
        for c = 1:length(intensity);
            for d = 1:length(position);
                for e = 1:length(drug);
                    if isempty(slices{a,control,1,d,e}) == 0;
                        if isempty(slices{a,b,c,d,e}) == 0;
                            date_cut{a,b,c,d,e} = dates{a,b,c,d,e}>cutoff{d};
                            slopesMean{a,b,c,d,e} = mean(slopesN{a,b,c,d,e}(:,date_cut{a,b,c,d,e}),2);
                            slopesStd{a,b,c,d,e} = std(slopesN{a,b,c,d,e}(:,date_cut{a,b,c,d,e}),0,2);
                            slopesEnd{a,b,c,d,e} = mean(slopesN{a,b,c,d,e}(end-9:end,(date_cut{a,b,c,d,e})),1);
                            slopesEndMean{a,b,c,d,e} = mean(slopesEnd{a,b,c,d,e},2);
                            slopesEndStd{a,b,c,d,e} = std(slopesEnd{a,b,c,d,e},0,2);
                            [slopesH{a,b,c,d,e},slopesPval{a,b,c,d,e}] = ttest2(slopesEnd{a,b,c,d,e},slopesEnd{a,control,1,d,e});
                            for g  = 2:length(indBlock)
                                if isempty(slopesN_extra{a,b,c,d,e}{g})==0
                                    slopes_extra_mean{a,b,c,d,e}{g}= mean(slopesN_extra{a,b,c,d,e}{g}(:,num_ind{a,b,c,d,e}>=g),2);
                                    slopes_extra_std{a,b,c,d,e}{g} = std(slopesN_extra{a,b,c,d,e}{g}(:,num_ind{a,b,c,d,e}>=g),0,2);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% plot slopes for each condition
for a = 1:length(induction);
    for b = 2:length(stim);
        for c = 1:length(intensity);
            for d = 1:length(position);
                for e = 1:length(drug);
                    if isempty(slices{a,control,1,d,e}) == 0;
                        if isempty(slices{a,b,c,d,e}) == 0;
                            figure;hold on
                            errorbar(t,slopesMean{a,b,c,d,e},slopesStd{a,b,c,d,e}/sqrt(size(slopesN{a,b,c,d,e},2)),...
                                '.','Color',stimcolor{b},'MarkerSize',30);
                            errorbar(t,slopesMean{a,1,1,d,e},slopesStd{a,1,1,d,e}/sqrt(size(slopesN{a,1,1,d,e},2)),...
                                '.','Color',stimcolor{1},'MarkerSize',30);
%                             title(strcat('TBS with ',stim{b},', ',num2str(intensity(c)),'V/m, ',position{d}));

                            % format figure
                            xSize = 15; ySize = 10;
                            set(gcf,'PaperUnits','centimeters')
                            xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
%                             set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
                            set(gcf,'Position',[0.5 0.5 xSize*50 ySize*50])
                            %set(gca,'TickDir','out','TickLength',[0.03 0.03]); % make the tick sizes bigger and outside
                            set(gcf, 'Color', 'w'); 
                            xlim([t(1) t(end)])
                            ylim([.6 1.8])
                            set(gca,'TickDir','out')
                            set(gca,'TickLength',[.03 .03],'FontSize',10,'FontWeight','bold','LineWidth',2)
                            xlabel('Time (min)','FontSize',30,'FontWeight','bold')
                            ylabel('Normalize fEPSP slope','FontSize',30,'FontWeight','bold')
                        end
                    end
                end
            end
        end
    end
end

%% plot slopes for conditions with multiple inductions
for a = 1:length(induction);
    for b = 2:length(stim);
        for c = 1:length(intensity);
            for d = 1:length(position);
                for e = 1:length(drug);
                    if isempty(slices{a,control,1,d,e}) == 0;
                        if isempty(slices{a,b,c,d,e}) == 0;
                            if isempty(slopes_extra_mean{a,b,c,d,e}{2})==0
                                figure;hold on
                                errorbar(t,slopesMean{a,b,c,d,e},slopesStd{a,b,c,d,e}/sqrt(size(slopesN{a,b,c,d,e},2)),'.','Color',stimcolor{b});
                                errorbar(t,slopesMean{a,1,1,d,e},slopesStd{a,1,1,d,e}/sqrt(size(slopesN{a,1,1,d,e},2)),'.','Color',stimcolor{1});
                                for g=2:length(indBlock)
                                    if isempty(slopes_extra_mean{a,b,c,d,e}{g})==0
                                        errorbar(t(end)+(g-2)*length(slopes_extra_mean{a,b,c,d,e}{g})+1:...
                                            t(end)+(g-1)*length(slopes_extra_mean{a,b,c,d,e}{g}),...
                                            slopes_extra_mean{a,b,c,d,e}{g},...
                                            slopes_extra_std{a,b,c,d,e}{g}/sqrt(length(slopes_extra_mean{a,b,c,d,e}{g})),...
                                            '.','Color',stimcolor{b});
                                    end
                                    if isempty(slopes_extra_mean{a,control,1,d,e}{g})==0
                                        errorbar(t(end)+(g-2)*length(slopes_extra_mean{a,control,1,d,e}{g})+1:...
                                            t(end)+(g-1)*length(slopes_extra_mean{a,control,1,d,e}{g}),...
                                            slopes_extra_mean{a,control,1,d,e}{g},...
                                            slopes_extra_std{a,control,1,d,e}{g}/sqrt(length(slopes_extra_mean{a,control,1,d,e}{g})),...
                                            '.','Color',stimcolor{control});
                                    end
                                end        
                                title(strcat('TBS with ',stim{b},', ',num2str(intensity(c)),'V/m, ',position{d}));
                                % format figure
                                xSize = 15; ySize = 10;
                                set(gcf,'PaperUnits','centimeters')
                                xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
                                set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
                                set(gcf,'Position',[0.5 0.5 xSize*50 ySize*50])
                                %set(gca,'TickDir','out','TickLength',[0.03 0.03]); % make the tick sizes bigger and outside
                                set(gcf, 'Color', 'w'); 
%                                 xlim([t(1) t(end)])
        %                         ylim([.8 2])
                                set(gca,'TickDir','out')
                                set(gca,'TickLength',[.03 .03],'FontSize',10,'FontWeight','bold','LineWidth',2)
                                xlabel('Time (min)','FontSize',30,'FontWeight','bold')
                                ylabel('Normalized slope','FontSize',30,'FontWeight','bold')
                            end
                        end
                    end
                end
            end
        end
    end
end

%% plasticity as a function dorsal ventral location
for d = 1:length(position)
    figure;hold on
    for a = 1:length(induction)
        for b = 1:length(stim)
            for c = 1:length(intensity)
                for e = 1:length(drug)
                    if isempty(slices{a,b,c,d,e}) == 0
                        plot(heights{a,b,c,d,e}(date_cut{a,b,c,d,e}),...
                            slopesEnd{a,b,c,d,e},...
                            '.','Color',stim_color{b},'MarkerSize',20);
                        xlabel('Position along dorsal ventral axis')
                        ylabel('Plasticity')
                        title(position{d})
                    end
                end
            end
        end
    end
end

%% plasticity as a function of baseline EPSP size
figure;hold on
for a = 1:length(induction);
    for b = 1:length(stim);
        for c = 1:length(intensity);
            for d = 2;%1:length(position);
                for e = 1:length(drug);
                    if isempty(slices{a,b,c,d,e}) == 0;
                        plot(mean(slopes{a,b,c,d,e}(1:tbase,date_cut{a,b,c,d,e}),1),slopesEnd{a,b,c,d,e},'x','Color',stimcolor{b},'MarkerSize',10)
%                         plot(mean(spikes{a,b,c,d,e}(1:tbase,:),1),slopesEnd{a,b,c,d,e},'x','Color',stimcolor{b},'MarkerSize',10)
                    end
                end
            end
        end
    end
end
xlabel('initial epsp')
ylabel('plasticity')

%%
figure;hold on
for a = 1:length(induction);
    for b = 1:length(stim);
        for c = 1:length(intensity);
            for d = 2;%1:length(position);
                for e = 1:length(drug);
                    if isempty(slices{a,b,c,d,e}) == 0;
                        plot(dates{a,b,c,d,e}(date_cut{a,b,c,d,e}),slopesEnd{a,b,c,d,e},'.','Color',stimcolor{b},'MarkerSize',10)
                    end
                end    
            end
        end
    end
end
