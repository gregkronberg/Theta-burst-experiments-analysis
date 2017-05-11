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