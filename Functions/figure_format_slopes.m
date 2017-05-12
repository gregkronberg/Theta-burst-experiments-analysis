%% figure format

xSize = 15; ySize = 10;
set(gcf,'PaperUnits','centimeters')
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
%                             set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[0.5 0.5 xSize*50 ySize*50])
%set(gca,'TickDir','out','TickLength',[0.03 0.03]); % make the tick sizes bigger and outside
set(gcf, 'Color', 'w'); 
xlim([t(1) t(end)])
% ylim([.6 1.8])
set(gca,'TickDir','out')
set(gca,'TickLength',[.03 .03],'FontSize',10,'FontWeight','bold','LineWidth',2)