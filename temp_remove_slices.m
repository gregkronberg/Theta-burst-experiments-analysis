%% reject slices
%==========================================================================
% notes
%==========================================================================
clear all
close all

%% file paths
%==========================================================================
 fpath_variables = 'D:\Google Drive\Work\Research Projects\Theta LTP\Matlab Variables\'; % variables
 
 
 load(strcat(fpath_variables,'slices.mat'));
 load(strcat(fpath_variables,'slopes.mat'));
 
 x = [0:59]'; % vector for fitting curves after induction
 % loop over conditions
 for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
                        cnt = 0;
                        for f = 1:length(slices{a,b,c,d,e})
% remove noise
% loop over points, if it is more than x std from mean replace with average
% of surrounding points
win = 10;
for g = 1:length(slopes{a,b,c,d,e}(f).slopes_norm(21:80))
    if g < win
        temp = slopes{a,b,c,d,e}(f).slopes_norm(21:20+g);
    else
        temp = slopes{a,b,c,d,e}(f).slopes_norm(20+g-win:20+g);
    end
    if abs(temp(end) - mean(temp))/std(temp) > 2.5
        slopes{a,b,c,d,e}(f).slopes_norm(20+g)  = mean(temp(1:end-1));
    end
end
y = abs(slopes{a,b,c,d,e}(f).slopes_norm(21:80)); % data to fit'
% linear fit
lin_fit = polyfit(x,y,1);
slopes{a,b,c,d,e}(f).lin_fit_coef = lin_fit';
slopes{a,b,c,d,e}(f).lin_fit = lin_fit(1)*x+lin_fit(2);
slopes{a,b,c,d,e}(f).lin_fit_error = y - slopes{a,b,c,d,e}(f).lin_fit;
slopes{a,b,c,d,e}(f).max = max(y);
% exponential fit
exp_fit = fit(x,y,'exp1');
slopes{a,b,c,d,e}(f).exp_fit = [exp_fit.a;exp_fit.b];
slopes{a,b,c,d,e}(f).exp_fit_post = exp_fit.a*exp(exp_fit.b*x);
slopes{a,b,c,d,e}(f).exp_fit_error = y - slopes{a,b,c,d,e}(f).exp_fit_post;
                        end
                    end
                end
            end
        end
    end
 end


 %% reject slices
  for a = 1:length(conditions{1})
    for b = 1:length(conditions{2})
        for c = 1:length(conditions{3})
            for d = 1:length(conditions{4})
                for e = 1:length(conditions{5})
                    if isempty(slices{a,b,c,d,e})==0
lin_coef =  [slopes{a,b,c,d,e}(:).lin_fit_coef];
exp_coef = [slopes{a,b,c,d,e}(:).exp_fit];
exp_error = mean([slopes{a,b,c,d,e}(:).exp_fit_error],1);
t80 = lin_coef(1,:)*90 + lin_coef(2,:);
slopes_norm{a,b,c,d,e} = [slopes{a,b,c,d,e}(:).slopes_norm];
slopes_end{a,b,c,d,e} = mean(slopes_norm{a,b,c,d,e}(end-9:end,:),1);
reject{a,b,c,d,e} = exp_coef(2,:) < -4e-3 | exp_coef(2,:)>0;
slopes_end_keep{a,b,c,d,e} = slopes_end{a,b,c,d,e}(~reject{a,b,c,d,e});
                    end
                end
            end
        end
    end
  end


%  dcs_color = {[0 0 0],[0 0 1],[1 0 0]};
%  cnt = 0 ;
%  for d = 1:length(conditions{4})
%      cnt = cnt+1;
%      figure(cnt);hold on
%      cnt = cnt+1;
%      figure(cnt);hold on
%      cnt = cnt+1;
%      figure(cnt);hold on
%   for a = 1:length(conditions{1})
%     for b = 1:length(conditions{2})
%         for c = 1:length(conditions{3})
%                 for e = 1:length(conditions{5})
%                     if isempty(slices{a,b,c,d,e})==0
%                         lin_coef =  [slopes{a,b,c,d,e}(:).lin_fit_coef];
%                         exp_coef = [slopes{a,b,c,d,e}(:).exp_fit];
%                         max{a,b,c,d,e} = [slopes{a,b,c,d,e}(:).max];
%                         coef_norm = lin_coef(1,:)./max{a,b,c,d,e};
%                         t80 = lin_coef(1,:)*90 + lin_coef(2,:);
%                         slopes_norm{a,b,c,d,e} = [slopes{a,b,c,d,e}(:).slopes_norm];
%                         slopes_end{a,b,c,d,e} = mean(slopes_norm{a,b,c,d,e}(end-9:end,:),1);
%                         figure(cnt-2)
%                         plot(coef_norm,slopes_end{a,b,c,d,e},'.','Color',dcs_color{b},'MarkerSize',25)
%                         figure(cnt-1)
%                         plot(t80,slopes_end{a,b,c,d,e},'.','Color',dcs_color{b},'MarkerSize',25)
%                         reject{a,b,c,d,e} = t80<1 | exp_coef(2,:)>0;
%                         slopes_end_keep{a,b,c,d,e} = slopes_end{a,b,c,d,e}(~reject{a,b,c,d,e});
%                         figure(cnt)
%                         plot(b*ones(length(slopes_end_keep{a,b,c,d,e}),1),slopes_end_keep{a,b,c,d,e},'.',...
%                             'MarkerSize',25,'Color',dcs_color{b})
% %                         figure(cnt)
% %                         plot(b*ones(length(max{a,b,c,d,e}),1),max{a,b,c,d,e},'.',...
% %                             'MarkerSize',25,'Color',dcs_color{b})
%                     end
%                 end
%             end
%         end
%     end
%   end
  
  
 