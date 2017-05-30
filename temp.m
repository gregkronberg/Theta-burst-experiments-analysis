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

a = 1;
b = 1;
c = 1;
d = 1;
e = 1;
f = 4;
load(strcat(fpath_processed,slices{a,b,c,d,e}(f).name))

baseD_norm = (baseD  - ones(size(baseD,1),1)*mean(baseD,1))./(ones(size(baseD,1),1)*var(baseD,[],1));
baseS_norm = (baseS  - ones(size(baseS,1),1)*mean(baseS,1))./(ones(size(baseS,1),1)*var(baseS,[],1));
figure;hold on
plot(baseD_norm(20:200,indBlock(1):indBlock(1)+60),baseS_norm(20:200,indBlock(1):indBlock(1)+60),'.')
x  = ginput(4);
