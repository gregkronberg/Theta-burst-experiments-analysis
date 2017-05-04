function [new,direct1,direct2] = folder_diff(directory1,directory2,file_type)

% new = cell containing a list of different files between the two
% directories

% direct1 = cell containing list of all files in directory 1

% direct2 = cell containing list of all files in directory 2

% this function compares two directories and returns a cell containing file
% names of files that are different between the directories.  Note that the
% search is directional, i.e.files in directory1 are searched for in
% directory2.  Any unique files in directory2 are not identified.

% list contents of each directory as a structure
direct1_struct = dir(strcat(directory1,'*',file_type));
direct2_struct = dir(strcat(directory2,'*',file_type));

% convert lists to cells
direct1 = cell(length(direct1_struct),1);
for a = 1:length(direct1_struct)
    direct1{a} = direct1_struct(a).name;
end
direct2 = cell(length(direct2_struct),1);
for a = 1:length(direct2_struct)
    direct2{a} = direct2_struct(a).name;
end

cnt = 0 ;
new = cell(0);
% loop over files in first directory
for a = 1:length(direct1)
    % check if file exists in second directory
    if sum(strcmp(direct2,direct1{a}))==0                         
        cnt = cnt+1;
        % if file is unique, store its name
        new{cnt} = direct1{a};
    end
end
