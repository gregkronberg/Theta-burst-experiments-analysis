function new = folder_diff(directory1,directory2,file_type)

% this function compares two directories and returns a cell containing file
% names of files that are different between the directories.  Note that the
% search is direction, i.e.files in directory1 are searched for in
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
new = cell(0)
for a = 1:length(direct1)
    if sum(strcmp(direct2,direct1{a}))==0                              % check if current slice is new
        cnt = cnt+1;
        new{cnt} = direct1{a};
    end
end
