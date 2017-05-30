function replace = remove_noise(x,win,std_tol)
%==========================================================================
% remove samples that deviate from neighboring samples by a given tolerance

% x = original data

% win = window size to compare each sample against neighboring samples

% std_tol = allowed deviation from neighboring samples (in standard
% deviations)

% replace = logical output vector containing indeces of samples that exceed
% the deviation tolerance and are to be removed
%==========================================================================

replace = logical(zeros(length(x),1));
for g = 1:length(x)
    if g < win
        temp = x(1:g);
    else
        temp = x(g-win+1:g);
    end
    if abs(temp(end) - mean(temp(1:end-1)))/std(temp(1:end-1)) > std_tol
       replace(g) = logical(1);   
    end
end