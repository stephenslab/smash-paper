function m = mad2(x)

% mad2: mean absolute deviation estimation function,
% where x is suposed to be a vector column

med= median(x);
m = median(abs(x-med))/0.6745; 