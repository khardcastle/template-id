function [index_final] = reject_neighbors(x,y,m)

% x is the index
% y is the value (higher is better)
% m is the snippet length

% find the smallest ones that are
[sort_y,sort_y_ind] = sort(y,'descend');

sort_y_ind_loop = sort_y_ind; % changing list of mp index
sort_y_loop = sort_y; % changing list of mp values

index = [];

while ~isempty(sort_y_ind_loop)
    
    % candidates ind
    candidate_index = sort_y_ind_loop(1);
    
    % put that in the final list
    index = [index; candidate_index];
    
    % delete all of the ones within m
    delete = find(abs(x(sort_y_ind_loop) - x(candidate_index)) < m);
    sort_y_ind_loop(delete) = [];
    sort_y_loop(delete) = [];
    
end

index_final = x(index);

return