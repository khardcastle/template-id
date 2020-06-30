function [fine_behav_accel_traces,fine_behav_accel_goodframe_ind,...
    fine_behav_accel_ts_ind,fine_behav_accel_fine_behav_ind,fine_behav_accel_ts_vector] = ...
    pull_out_accel_traces_fcn(fine_labels_here,fine_labels_good_frames,good_frames,ts_ds,...
    signal,snip_length,buffer)


num_fine_labels = numel(fine_labels_here);
ds = ts_ds(2)-ts_ds(1);
%% for a list of potential timepoints, pull out the accelerometer traces

fine_behav_accel_traces = cell(num_fine_labels,1);
fine_behav_accel_goodframe_ind = cell(num_fine_labels,1);
fine_behav_accel_ts_ind = cell(num_fine_labels,1);
fine_behav_accel_fine_behav_ind = cell(num_fine_labels,1);
fine_behav_accel_ts_vector = cell(num_fine_labels,1);

for k = 1:num_fine_labels
    
    % find all the places where this fine behavior occurred
    behav_k = find(fine_labels_good_frames == fine_labels_here(k));
    behav_k_ind = good_frames(behav_k); %indexed into ts
    behav_k_ind_ts_ds = ceil(behav_k_ind/ds);
    
    behav_k_accel_mat = nan(numel(behav_k),snip_length+buffer*2,3);
    reject = [];
    index_j_mat = [];
    for j = 1:numel(behav_k)
        
        % pull out the snip length, plus maxlag buffer on each side
        index_j = behav_k_ind_ts_ds(j)-snip_length/2+1-buffer:behav_k_ind_ts_ds(j)+snip_length/2+buffer;
        if min(index_j) > 1 && max(index_j) < numel(ts_ds)
            behav_k_accel_mat(j,1:snip_length+buffer*2,:) = signal(:,index_j)';
            index_j_mat = [index_j_mat; index_j];
        else
            reject = [reject j];
        end
    end
    fine_behav_accel_goodframe_ind{k} = behav_k(setdiff(1:numel(behav_k),reject));
    fine_behav_accel_fine_behav_ind{k} = fine_labels_here(k)*ones(numel(behav_k)-numel(reject),1);
    behav_k_accel_mat(reject,:,:) = [];
    fine_behav_accel_ts_ind{k} = behav_k_ind(setdiff(1:numel(behav_k),reject));
    fine_behav_accel_traces{k} = behav_k_accel_mat;
    fine_behav_accel_ts_vector{k} = index_j_mat;
    k
end

return