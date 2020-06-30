

% promising behaviors: wet dog shake, groom, scratch, rear
promising_medium_fine_behaviors = [];
for k = 1:numel(promising_behaviors)
    promising_medium_fine_behaviors = [ promising_medium_fine_behaviors find(contains(medium_fine_behav,promising_behaviors{k}))];
end
promising_medium_fine_behaviors = unique(promising_medium_fine_behaviors);
promising_medium_fine_behaviors_index = find(ismember(fine_labels_here,promising_medium_fine_behaviors));
keep = intersect(promising_medium_fine_behaviors_index,find(cellfun('length',fine_behav_accel_ts_ind) > min_iter & cellfun(@range,fine_behav_accel_ts_ind) > min_range));

fine_behav_accel_traces = fine_behav_accel_traces(keep);
fine_behav_accel_goodframe_ind = fine_behav_accel_goodframe_ind(keep);
fine_behav_accel_ts_ind = fine_behav_accel_ts_ind(keep);
fine_behav_accel_fine_behav_ind = fine_behav_accel_fine_behav_ind(keep);
fine_behav_accel_ts_vector = fine_behav_accel_ts_vector(keep);