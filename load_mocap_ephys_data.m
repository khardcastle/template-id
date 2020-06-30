%% load the accelerometer data

% Structure that contains accelerometer data
load([ephys_path,ephys_file])

% the full 3d accelerometer
accel_3d = agg_unitstructure.acc_full; % this is reported in volts, corresponds to force
T = size(accel_3d,2); % time units are the same as in the motion capture data - 300 hz, total (18.5 hours)

%% lightly pre-process the accelerometer data
fs = 300;
ts = 1:T;
filter = gausswin(5)./sum(gausswin(5));
accel_3d_smooth = nan(size(accel_3d));
for k = 1:3
    accel_3d_smooth(k,:) = conv(accel_3d(k,:),filter,'same');
end
fs_ds = fs/downsample_ratio;
ts_ds = ts(1:downsample_ratio:end);
accel_3d_ds = accel_3d_smooth(:,ts_ds);
accel_3d_ds_z = zscore(accel_3d_ds,[],2);
T_ds = size(accel_3d_ds,2); % time units are the same as in the motion capture data - 300 hz, total (18.5 hours

%% load the ephys data
%% get the units that fire enough during this session
all_units = agg_unitstructure.units;
num_all_units = numel(all_units);
range_k = nan(num_all_units,1);
num_spikes = nan(num_all_units,1);
for k = 1:num_all_units
    if numel(all_units{k}{1}) > 1
        range_k(k) = range(all_units{k}{1});
        num_spikes(k) = numel(all_units{k}{1});
    end
end

good_units = intersect(find(range_k > T/2),find(num_spikes > T/300*0.1));
num_good_units = numel(good_units);
%% create spike matrix, downsampled to match accel data

spike_matrix = nan(num_good_units,T_ds);
bin_edges = linspace(1,T,T_ds+1);
for k = 1:num_good_units
    spike_matrix(k,:) = histcounts(all_units{good_units(k)}{1},bin_edges);
end
spike_matrixs_sparse = sparse(spike_matrix);

%% load behavioral annotation
% Behavioral Annotation
% This contains information about behavioral labels, as the fields
% 'annot_reordered', which has the continuous time categorization of the
% behavior, or 'Zvalues', which has an analog representation in the tsne
% space. These labels are only made at 6 Hz, and in frames in which the
% animal is (1) moving, and (2) all markers are tracked. So the frame
% correspondence here is determined by the field 'frameswithgoodtracking'.
% The annotations for each of these behaviors is given in the
% 'clusternames' field. There are different levels of coarseness to the
% clusternames that may be useful in the future.
load(strcat(analysis_struct_path,'analysisstruct_wnames_rat',rathere,'_days_',num2str(dayhere),'.mat'));

good_frames = analysisstruct.frames_with_good_tracking{1}; % frames with good tracking + animal is moving
Zvals = analysisstruct.zValues; % z-values in t-sne space

% pull out the three labels: watershed region, med-fine labels, coarse
% labels
fine_labels_good_frames = analysisstruct.annot_reordered{1}; % behavioral labels
fine_labels_here = unique(fine_labels_good_frames); % list of the possible fine behav frames
num_fine_labels = numel(fine_labels_here);
coarse_behav = analysisstruct.clusternames; % 13 behaviors (coarse)
medium_fine_behav = analysisstruct.clusternames_fine;