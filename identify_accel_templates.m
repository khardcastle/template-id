%% Load ephys and motion capture labels

%% clear everything
clear all; close all; clc

% add paths
addpath('C:\Users\kiahh\Dropbox\acceler-analysis\toolbox')
addpath('C:\Users\kiahh\Dropbox\acceler-analysis\mocap_based_template_detection')


%% load the data
disp('loading data')

% set the data to download
% right now the data is a bit dispersed... not ideal. can prob consolidate
rathere = 'JDM36';
dayhere = 1;
ephys_path = 'Y:\Jesse\Data\Dropbox_curated_sharefolders\ephys_data_diego\';
analysis_struct_path = strcat('Y:\Kiah\mocap_data_from_jesse\',rathere,'\day',num2str(dayhere),'\');
ephys_file = ['diego_ephys_files_rat_',rathere,'_day_',num2str(dayhere)];
downsample_ratio = 3; % original is 300 Hz. 3 = downsample to 100Hz
dim = 3;

% load the data
load_mocap_ephys_data

% set some initial parameters:
snip_length = fs/downsample_ratio; % 1 second
buffer = snip_length;
signal = accel_3d_ds_z; % the signal used to find accel traces

%% pull out all the accelerometer traces
% this finds the traces for each "super fine" behavior
% this does not do any alignment
disp('pulling out accel traces')
[fine_behav_accel_traces,fine_behav_accel_goodframe_ind,...
    fine_behav_accel_ts_ind,fine_behav_accel_fine_behav_ind,fine_behav_accel_ts_vector] = ...
    pull_out_accel_traces_fcn(fine_labels_here,fine_labels_good_frames,good_frames,...
    ts_ds,signal,snip_length,buffer);

%% refine the traces
% this will reset fine_behav_accel_traces, etc.. (from above)
% refine: look for fine behaviors that occur > 10 times
%                                 that occur with range > 1 hour
%                                 are from subset of promising behaviors
disp('refining accel traces')

min_iter = 20;
promising_behaviors = {'Shake','shake','groom','Groom','scratch','Scratch','rear','Rear','Squat','squat','Wipes','wipes','Head','Drinking'};
min_range = fs_ds*60*60;
refine_traces_for_mocap_template_id

%% align the traces for each super-fine behavior
% this is a time-intensive step
disp('aligning accel traces')
fine_behav_accel_traces_aligned = cell(size(fine_behav_accel_traces));
time_aligned = cell(size(fine_behav_accel_traces));
template_accel_traces = nan(numel(keep),snip_length,3);
iteration = 2;
maxlag = snip_length;
scalevec = linspace(0.8,1.25,15);
num_templates = numel(fine_behav_accel_goodframe_ind);
for k = 1:numel(keep)
    tic
    size(fine_behav_accel_traces{k},1)
    [fine_behav_accel_traces_aligned{k},template_k,~,~,time_aligned{k}] = align_ltw_pad(fine_behav_accel_traces{k},iteration,maxlag,scalevec);
    
    % save the aligned accelerometer templates (averaged traces)
    template_accel_traces(k,:,:) = template_k';
    k
    toc
end

%% find all the templates that have well-correlated traces
correlation_threshold = 0.7;
high_corr = nan(num_templates,1);
medium_fine_behav_name = cell(num_templates,1);
watershed_label = nan(num_templates,1);
corr_to_temp = cell(num_templates,1);
for k = 1:num_templates
    
    watershed_label(k) = fine_behav_accel_fine_behav_ind{k}(1);
    medium_fine_behav_name_k = medium_fine_behav{watershed_label(k)};
    medium_fine_behav_name{k} = medium_fine_behav_name_k;
    corr_to_temp_k = nan(size(fine_behav_accel_traces_aligned{k},1),size(fine_behav_accel_traces_aligned{k},3));
    for j = 1:3
        corr_to_temp_k(:,j) = corr(squeeze(fine_behav_accel_traces_aligned{k}(:,:,j))',squeeze(template_accel_traces(k,:,j))');
    end
    
    corr_to_temp{k} = mean(corr_to_temp_k,2);
    high_corr(k) = numel(find(corr_to_temp{k}>correlation_threshold));
    k
end

% keep only the ones that have enough traces with a high correlation
keep = find(high_corr>10);
almost_final_traces_aligned = fine_behav_accel_traces_aligned(keep);
almost_final_accel_ts_vector = fine_behav_accel_ts_vector(keep);
almost_final_templates = template_accel_traces(keep,:,:);
almost_final_medium_fine_behav_name = medium_fine_behav_name(keep);
almost_final_watershed_label = watershed_label(keep);
almost_final_time_aligned = time_aligned(keep);
corr_to_template = corr_to_temp(keep);
almost_final_ts_ind = fine_behav_accel_ts_ind(keep);
num_behavs = numel(keep);

%% second pass - doing this for templates, on the whole session

candidate_threshold = 0.65;
keep = [];
min_num = 30;
min_range = fs_ds*60*60*3;
%for k = 1:num_behavs
randvec = randperm(num_behavs);
for k = randvec
    
    % compute the correlation bw the motif and the full signal
    template = squeeze(almost_final_templates(k,:,:))';
    [motif_ind_k,identified_traces_k,scales_k,mean_corr_mp_k,mean_corr_scale_k,...
        identified_trace_time_k] = find_motifs_in_signal(template,signal,scalevec,candidate_threshold);

    if numel(motif_ind_k)>min_num && range(motif_ind_k)>min_range
        
        motif_ind_all{k} = motif_ind_k;
        identified_traces_all{k} = identified_traces_k;
        identified_trace_time_all{k} = identified_trace_time_k;
        
        % plot for the behavior
        %compute_plot_template_signal_evaluation
        
        % plot the spikes
        %keyboard
        
    end
    k
end

%%%%%%% extra code to do some filtering on the templates

%% filter some of the resulting templates to get a cleaner list

% compute the median spread on the z-map
median_distance_all = nan(num_behavs,1);
for k = 1:num_behavs
    if numel(motif_ind_all{k}) > 0
        median_distance_all(k) = compute_dist_zvals(Zvals,motif_ind_all{k},snip_length,good_frames,ts_ds);
    end
end
compact = find(median_distance_all < 6);

final_template_ind = intersect(find(~cellfun(@isempty,motif_ind_all)),compact);
final_templates = almost_final_templates(final_template_ind,:,:);
final_motif_ind = motif_ind_all(final_template_ind);
final_motif_ind = cellfun(@transpose,final_motif_ind,'UniformOutput',false);
final_identified_traces = identified_traces_all(final_template_ind);
final_trace_time = identified_trace_time_all(final_template_ind);
num_final_behav = numel(final_template_ind);

% get combined lists of traces, trace indexes, motif indices
all_motif_ind_together = cell2mat(final_motif_ind');
all_traces_together = cell2mat(final_identified_traces');
all_traces_together_ind = nan(size(all_traces_together,1),1);
count = 0;
for k = 1:numel(final_template_ind)
    num_traces_k = size(final_identified_traces{k},1);
    count = count + num_traces_k;
    countvec = count(end)-num_traces_k+1:count(end);
    all_traces_together_ind(countvec) = k;
end

%% CLUSTER THE TEMPLATES

index = ts_ds(all_motif_ind_together)+snip_length; % the index of the trace in ts
[behav_index_k1, distance_to_label_k] = knnsearch(good_frames,index');
close_enough = find(distance_to_label_k<25);

behav_index_k1 = behav_index_k1(close_enough);
all_traces_together_close = all_traces_together(close_enough,:,:);
all_motif_ind_together_close = all_motif_ind_together(close_enough);
all_traces_together_ind_close = all_traces_together_ind(close_enough);

% compute median for each behavior
median_Zval_behaviors = nan(num_final_behav,2);
for k = 1:num_final_behav
    median_Zval_behaviors(k,:) = median(Zvals(behav_index_k1(all_traces_together_ind_close==k),:));
end

[ws,ws_ind,density_plot] = compute_watershed(median_Zval_behaviors(:,1),median_Zval_behaviors(:,2));
numclusters = numel(unique(ws_ind));

behav_cluster_ind = nan(size(all_traces_together_ind));
for k = 1:num_final_behav
    behav_cluster_ind(all_traces_together_ind==k) = ws_ind(k);
end


%% plot the traces for each of the clusters
for k = 1:numclusters
    
    behavs_k = find(ws_ind==k);
    
    figure(2)
    plot(Zvals(:,1),Zvals(:,2),'.k')
    hold on
    scatter(median_Zval_behaviors(:,1),median_Zval_behaviors(:,2),25,ws_ind,'filled')
    plot(median_Zval_behaviors(behavs_k,1),median_Zval_behaviors(behavs_k,2),'*r')
    hold off
    
    figure(1)
    for j = 1:numel(behavs_k)
        subplot(2,ceil(numel(behavs_k)/2),j)
        for i = 1:3
            % plot the template and the traces
            plot((1:snip_length)+snip_length*(i-1)+20,squeeze(final_identified_traces{behavs_k(j)}(:,:,i)),'color',0.5*ones(1,3))
            hold on
            plot((1:snip_length)+snip_length*(i-1)+20,final_templates(behavs_k(j),:,i),'k','linewidth',1.5)
            title(almost_final_medium_fine_behav_name{compact(behavs_k(j))})
        end
    end
    pause
    close all
end
