%% plot information about the template
close all
motif_ind_k = motif_ind_all{k};
num_trials = size(motif_ind_k,2);
template = squeeze(almost_final_templates(k,:,:))';
traces = identified_traces_all{k};
identified_trace_time_k = identified_trace_time_all{k};

%%%%%% find the behavioral labels associated with this trace %%%%%%
index = ts_ds(motif_ind_k)+snip_length; % the index of the trace in ts
[behav_index_k1, distance_to_label_k] = knnsearch(good_frames,index');
too_far = find(distance_to_label_k>25);
behav_index_k = max(fine_labels_good_frames(behav_index_k1),1);
behav_label_k = medium_fine_behav(behav_index_k);
if numel(too_far) > 0
    behav_label_k(too_far) = {'NotLabeled'};
end
[unique_labels,~,c] = unique(behav_label_k,'stable');
num_unique = numel(unique_labels);
colorlabels = parula(num_unique);
[label_freq,~] = hist(c,1:num_unique);
[unique_label_n,unique_label_ind] = sort(label_freq,'descend');

% find the most common medium fine label
% find the other times this happened in the recording

most_common_label = unique_labels{mode(c)};
common_label_fine_ind = find(strcmp(medium_fine_behav,unique_labels{mode(c)}));
common_label_ind = good_frames(find(ismember(fine_labels_good_frames,common_label_fine_ind)));
common_template = find(strcmp(behav_label_k,most_common_label));
notcommon_template = setdiff(1:numel(motif_ind_k),find(strcmp(behav_label_k,most_common_label)));


figure(1)
subplot(3,5,1)
plot(Zvals(:,1),Zvals(:,2),'.k');
hold on
medZval1 = median(Zvals(behav_index_k1(:),1));
medZval2 = median(Zvals(behav_index_k1(:),2));
plot(median(Zvals(behav_index_k1(:),1)),median(Zvals(behav_index_k1(:),2)),'*r')
% plot the traces for the template

for i = 1:num_unique
    labeltype = unique_labels(unique_label_ind(i));
    if ~strcmp(labeltype,'NotLabeled')
    index_i = find(strcmp(behav_label_k,unique_labels(unique_label_ind(i))));
    figure(1)
    subplot(3,5,1)
    plot(Zvals(behav_index_k1(index_i),1),Zvals(behav_index_k1(index_i),2),'*','color',colorlabels(i,:))
    axis off
    end
end
unique_vec = 1:num_unique;
for i = 1:min(num_unique,6)
    figure(1)
    subplot(3,5,2)
    h = bar(i,unique_label_n(i));
    set(h,'FaceColor',colorlabels(unique_label_ind(i),:))
    hold on
end
figure(1)
subplot(3,5,2)
xticks(1:numel(unique_labels))
xticklabels(unique_labels(unique_label_ind))
box off
xtickangle(30)
title([num2str(numel(common_template)),'/',num2str(numel(motif_ind_k))])

% plot the traces and template

for d = 1:3
    subplot(3,5,2+d)
    plot(squeeze(traces(:,:,d))','color',0.5*ones(1,3));
    hold on
    plot(template(d,:),'k','linewidth',2);
    hold off
    box off
    title(almost_final_medium_fine_behav_name{k})
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the neural responses that are the most modulated
spike_mat_k = nan(num_good_units,snip_length,numel(motif_ind_k));
for j = 1:numel(motif_ind_k)
    spike_mat_k(:,:,j) = spike_matrix(:,round(identified_trace_time_k(j,:)));
end

spikes_per_cell = sum(sum(spike_mat_k,2),3);
enough_spikes = find(spikes_per_cell>0.25*num_trials);
tuning_stability = nan(num_good_units,1);
modulation = nan(num_good_units,1);
p_mod = nan(num_good_units,1);

for i = enough_spikes'
    tuning_strength_all = compute_trace_similarity(squeeze(spike_mat_k(i,:,:))');
    tuning_stability(i) = nanmedian(tuning_strength_all);
    [modulation(i),p_mod(i)] = compute_modulation(squeeze(spike_mat_k(i,:,:))');
end

% plot the top 10
[val,ind] = sort(tuning_stability,'descend');
ind(isnan(val)) = [];
num_example = 10;
count = 0;
for i = ind(1:num_example)'
    count = count + 1;
    r = smoothdata(sum(squeeze(spike_mat_k(i,:,:)),2),'gaussian',10)/numel(motif_ind_k);
    r_sem = smoothdata(std(squeeze(spike_mat_k(i,:,:)),[],2)./sqrt(numel(motif_ind_k)),'gaussian',10);
    
    figure(1)
    subplot(3,5,5+count)
    for h = 1:numel(motif_ind_k)
        
        spikes = find(spike_mat_k(i,:,h));
        plot([spikes;spikes],(h-1)+[ones(size(spikes));zeros(size(spikes))],'k-','linewidth',2)
        hold on
  
    end
    
    scalefactor = numel(motif_ind_k)/5;
    scaled_rsem = (r_sem - min(r))/range(r)*scalefactor;
    scaled_r = (r - min(r))/range(r)*scalefactor;
    fill([(1:snip_length)'; flipud((1:snip_length)')],h+[scaled_r-scaled_rsem; flipud(scaled_r+scaled_rsem)],0.9*ones(1,3),'linestyle','none')
    plot(1:snip_length,h+scaled_r,'k','linewidth',2)
    box off
    axis tight
    
    maxval = max(scaled_r+scaled_rsem)+h;
    for j = 1:3
        template_j =(template(j,:)  - min(template(j,:)))/range(template(j,:))*scalefactor;
        plot(1:snip_length,template_j+maxval+scalefactor*(j-1),'k','linewidth',2)
    end
    hold off
    title(num2str(tuning_stability(i)))
end

