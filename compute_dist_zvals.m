function [median_distance] = compute_dist_zvals(Zvals,motif_ind_k,snip_length,good_frames,ts_ds)


index = ts_ds(motif_ind_k)+snip_length; % the index of the trace in ts
[behav_index_k1, distance_to_label_k] = knnsearch(good_frames,index');
behav_index_k1(distance_to_label_k>100) = [];

medZval1 = median(Zvals(behav_index_k1(:),1));
medZval2 = median(Zvals(behav_index_k1(:),2));

median_distance = median(sqrt((Zvals(behav_index_k1(:),1) - medZval1).^2 + (Zvals(behav_index_k1(:),2) - medZval2).^2));

return