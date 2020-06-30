function [motif_ind,identified_traces,scales,mean_corr_mp_k,mean_corr_scale,...
    identified_trace_ind] = find_motifs_in_signal(template,signal,scalevec,candidate_threshold)


mean_corr_mp_k = zeros(1,length(signal));
mean_corr_scale = zeros(1,length(signal));
D = size(template,1);
m = size(template,2);
m_scalevec =  round(m*scalevec);

for k = 1:numel(scalevec)
    
    % convolute the motif with the signal
    m_k = m_scalevec(k);
    template_k = nan(D,m_k);
    for j = 1:D
        template_k(j,:) = interp1(1:m,template(j,:),linspace(1,m,m_k));
    end
    
    corr_mp_k = nan(D,size(signal,2));
    for j = 1:D
        d1 = real(mass(signal(j,:),template_k(j,:)));
        d1 = [d1 20*ones(1,length(signal)-length(d1))];
        corr_mp_k(j,:) = 1-d1.^2/(2*m_k);

    end
    
    higher  = find(mean(corr_mp_k)>mean_corr_mp_k);
    mean_corr_scale(higher) = k;
    mean_corr_mp_k = max(mean(corr_mp_k),mean_corr_mp_k);
    
end

above_thresh_k = find(mean_corr_mp_k>candidate_threshold);

% find the threshold crossings, reject neigbors
[motif_ind] = reject_neighbors(above_thresh_k,mean_corr_mp_k(above_thresh_k),m);


motif_ind(motif_ind>length(signal)-max(m_scalevec)) = [];
scales = scalevec(mean_corr_scale(motif_ind));
m_k_motif = m_scalevec(mean_corr_scale(motif_ind));
identified_traces = nan(numel(motif_ind),m,D);
identified_trace_ind = nan(numel(motif_ind),m);
for k = 1:numel(motif_ind)
    
    temp = zscore(signal(:,motif_ind(k):motif_ind(k)+m_k_motif(k)-1)');
    identified_trace_ind(k,:) = interp1(1:m_k_motif(k),motif_ind(k):motif_ind(k)+m_k_motif(k)-1,linspace(1,m_k_motif(k),m));
    interp1(1:m_k_motif(k),temp(:,j),linspace(1,m_k_motif(k),m))
    
    for j = 1:D
        identified_traces(k,:,j) = interp1(1:m_k_motif(k),temp(:,j),linspace(1,m_k_motif(k),m));
    end
    
end

% visualize
%{
    figure(1)
    for j = 1:D
        subplot(1,3,j)
        plot(identified_traces(:,:,j)','color',0.5*ones(1,3))
        hold on
        plot(zscore(template(j,:)),'k','linewidth',2)
        hold off
    end
%}


return