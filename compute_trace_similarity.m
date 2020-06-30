function [trace_similarity] = compute_trace_similarity(spikemat)

num_trials = size(spikemat,1);
r = smoothdata(sum(spikemat),'gaussian',10)/num_trials;
smoothed_spikemat = smoothdata(spikemat,2,'gaussian',10);

trace_similarity = corr(r',smoothed_spikemat');
trace_similarity(isnan(trace_similarity)) = 0;
return