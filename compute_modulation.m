function [modulation_data,p_mod] = compute_modulation(spikemat)

num_trials = size(spikemat,1);
r = smoothdata(sum(spikemat),'gaussian',10)/num_trials;
modulation_data = range(r)/nanmean(r);

iter = 100;
modulation_shuffle = nan(iter,1);
for i = 1:iter
    id = randi(num_trials,num_trials,1);
    spikemat_i = cell2mat(arrayfun(@(x) circshift(spikemat(x,:),[1 id(x)]),(1:numel(id))','un',0));
    r_i = smoothdata(sum(spikemat_i),'gaussian',10)/num_trials;
    modulation_shuffle(i) = range(r_i)/nanmean(r_i);
end

p_mod = sum(modulation_shuffle > modulation_data)/iter;

return