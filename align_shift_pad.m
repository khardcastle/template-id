function [aligned_X_mat,template,delay_all] = align_shift_pad(X,iteration,maxlag)

% basic idea:
% align just by shifting trials in X
% X = matrix of M trials, of length T, with dim D (M x T x D)

% get dimensions
M = size(X,1); % trial number
T = size(X,2); % trial length
D = size(X,3); % dim number

X_inside = X(:,maxlag+1:T-maxlag,:);
T_inside = size(X_inside,2);

% zscore X
X = zscore(X,[],2);

% build Xcell = M cells, each D x T
Xcell = cell(M,1);
X_4corr = nan(M,T_inside*D);
for k = 1:M
    Xcell{k} = squeeze(X(k,:,:))';
    X_4corr(k,:) = reshape((X_inside(k,:,:)),1,D*T_inside);
end


% initialize
aligned_X = Xcell;
aligned_X_4corr = X_4corr;

% find two the closest, choose one as initial template
c = corr(X_4corr');
[val,~] = max(triu(c,1));
[~,row] = max(val);
template = Xcell{row}(:,maxlag+1:T-maxlag);

total_mean_dist = [];
delay_all = zeros(M,1);

for j = 1:iteration
    
    if j == 1
        rem_trials = setdiff(1:M,row);
        maxm = M-1;
    else
        rem_trials = 1:M;
        maxm = M;
    end
    count = 0;
    for k = 1:maxm
        
        % find next sequence closest to template
        [distance,index1] = pdist2(X_4corr(rem_trials,:),reshape(template',1,D*T_inside),'correlation','Smallest',1);
        index = rem_trials(index1);
        
        % pull out trial
        trial = Xcell{index};
        
        % align trial to template
        [trial_k,delay] = alignmatrix_pad(template,trial,maxlag);

        % put back in units of template
        aligned_X_4corr(k,:) = reshape(trial_k',1,D*T_inside);
        
        % save aligned trial
        aligned_X{index} = trial_k;
        delay_all(index) = delay;
        
        % re-compute template, but only if the current trial is correlated
        % enough with the template
        correlation_k = corr(reshape(trial_k',1,D*T_inside)',reshape(template',1,D*T_inside)');
        if correlation_k > 0
            count = count + 1;
            template = (template.*count + trial_k)./(count+1);
        end
        
        % take away remaining trials
        rem_trials(index1) = [];
        
    end
    
end


% re-arrange
aligned_X_mat = nan(size(X_inside));
for k = 1:M
    aligned_X_mat(k,:,:) = aligned_X{k}';
end

% plot the results
%{
figure(1)
for k = 1:3
    
    subplot(2,3,k)
    plot(X_inside(:,:,k)','color',[0.5 0.5 0.5]);
    
    subplot(2,3,k+3)
    plot(aligned_X_mat(:,:,k)','color',[0.5 0.5 0.5]);
    hold on
    plot(template(k,:),'r','linewidth',0.5)
end
%}


return