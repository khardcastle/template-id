function [shifted_Y,delay,maxcorr,index] = alignmatrix_pad(X,Y,maxlag)

d = size(X,1);
l = size(X,2);
numlags = maxlag*2+1;
allcorr = nan(d,numlags);

% create a Y matrix to correlate with X
Ymat = nan(d,l,numlags);
index_mat = nan(numlags,l);
for k = 1:numlags
    Ymat(:,:,k) = Y(:,k:k+l-1);
    index_mat(k,:) = k:k+l-1;
end

for k = 1:d
    allcorr(k,:) = 1-pdist2(X(k,:),squeeze(Ymat(k,:,:))','correlation');
end
avgcorr = nanmean(allcorr);
lags = -maxlag:maxlag;
[maxcorr,val] = max(avgcorr);
delay = lags(val);
index = index_mat(val,:);
shifted_Y = Ymat(:,:,val);


return
