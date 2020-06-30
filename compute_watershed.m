function [ws,ws_ind,density_plot_smooth] = compute_watershed(x1,x2)

% compute the density plot
xaxis = linspace(min(x1),max(x1),100);
yaxis = linspace(min(x2),max(x2),100);
density_plot = zeros(numel(xaxis),numel(yaxis));
for k = 1:numel(xaxis)-1
    indx = find(x1 > xaxis(k) & x1 <= xaxis(k+1));
    for j = 1:numel(yaxis)-1
        indy = find(x2 > yaxis(j) & x2 <= yaxis(j+1));
        density_plot(numel(yaxis)-j+1,k) = numel(intersect(indx,indy));
    end
end

density_plot_smooth = imgaussfilt(density_plot,[3 3]);

%% determine whether identified features correspond at all to mocap behaviors

ws = watershed(-density_plot_smooth);
ws(density_plot==0) = 0;

ws_ind = nan(numel(x1),1);
for k = 1:numel(xaxis)-1
    indx = find(x1 >= xaxis(k) & x1 <= xaxis(k+1));
    for j = 1:numel(yaxis)-1
        indy = find(x2 >= yaxis(j) & x2 <= yaxis(j+1));
        ws_ind(intersect(indx,indy)) = ws(numel(yaxis)-j+1,k);
    end
end

ws_nonzero = find(ws_ind~=0);
for k = 1:numel(x1)
    if ws_ind(k) == 0
        index = knnsearch([x1(ws_nonzero) x2(ws_nonzero)],[x1(k) x2(k)]);
        ws_ind(k) = ws_ind(ws_nonzero(index));
    end
end


return