%% Load in the data
clear all

% Path to LabDir on Holyoke
clusterPath = 'Y:/';
% clusterPath = '/n/holylfs02/LABS/olveczky_lab/';

% Get dependencies
addpath(genpath(fullfile(clusterPath,'Diego/code/Label3D')))

% Load data (Might take a second over network connection)
predFolder = fullfile(clusterPath, 'Diego/tdata/dannce/paw_labeling/JDM31/DANNCE/predict12/');
data = load_multi_gpu_data(predFolder);
skeleton = load('jesse_skeleton_paws.mat');

%% Align the data
marker_F = squeeze(data(:,:,4));
marker_M = squeeze(data(:,:,5));
rotangle = atan2(-(marker_F(:,2)-marker_M(:,2)),(marker_F(:,1)-marker_M(:,1)));
global_rotmatrix = zeros(2,2,numel(rotangle));
global_rotmatrix(1,1,:) = cos(rotangle);
global_rotmatrix(1,2,:) = -sin(rotangle);
global_rotmatrix(2,1,:) = sin(rotangle);
global_rotmatrix(2,2,:) = cos(rotangle);

for nMarker = 1:size(data,3)
    data(:,:, nMarker) = squeeze(data(:,:,nMarker)) - marker_M;
end

for nMarker = 1:size(data,3)
    for nFrame = 1:size(global_rotmatrix,3)
        marker = squeeze(data(nFrame, 1:2, nMarker));
        data(nFrame, 1:2, nMarker) = squeeze(global_rotmatrix(:,:,nFrame))*marker';
    end
end

% %% Look at the data to check it is aligned. 
% close all
% a = Keypoint3DAnimator(data, skeleton);
%% PCA
[coeff, score, latent, tsquared, explained] = pca(reshape(data, size(data, 1), []));

%% Build a stacked trace viewer with keypoints
% Move with the LR arrow keys. increase/decrease speed with Up/Down. 
nPCs = 5;
close all;
h = cell(2,1);
h{1} = StackedTraceAnimator((1:size(score,1))', score(:,1:nPCs));
h{2} = Keypoint3DAnimator(data, skeleton);
set(h{2}.Axes,'XLim',[-140 140],'YLim',[-75 75],'ZLim',[-120 120],...
    'CameraPosition',[985.0624 744.2830 332.1425])
Animator.tileAnimators(h)
Animator.linkAll(h)