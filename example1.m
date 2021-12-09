%% Example 1: Run analyzePRF on an example dataset

%% Download dataset (if necessary) and add analyzePRF to the MATLAB path

setup;

%% Load in the data

% Load in the data
load('exampledataset.mat');

% Check the workspace
whos
%%

%% Inspect the data

% Check dimensionality of the data
data
%%

% There are four runs of data; each run consists of 150 time points (TR = 2 s).
% The data have already been pre-processed for slice timing correction, motion 
% correction, and spatial undistortion.  For simplicity, we have selected
% 10 example voxels from the left hemisphere.  Let's visualize the time-series 
% data for the second voxel.
temp = cellfun(@(x) x(2,:),data,'UniformOutput',0);
figure; hold on;
set(gcf,'Units','points','Position',[100 100 600 150]);
plot(cat(2,temp{:}),'r-');
straightline(150*(1:4)+.5,'v','g-');
xlabel('TR');
ylabel('BOLD signal');
ax = axis;
axis([.5 600+.5 ax(3:4)]);
title('Time-series data');
%%

%% Inspect the stimuli

% Check dimensionality of the stimuli
stimulus
%%

% The stimulus images have been prepared at a resolution of 100 pixels x 100 pixels.
% There are 300 images per run because we have prepared the images at a time resolution
% of 1 second.  (Note that this is faster than the data sampling rate.  When analyzing
% the data, we will resample the data to match the stimulus rate.)  Let's inspect a 
% few of the stimulus images in the first run.
figure;
set(gcf,'Units','points','Position',[100 100 700 300]);
for p=1:3
  subplot(1,3,p); hold on;
  num = 239+2*p;
  imagesc(stimulus{1}(:,:,num),[0 1]);
  axis image tight;
  set(gca,'YDir','reverse');
  colormap(gray);
  title(sprintf('Image number %d',num));
end
%%

% Notice that the range of values is 0 to 1 (0 indicates that the gray background was
% present; 1 indicates that the stimulus was present).  For these stimulus images,
% the stimulus is a bar that moves downward and to the left.

%% Analyze the data

% Start parallel MATLAB to speed up execution.
if isempty(gcp)
  parpool;
end

% We need to resample the data to match the temporal rate of the stimulus.  Here we use 
% cubic interpolation to increase the rate of the data from 2 seconds to 1 second (note 
% that the first time point is preserved and the total data duration stays the same).
data = tseriesinterp(data,2,1,2);

% Finally, we analyze the data using analyzePRF.  The third argument is the TR, which 
% is now 1 second.  Here, we set 'seedmode' to [0 1] which means to just use 
% two generic initial seeds.  We suppress command-window output by 
% setting 'display' to 'off'.
%
% Note that you may want to try using 'seedmode' set to 2, which uses 
% a single initial seed for each given voxel based on a very large
% "super grid" collection of possible parameter combinations. Also, you may want to
% try using 'seedmode' set to -2 which skips costly nonlinear optimization and
% just returns the best of the "supergrid" found for each voxel (this is useful
% for getting a pretty good answer in very little time).
results = analyzePRF(stimulus,data,1,struct('seedmode',[0 1],'display','off'));
%%

% Note that because of the use of parfor, the command-window output for different
% voxels will come in at different times (and so the text above is not really
% readable).

%% Inspect the results

% The stimulus is 100 pixels (in both height and weight), and this corresponds to
% 10 degrees of visual angle.  To convert from pixels to degreees, we multiply
% by 10/100.
cfactor = 10/100;

% Visualize the location of each voxel's pRF
figure; hold on;
set(gcf,'Units','points','Position',[100 100 400 400]);
cmap = jet(size(results.ang,1));
for p=1:size(results.ang,1)
  xpos = results.ecc(p) * cos(results.ang(p)/180*pi) * cfactor;
  ypos = results.ecc(p) * sin(results.ang(p)/180*pi) * cfactor;
  ang = results.ang(p)/180*pi;
  sd = results.rfsize(p) * cfactor;
  h = drawellipse(xpos,ypos,ang,2*sd,2*sd);  % circle at +/- 2 pRF sizes
  set(h,'Color',cmap(p,:),'LineWidth',2);
  set(scatter(xpos,ypos,'r.'),'CData',cmap(p,:));
end
drawrectangle(0,0,10,10,'k-');  % square indicating stimulus extent
axis([-10 10 -10 10]);
straightline(0,'h','k-');       % line indicating horizontal meridian
straightline(0,'v','k-');       % line indicating vertical meridian
axis square;
set(gca,'XTick',-10:2:10,'YTick',-10:2:10);
xlabel('X-position (deg)');
ylabel('Y-position (deg)');
%%

% Please see the example2.m script for an example of how to inspect the model fit 
% and compare it to the data.
