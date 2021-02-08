
close all;
clc;
clear;

matrix_size= 256;
% this is a standard shepp-logan phantom often used in MRI 
addpath 'Jenny_Code'
load('data.mat')
P = brain';

subplot(1,2,1); 
imagesc(P); 
axis('square'); 
% Flip the y-axis
set(gca,'YDir','normal')
colormap('gray'); 
title('Brain Image')

pause(1);
% push to k-space (by doing FT in y and x directions)
k_P = itok(P);
subplot(1,2,2); 
imagesc(abs(k_P), [0 500]); 
axis('square'); 
% Flip the y-axis
set(gca,'YDir','normal')
% Apply the custom colormap.
colormap('gray'); 
% Set a log scale
set(gca,'ColorScale','log')
title('k-space: Brain Image')
% This part ^^^^^^^   just scaled the image so that we can
% see the contrast better

% add a breakpoint here

%% -----------------------------------------------------------------------
% % what information does the centre of k-space contain?
% 
k_P_inside = zeros(matrix_size,matrix_size);
% just take 40 lines from the middle of k-space
k_P_inside((matrix_size/2)-20:(matrix_size/2)+20, (matrix_size/2)-20:(matrix_size/2)+20) = k_P((matrix_size/2)-20:(matrix_size/2)+20, (matrix_size/2)-20:(matrix_size/2)+20);

subplot(1,3,1); 
imagesc(abs(k_P_inside), [0 500]); 
axis('square');
colormap('gray');
set(gca,'ColorScale','log')
title('Middle of k-space: Brain Image');
hold on

% now push back to image space
i_P_inside = ktoi(k_P_inside);

subplot(1,3,2); 
imagesc(abs(i_P_inside)); 
axis('square'); 
colormap('gray'); 
title('Image from middle of k-space: Brain Image')
% it is a low resolution image (the rings are caused by truncartion artefacts)

% add a breakpoint here

% you could remove these by using a guassian type filter (here i have used a 'hanning' filter as you can control the parameters)
transition_width = 8;
flat_radius = 12;
[filter] = hanning2d(transition_width, flat_radius, matrix_size);
k_P_inside_filtered = filter .* k_P;
% now push back to image space
i_P_inside_filtered = ktoi(k_P_inside_filtered);

subplot(1,3,3); 
imagesc(abs(i_P_inside_filtered)); 
axis('square'); 
colormap('gray'); 
title('Image from FILTERED middle of k-space: Brain Image')
% the rings have been reduced, if not removed(!) but the resolution is also
% further reduced (as i still only kept the middle 40 lines, but only 24 - flat_radius*2 - were fully maintained)

% add a breakpoint here

%% -----------------------------------------------------------------------
% what information does the outside of k-sapce contain?

k_P_outside = k_P;
% zero all but the 40 lines around the outside of k-space
k_P_outside(40:matrix_size-40, 40:matrix_size-40) = 0;
subplot(1,2,1); 
imagesc(abs(k_P_outside), [0 200]); 
axis('square'); 
colormap('gray'); 
title('Outer part of k-space: Brain Image');

% now push back to image space
i_P_outside = ktoi(k_P_outside);
subplot(1,2,2); 
imagesc(abs(i_P_outside)); 
axis('square'); colormap('gray'); 
title('Image from middle of k-space: Brain Image')
% this image contains the high frequency parts - the edges

% add a breakpoint here

%% -----------------------------------------------------------------------
% THE RELATIONSHIP

% pixel_size = 1/k_FOV

%       field of view in k-space: k_FOV

% so by reducing the number of lines acquired in k-space (the matrix size) we reduce
% the spatial resolution in image space

% see: http://mriquestions.com/field-of-view-fov.html

% IN ML: in super-resolution type studies we acquire low
% resoltuion data and try to recover high resolution images

%% -----------------------------------------------------------------------
% what if we only acquired alternate lines in k-space?

k_P_x2 = zeros(matrix_size,matrix_size);
% collect only odd lines of k-space (the other direction remains fully sampled)
k_P_x2(:,1:2:end) = k_P(:,1:2:end);
figure; imagesc(abs(k_P_x2), [0 200]); axis('square'); colormap('jet'); title('Only odd lines of k-space: Shepp-Logan phantom');


% now push back to image space
i_P_x2 = ktoi(k_P_x2);
figure; 
imagesc(abs(i_P_x2)); 
axis('square'); 
colormap('gray');
title('Image from odd lines of k-space: Brain Image')
% this image contains aliasing (or wrap) in the direction which we missed out the alternative lines
% in-fact because this is Cartesian data, the image just folds over onto
% itself

% add a breakpoint here

% if we did the undersampling in the other direction
k_P_x2b = zeros(matrix_size,matrix_size);
% collect only odd lines of k-space (the other direction remains fully sampled)
k_P_x2b(1:2:end,:) = k_P(1:2:end, :);
figure; 
imagesc(abs(k_P_x2b), [0 200]); 
axis('square'); 
colormap('jet'); 
title('Only odd lines of k-space: Brain Image');

pause(1);

% now push back to image space
i_P_x2b = ktoi(k_P_x2b);
figure; 
imagesc(abs(i_P_x2b)); 
axis('square'); 
colormap('gray'); 
title('Image from odd lines of k-space: Brain Image')
% this image contains aliasing (or wrap) in the other direction 

%%
% add a breakpoint here

%if we only acauired one in every three lines what do you think would
%happen?
k_P_x3 = zeros(matrix_size,matrix_size);
% collect only odd lines of k-space (the other direction remains fully sampled)
k_P_x3(:,1:3:end) = k_P(:,1:3:end);
figure; imagesc(abs(k_P_x3), [0 200]); axis('square'); colormap('jet'); title('Only one in every three lines of k-space: Shepp-Logan phantom');

pause(1);

% now push back to image space
i_P_x3 = ktoi(k_P_x3);
figure; imagesc(abs(i_P_x3)); axis('square'); colormap('gray'); title('Image from one in every three lines of k-space: Shepp-Logan phantom')
% this image contains aliasing (or wrap) in the direction which we missed out the alternative lines
% this time there are 3 replicas of the image wrapped onto each other

% add a breakpoint here

%% -----------------------------------------------------------------------
% THE RELATIONSHIP

% distance_k = 1/FOV

%       distance between the lines in k-space: distance_k
%       field of view of the image: FOV

% so by increasing the distancde between the lines in k-space (but keeping 
% the k-space FOV the same, i.e. still going out the same distance in k-space) 
% we are effectively reducing the imaging FOV

% see: http://mriquestions.com/field-of-view-fov.html

% IN ML: in de-alisaing type studies we acquire undersampled k-data (so the 
% images contain artefacts) and try to recover artefact-free  images

% -----------------------------------------------------------------------
% so what about other trajectories?
% See my thesis and movie, but briefly in Cartesian imaging only 1 lines of k-space
% is acquired afte each RF excitation, so the time acquiring data is short
% compared to the repetition time. Other trajectories may enable more
% efficient fillig of k-space where large proportions of k-space are
% acquired in each repetition time.This includes Echo Planar Imaging (EPI)
% - which is essentially Cartesian - and spirals
% 
% See: http://mriquestions.com/k-space-trajectories.html
% 
% SPIRALS
% this is non-Cartesian trajectory which means we cannot use the FFT to
% reconstruct the data. So we use 'gridding' to put the non-Cartesian data
% back onto a Cartesian grid (see a part of my theses and 
% http://mriquestions.com/uploads/3/4/5/7/34572113/pauly-non-cartesian_recon.pdf
% ) . 

% -----------------------------------------------------------------------
% here is a spiral example

Nints    = 24;               % number of spiral interleaves required to fully sample k-space
acc_fact = 1;
[trajectory, weights] = CalculateSpiralTrajectoryDL(1, Nints, acc_fact);
trajectory = trajectory(:,1:2);
        
k_P_spiral   = grid_data_bck(trajectory, k_P, [0 0]); 

% now push spiral data back to Cartesian k-space
k_P_spiralGridded = grid_data(trajectory, k_P_spiral, weights, [matrix_size, matrix_size], [0 0]);
figure; imagesc(abs(k_P_spiralGridded), [0 200]); axis('square'); colormap('jet'); title('Spiral trajectory in k-space: Shepp-Logan phantom');

pause(1);

% push gridded datato image space
i_P_spiral_gridded = ktoi(k_P_spiralGridded);
figure; imagesc(abs(i_P_spiral_gridded)); axis('square'); colormap('gray'); title('Image from spiral lines of k-space: Shepp-Logan phantom')
% so you can recover the same data with the same resolution using a spiral
% trajectory

% add a breakpoint here

%% -----------------------------------------------------------------------
% what if you undersample a spiral trajectory?

% just acquire odd spokes
acc_fact = 2;
[trajectoryx2, weightsx2] = CalculateSpiralTrajectoryDL(1, Nints, acc_fact);
trajectoryx2 = trajectoryx2(:,1:2);
        
k_P_spiralx2   = grid_data_bck(trajectoryx2, k_P, [0 0]); 

% now push spiral data back to Cartesian k-space
k_P_spiralGridded_x2 = grid_data(trajectoryx2, k_P_spiralx2, weightsx2, [matrix_size, matrix_size], [0 0]);
figure; imagesc(abs(k_P_spiralGridded_x2), [0 200]); axis('square'); colormap('jet'); title('undersampeld x2 Spiral trajectory in k-space: Shepp-Logan phantom');

pause(1);

% push gridded datato image space
i_P_spiral_gridded_x2 = ktoi(k_P_spiralGridded_x2);
figure; imagesc(abs(i_P_spiral_gridded_x2)); axis('square'); colormap('gray'); title('Image from odd spiral lines of k-space: Shepp-Logan phantom')
% the artefact pattern is much more complex

% add a breakpoint here

%% -----------------------------------------------------------------------
% just acquire 1 in every 4 spokes
% just acquire odd spokes
acc_fact = 4;
[trajectoryx4, weightsx4] = CalculateSpiralTrajectoryDL(1, Nints, acc_fact);
trajectoryx4 = trajectoryx4(:,1:2);
        
k_P_spiralx4   = grid_data_bck(trajectoryx4, k_P, [0 0]); 

% now push spiral data back to Cartesian k-space
k_P_spiralGridded_x4 = grid_data(trajectoryx4, k_P_spiralx4, weightsx4, [matrix_size, matrix_size], [0 0]);
figure; imagesc(abs(k_P_spiralGridded_x4), [0 200]); axis('square'); colormap('jet'); title('undersampled x4 Spiral trajectory in k-space: Shepp-Logan phantom');

pause(1);

% push gridded datato image space
i_P_spiral_gridded_x4 = ktoi(k_P_spiralGridded_x4);
figure; imagesc(abs(i_P_spiral_gridded_x4)); axis('square'); colormap('gray'); title('Image from every fourth spiral lines of k-space: Shepp-Logan phantom')
% the artefact pattern is much more complex

% add a breakpoint here

%% -----------------------------------------------------------------------
%% -----------------------------------------------------------------------
% here is a Radial example

[trajectory, weights] = CalculateRadialTrajectoryDL_sortedGA(matrix_size, 1, 1);
% you do not need to look into the 'CalculateRadialTrajectoryDL_sortedGA' code - just accept it as it is for
% now
trajectory = trajectory(:,1:2);

% now we have to sample the image on this trajectory
% we do back gridding, from the Cartesian data onto the radial trajectory
        
k_P_radial   = grid_data_bck(trajectory, k_P, [0 0]); 

% now push spiral data back to Cartesian k-space
k_P_radialGridded = grid_data(trajectory, k_P_radial, weights, [matrix_size, matrix_size], [0 0]);
subplot(1,2,1); 
imagesc(abs(k_P_radialGridded), [0 2000]); 
axis('square'); colormap('jet'); 
title('US Radial trajectory in k-space: Brain Image');


% push gridded datato image space
i_P_radial_gridded = ktoi(k_P_radialGridded);
subplot(1,2,2); 
imagesc(abs(i_P_radial_gridded)); 
axis('square'); 
colormap('gray'); 
title('Image from US radial lines of k-space: Brain Image')
% so you can recover the same data with the same resolution using a radial
% trajectory

% add a breakpoint here

%% -----------------------------------------------------------------------
% what if you undersample a radial trajectory?

% just acquire odd spokes
[trajectory_x2, weights_x2] = CalculateRadialTrajectoryDL_sortedGA(matrix_size, 1, 2);
% you do not need to look into the 'CalculateRadialTrajectoryDL_sortedGA' code - just accept it as it is for
% now
trajectory_x2 = trajectory_x2(:,1:2);

% now we have to sample the image on this trajectory
% we do back gridding, from the Cartesian data onto the radial trajectory
        
k_P_radial_x2   = grid_data_bck(trajectory_x2, k_P, [0 0]); 

% now push spiral data back to Cartesian k-space
k_P_radialGridded_x2 = grid_data(trajectory_x2, k_P_radial_x2, weights_x2, [matrix_size, matrix_size], [0 0]);
subplot(1,2,1); 
imagesc(abs(k_P_radialGridded_x10), [0 2000]); 
axis('square'); colormap('jet'); 
title('US x2 Radial trajectory in k-space: Brain Image');

% push gridded datato image space
i_P_radial_gridded_x2 = ktoi(k_P_radialGridded_x2);
subplot(1,2,2); 
imagesc(abs(i_P_radial_gridded_x2)); 
axis('square'); 
colormap('gray'); 
title('Image from x2 US radial lines of k-space: Brain Image')

% the artefact pattern here is much more subtle - in fact radials are
% relatively robust to underampleing because they are oversampled in the
% centre of k-space

% add a breakpoint here

%% -----------------------------------------------------------------------
% just acquire 1 in every 10 spokes

[trajectory_x10, weights_x10] = CalculateRadialTrajectoryDL_sortedGA(matrix_size, 1, 10);
% you do not need to look into the 'CalculateRadialTrajectoryDL_sortedGA' code - just accept it as it is for
% now
trajectory_x10 = trajectory_x10(:,1:2);

% now we have to sample the image on this trajectory
% we do back gridding, from the Cartesian data onto the radial trajectory
        
k_P_radial_x10   = grid_data_bck(trajectory_x10, k_P, [0 0]); 

% now push spiral data back to Cartesian k-space
k_P_radialGridded_x10 = grid_data(trajectory_x10, k_P_radial_x10, weights_x10, [matrix_size, matrix_size], [0 0]);
subplot(1,2,1); 
imagesc(abs(k_P_radialGridded_x10), [0 2000]); 
axis('square'); colormap('jet'); 
title('US x10 Radial trajectory in k-space: Brain Image');


% push gridded datato image space
i_P_radial_gridded_x10 = ktoi(k_P_radialGridded_x10);
subplot(1,2,2); 
imagesc(abs(i_P_radial_gridded_x10)); 
axis('square'); 
colormap('gray'); 
title('Image from x10 US radial lines of k-space: Brain Image')
% the artefact pattern here is much more subtle - in fact radials are
% relatively robust to underampleing because they are oversampled in the
% centre of k-space