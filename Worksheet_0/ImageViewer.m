function ImageViewer(filename, dimensions, location)
% function ImageViewer(filename, dimensions, location)
% 
% This is a simple matlab function to visualise raw byte filethat stores 3-D
% image volumes.
% 
% INPUTS:
%     1. filename: the raw byte file
%     2. dimensions: the dimensions of the image volume specified as a 
%     3-by-1 array.
%     3. spacings: the voxel specings specified as a 3-by-1 array.
%     4. location: the coordinates of the point where the axial, coronal and
%     sagittal pLanes intersect
% OUTPUTS:

%% 1. Load the data from file

% open the raw byte file
fileID = fopen(filename,'r','l');

% load the data from the file
data = fread(fileID);

% close the file
fclose(fileID);

%% 2. Convert the data, 1-D array into a 3-D array
volume = reshape(data,dimensions);

% Clear the data
clear data

%% 3. Extract the three 2-D slices
axial = volume(:,:, location(3));

coronal = volume(:,location(2),:);

sagittal = volume(location(1),:,:);

% Remove the singleton dimension
axial = squeeze(axial);
coronal = squeeze(coronal);
sagittal = squeeze(sagittal);

% Plot the slices
imagescxy(permute(axial,[2 1]));

axis equal;

axis xy;

title('axial');

imagescxy(permute(sagittal,[2,1]));

axis equal;

axis xy;

title('sagittal');

imagescxy(permute(coronal,[2 1]));

axis equal;

axis xy;

title('coronal');

end

