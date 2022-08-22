imagefiles = dir('*.png');      % Get All raw files
nfiles = length(imagefiles);    % Number of files found
for ii=1:nfiles       % For all raw files in folder
   currentfilename = imagefiles(ii).name;  % Get current file name
   currentimage = imread(currentfilename);  % read file
   images{ii} = currentimage; 
   a = uint8(zeros(4000, 6000)); % initialize blue channel variable
  
   a(:, :) = currentimage(:, :, 3) + 50; % acquire blue channel and brighten image

   Jc = imadjust(a, [0.35, 0.60]); % Improve image contrast for ilastik

   imwrite(Jc, strcat('C:\Users\rft2\Pictures\500 mL Flasks Circles\', currentfilename, '_B.tif')); % Save Blue Channel
end 