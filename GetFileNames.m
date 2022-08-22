imagefiles = dir('*.tif');      
nfiles = length(imagefiles);    % Number of files found
namelist = {};
for ii=1:nfiles
   currentfilename = imagefiles(ii).name;
   disp(currentfilename);
   namelist{end+1} = currentfilename;
   x = imread(currentfilename);
   x = logical(x-1);
   x = bwareaopen(x, 10);
   imshow(x);
   pause;
end 