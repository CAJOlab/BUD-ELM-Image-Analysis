%% Initialize
p = 995.67;  %Fluid Density (kg/m^3)
u = 0.000797;   % Dynamic Viscosity (Pa*s) for water at 30 C
d0 = 0.05;  % Orbital size (m)
T = 30; % Temperature (C)
v = 8.005*10^-7; % Kinematic Viscosity (m^2/s) for water at 30 C

imagefiles = dir('*.tif');      %Parse in all .tif files in folder
Conv = readtable('Flask Diam Conversions mm_px.csv');  % List of mm/px conversions for each sample
nfiles = length(imagefiles);    % Number of files found
PCAArray = zeros(nfiles,8);  % Initialize data array
areaslist = {}; % Initialize array for each set of areas

%% Process files 

for ii = 1:nfiles  % for all files
   currentfilename = imagefiles(ii).name; % Get current file
   [d,n,V] = ConvertFlask(currentfilename); % Read in flask parameters from file name
   conv_factor = px2cm(Conv, currentfilename); % Read in conversion factor table
   mask = imread(currentfilename);  % Read in segmentation file
   mask = logical(mask - 1); % Convert Segmentation file to mask
   
   output_props = regionprops(mask, 'Area'); %Get list of all areas from mask
   areas = zeros(size(output_props)); %initialize areas matrix
   for i = 1:length(areas)  %for each area
       areas(i) = output_props(i).Area;  % add current area to areas matrix
   end
   sorted_areas = sort(areas); % Sort areas smallest to largest
   
   percindex = round(0.95*length(sorted_areas));   %Top 5 Percentile
   areaslist{1,ii} = sorted_areas(percindex:end);  % Get all areas that are within the top 5 Percentile
   Top5Perc = mean(sorted_areas(percindex:end)); % Average areas in top 5 percentile
   StDTop5Perc= std(sorted_areas(percindex:end))/sqrt(length(sorted_areas(percindex:end))); % Standard error of said Areas
   
   PCAArray(ii,1) = d;  % Flask Diameter
   PCAArray(ii,2) = n; % Shaking Speed
   PCAArray(ii,3) = V; % Culture Volume
   PCAArray(ii,4) = GetVolumetricPower(p, n, d, u, V); % Volumetric Power Input
   PCAArray(ii,5) = GetkLa(d, n, d0, V, v, T, u); % kLa
   PCAArray(ii,6) = (d^5).*PCAArray(ii,5).*PCAArray(ii,4)*1000; % Modified Volumetric Power
   PCAArray(ii,7) = Top5Perc*conv_factor; % Size of materials
   PCAArray(ii,8) = PCAArray(ii,5).*PCAArray(ii,4)*1000; % PV * kLa
end

%% Modified Volumetric Power
UniquePVA = sort(unique(PCAArray(:,6)));  % Sort all unique Modified Volumetric Power values
AvgPerc5 = zeros(1,length(UniquePVA)); % Initialize Average Array
StDPerc5 = AvgPerc5; % Initialize Standard Error Array
for ii = 1:length(UniquePVA) % for number of Unique PVA values
    PVindex = find(PCAArray(:,6) == UniquePVA(ii));  %Find all flasks with that PVA Value
    
    AvgPerc5(ii) = mean(PCAArray(PVindex,7));  % Average the Material Size from those flasks
    StDPerc5(ii) = std(PCAArray(PVindex,7))/sqrt(length(PCAArray(PVindex,7))); % Calculate Standard Error of Average
end
errorbar(UniquePVA, AvgPerc5, StDPerc5) % Plot

%% PV
UniquePV = sort(unique(PCAArray(:,4))); % Sort all unique Volumetric Power Input values
AvgPerc5 = zeros(1,length(UniquePVA)); % Initialize Average Array
StDPerc5 = AvgPerc5; % Initialize Standard Error Array
for ii = 1:length(UniquePV) % for number of Unique PV values
    PVindex = find(PCAArray(:,4) == UniquePV(ii)); %Find all flasks with that PV Value
    
    AvgPerc5(ii) = mean(PCAArray(PVindex,7)); % Average the Material Size from those flasks
    StDPerc5(ii) = std(PCAArray(PVindex,7))/sqrt(length(PCAArray(PVindex,7)));  % Calculate Standard Error of Average
end
errorbar(UniquePV, AvgPerc5, StDPerc5) %Plot 

%% kLa
UniquekLa = sort(unique(PCAArray(:,5))); % Sort all unique kLa values
AvgPerc5 = zeros(1,length(UniquePVA)); % Initialize Average Array
StDPerc5 = AvgPerc5; % Initialize Standard Error Array
for ii = 1:length(UniquekLa) % for number of Unique kLa values
    kLaindex = find(PCAArray(:,5) == UniquekLa(ii)); %Find all flasks with that kLa Value
    
    AvgPerc5(ii) = mean(PCAArray(kLaindex,7)); % Average the Material Size from those flasks
    StDPerc5(ii) = std(PCAArray(kLaindex,7))/sqrt(length(PCAArray(kLaindex,7)));  % Calculate Standard Error of Average
end
errorbar(UniquekLa, AvgPerc5, StDPerc5) % Plot

%% PV * kLa
UniquePVkLa = sort(unique(PCAArray(:,8))); % Sort all unique PV*kLa values

%% Functions
function [d,n,V] = ConvertFlask(currentfilename) %Gets Flask Parameters from File name
   flask = double(string(currentfilename(1:3)));
   rpm = double(string(currentfilename(5:7)));
   volume = double(string(currentfilename(9:11)));
   
   d = FlaskSize(flask);
   n = ShakingSpeed(rpm);
   V = GetVolume(volume);
end

function f = FlaskSize(flask)  % Get flask diameter in m from flask size
    if flask == 125
        f = 0.06011;
    elseif flask == 250
        f = 0.0744;
    elseif flask == 500
        f = 0.09553;
    end
end

function s = ShakingSpeed(rpm) % Convert Shaking Speed from rpm to s^-1
    conversion = 0.01666666666;
    s = conversion*rpm;
end

function Vol = GetVolume(volume)  % Convert Volume from mL to m^3
    Vol = 0.000001*volume;
end

function factor = px2cm(Conversion_table, file)   % Calculate pixel to cm^2 conversion factor for a particular sample
    [index,~] = find(Conversion_table{:,1}==string(file));
    convert = 0.1*Conversion_table{index,2};
    factor = convert(1)^2;
end
function Re = GetReynolds(p, n, d, u) % Calculate Reynolds Number
    Re = p*n*(d^2)/u;
end

function Ne = GetNewtons(p, n, d, u) % Calculate Newtons Number
    Re = GetReynolds(p, n, d, u);
    Ne = 70*(Re^-1) + 25*(Re^-0.6) + 1.5*(Re^-0.2);
end

function PV = GetVolumetricPower(p, n, d, u, V) % Calculate Volumetric Power, in W/(m^3)
    Ne = GetNewtons(p, n, d, u);
    PV = Ne*p*(n^3)*(d^4)/(V^(2/3));
end

function DO2 = GetDiffuseConstant(T, u)  % Calculate estimate for Diffusion Constant
    TK = T + 273.15;
    Y = 2.26;
    MH20 = 18;
    uCP = 1000*u;
    VO2 = 25.6;
    K = 7.4*10^-8;
    DO2 = (K*TK*((Y*MH20)^0.5))/(uCP*(VO2^0.6));
end

function kLa = GetkLa(d, n, d0, VL, v, T, u) % Calculate kLa in s^-1
    D = GetDiffuseConstant(T, u);
    g = 9.807;
    kLa = 0.5 * (d^(73/36)) * n * (d0^(1/4)) * (VL^(-8/9)) * (D^(1/2)) * (v^(-13/54)) * (g^(-7/54));
end   

