%% Initialize
p = 995.67;  %Fluid Density (kg/m^3)
u = 0.000797;   % Dynamic Viscosity (Pa*s)  0.000797 for water at 30
d0 = 0.05;  % Orbital size (mm)
T = 30; % Temperature (C)
v = 8.005*10^-7; % constant for ???

imagefiles = dir('*.tif');      %Parse in all .tif files in folder
Conv = readtable('Flask Diam Conversions mm_px.csv');  % List of mm/px conversions for each sample
nfiles = length(imagefiles);    % Number of files found
PCAArray = zeros(nfiles,11);  % Initialize data array
areaslist = {}; % Initialize array for each set of areas
TopPercs75to99 = zeros(nfiles, 6);
%% 

for ii = 1:nfiles  %
   currentfilename = imagefiles(ii).name;
   [d,n,V] = ConvertFlask(currentfilename);
   conv_factor = px2cm(Conv, currentfilename);
   mask = imread(currentfilename);
   mask = logical(mask - 1);
   
   output_props = regionprops(mask, 'Area');
   areas = zeros(size(output_props));
   for i = 1:length(areas)
       areas(i) = output_props(i).Area;
   end
   sorted_areas = sort(areas);
   
   percindex = round(0.95*length(sorted_areas));   %Top 5 Percentile
   areaslist{1,ii} = sorted_areas(percindex:end);  
   Top5Perc = mean(sorted_areas(percindex:end));
   StDTop5Perc= std(sorted_areas(percindex:end))/sqrt(length(sorted_areas(percindex:end)));
   
   
  % [TestTop, TestStD] = TopPercentileAreas(0.95, sorted_areas);
  
   PCAArray(ii,1) = d;
   PCAArray(ii,2) = n;
   PCAArray(ii,3) = V;
   PCAArray(ii,4) = GetVolumetricPower(p, n, d, u, V);
   PCAArray(ii,5) = GetkLa(d, n, d0, V, v, T, u);
   PCAArray(ii,6) = (d^5).*PCAArray(ii,5).*PCAArray(ii,4)*1000;
   PCAArray(ii,7) = Top5Perc*conv_factor;
   PCAArray(ii,8) = StDTop5Perc*conv_factor;
   PCAArray(ii,9) = PCAArray(ii,5).*PCAArray(ii,4)*1000;
   
   
end

UniquePVkLa = sort(unique(PCAArray(:,6)));
PVkLaAvgAreas = zeros(1,length(UniquePVkLa));
PVkLaNumSamples = zeros(1,length(UniquePVkLa));
PVkLaStD = zeros(1,length(UniquePVkLa));
AvgPerc5 = zeros(1,length(UniquePVkLa));
StDPerc5 = AvgPerc5;
lenNum = zeros(1,20);
for ii = 1:length(UniquePVkLa)
    PVindex = find(PCAArray(:,6) == UniquePVkLa(ii));
    
    tempAreas =[];
    for kk = 1:length(PVindex)
        index = PVindex(kk);
        tempAreas = [tempAreas;(cell2mat(areaslist(index)))];
    end
    %PVkLaAvgAreas(ii) = mean(tempAreas)*0.00001032;
    PVkLaNumSamples(ii) = length(tempAreas);
    %PVkLaStD(ii) = std(tempAreas)/PVNumSamples(ii)*0.00001032;
    AvgPerc5(ii) = mean(PCAArray(PVindex,7));
    StDPerc5(ii) = std(PCAArray(PVindex,7))/sqrt(length(PCAArray(PVindex,8)));
    lenNum(1,ii) = length(PCAArray(PVindex,8));
end
errorbar(UniquePVkLa, AvgPerc5, StDPerc5)

%% PV
UniquePV = sort(unique(PCAArray(:,4)));
for ii = 1:length(UniquePV)
    PVindex = find(PCAArray(:,4) == UniquePV(ii));
    
    tempAreas =[];
    for kk = 1:length(PVindex)
        index = PVindex(kk);
        tempAreas = [tempAreas;(cell2mat(areaslist(index)))];
    end
    PVAvgAreas(ii) = mean(tempAreas)*0.00001032;
    PVNumSamples(ii) = length(tempAreas);
    PVStD(ii) = std(tempAreas)/PVNumSamples(ii)*0.00001032;
    AvgPerc5(ii) = mean(PCAArray(PVindex,8));
    StDPerc5(ii) = std(PCAArray(PVindex,8))/length(PCAArray(PVindex,8));
end
errorbar(UniquePV, AvgPerc5, StDPerc5)
%% kLa
for ii = 1:length(UniquePVkLa)
    kLaindex = find(PCAArray(:,5) == UniquekLa(ii));
    
    tempAreas =[];
    for kk = 1:length(kLaindex)
        index = kLaindex(kk);
        tempAreas = [tempAreas;(cell2mat(areaslist(index)))];
    end
    PVAvgAreas(ii) = mean(tempAreas)*0.00001032;
    PVNumSamples(ii) = length(tempAreas);
    PVStD(ii) = std(tempAreas)/PVNumSamples(ii)*0.00001032;
    AvgPerc5(ii) = mean(PCAArray(PVindex,8));
    StDPerc5(ii) = std(PCAArray(PVindex,8))/length(PCAArray(PVindex,8));
end
errorbar(UniquekLa, AvgPerc5, StDPerc5)
%% Functions

function [d,n,V] = ConvertFlask(currentfilename)
   flask = double(string(currentfilename(1:3)));
   rpm = double(string(currentfilename(5:7)));
   volume = double(string(currentfilename(9:11)));
   
   d = FlaskSize(flask);
   n = ShakingSpeed(rpm);
   V = GetVolume(volume);
end

function f = FlaskSize(flask)
    if flask == 125
        f = 0.06011;
    elseif flask == 250
        f = 0.0744;
    elseif flask == 500
        f = 0.09553;
    end
end

function s = ShakingSpeed(rpm)
    conversion = 0.01666666666;
    s = conversion*rpm;
end

function Vol = GetVolume(volume)
    Vol = 0.000001*volume;
end

function factor = px2cm(Conversion_table, file)
    [index,~] = find(Conversion_table{:,1}==string(file));
    convert = 0.1*Conversion_table{index,2};
    factor = convert(1)^2;
end
function Re = GetReynolds(p, n, d, u)
    Re = p*n*(d^2)/u;
end

function Ne = GetNewtons(p, n, d, u)
    Re = GetReynolds(p, n, d, u);
    Ne = 70*(Re^-1) + 25*(Re^-0.6) + 1.5*(Re^-0.2);
end

function PV = GetVolumetricPower(p, n, d, u, V)
    Ne = GetNewtons(p, n, d, u);
    PV = Ne*p*(n^3)*(d^4)/(V^(2/3));
end

function DO2 = GetDiffuseConstant(T, u)
    TK = T + 273.15;
    Y = 2.26;
    MH20 = 18;
    uCP = 1000*u;
    VO2 = 25.6;
    K = 7.4*10^-8;
    DO2 = (K*TK*((Y*MH20)^0.5))/(uCP*(VO2^0.6));
end

function kLa = GetkLa(d, n, d0, VL, v, T, u)
    D = GetDiffuseConstant(T, u);
    g = 9.807;
    kLa = 0.5 * (d^(73/36)) * n * (d0^(1/4)) * (VL^(-8/9)) * (D^(1/2)) * (v^(-13/54)) * (g^(-7/54));
end   

function [TopPerc, StDTopPerc] = TopPercentileAreas(percentile, sorted_areas_list) 
   percindex = round(percentile*length(sorted_areas_list));   %Top 5 Percentile
   areaslist = sorted_areas_list(percindex:end);  
   TopPerc = mean(areaslist);
   StDTopPerc= std(areaslist)/sqrt(length(areaslist));
end