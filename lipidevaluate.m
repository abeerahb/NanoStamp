function [outputMatrix,TumourSp,GreySp,WhiteSp,VolcanoPlots]=lipidevaluate(inputfile,p_val_thres,fold_thres,skipFiltering)
% Example usage: lipidevaluate('feature_matrix_brain_B3_T1_Day5.xlsx',0.05,2,true)
if nargin < 1
    error('Input file name is required.');
end
if nargin < 2 || isempty(p_val_thres), p_val_thres = 0.05; end
if nargin < 3 || isempty(fold_thres), fold_thres = 2; end
if nargin < 4 || isempty(skipFiltering), skipFiltering = true; end

% e.g.,
% [B3_D0T1_TopMax_outputMatrix,B3_D0T1_TopMax_TumourSp,B3_D0T1_TopMax_GreySp,B3_D0T1_TopMax_WhiteSp,B3_D0T1_TopMax_VolcanoPlots_0p05_FC2]=lipidevaluate('20240423_B3_T_Day0_Rep3_600-900_TopMax.xlsx',0.05,2);%

%input file name in ''. File must be the xls in the same directory as the
%script
%the script outputs the following matrices: 
% OutputMatrix in the format: lipid m/z, tumour avg spectrum, tumour spetrum std, tumour Number of pixels, grey avg, grey stg, grey N, white avg, white std, white N
% The spectra of the segemented tumour, grey and white matter regions.
% VolcanoPlots: lipid m/z, log2 of fold change Tumour vs Grey; -log10 of p_value Tumour vs Grey;og2 of fold change Tumour vs White; -log10 of p_value Tumour vs White;

filename = inputfile;
outDir = fileparts(mfilename('fullpath'));

if isempty(outDir)
    outDir = pwd;
end

if ~isfile(filename)
    error('File not found: %s', filename);
end

try
    T = readtable(filename);
catch ME
    error('Unable to open workbook ''%s'': %s', filename, ME.message);
end

coord = [T.x T.y];
rotated_coord = coord;
spectra = table2array(T(:,5:end));

lipidNames = T.Properties.VariableNames(5:end);

lipidNames = erase(lipidNames,"x");
lipidNames = strrep(lipidNames,"_", ".");

lipids = str2double(lipidNames);

disp(T.Properties.VariableNames(1:10))

disp(lipids(1:10))

MatLength = height(T);
MatWidth  = width(T);

%calculate dimensions for resizing
unique_x = unique(coord(:,1));
unique_y = unique(coord(:,2));

Nx = numel(unique_x);
Ny = numel(unique_y);

hyperimage = zeros(Nx,Ny,size(spectra,2));

for i = 1:size(coord,1)

    x = find(unique_x == coord(i,1));
    y = find(unique_y == coord(i,2));

    hyperimage(x,y,:) = spectra(i,:);

end

disp(max(hyperimage(:)))
disp(min(hyperimage(:)))


%identify the matrices relative to the representative peaks
TumourPeak=682.59;
GreyPeak=600.51;
WhitePeak=888.62;
[~,TumourPeakPos]=min(abs(lipids - TumourPeak));
[~,WhitePeakPos]=min(abs(lipids - WhitePeak));
[~,GreyPeakPos]=min(abs(lipids - GreyPeak));

disp("Tumour")
disp(lipids(TumourPeakPos))

disp("Grey")
disp(lipids(GreyPeakPos))

disp("White")
disp(lipids(WhitePeakPos))

%normalise the three peaks for image display and segmentation
TIC=sum(hyperimage,3);
hyperimageN=hyperimage./TIC;
hyperimageN(isnan(hyperimageN))=0;

TumourImN=hyperimageN(:,:,TumourPeakPos)/max(max(hyperimageN(:,:,TumourPeakPos)));
WhiteImN=hyperimageN(:,:,WhitePeakPos)/max(max(hyperimageN(:,:,WhitePeakPos)));
GreyImN=hyperimageN(:,:,GreyPeakPos)/max(max(hyperimageN(:,:,GreyPeakPos)));

disp(max(TumourImN(:)))
disp(max(GreyImN(:)))
disp(max(WhiteImN(:)))

%display three-colour peaks image
imageBrain(:,:,1)=TumourImN;
imageBrain(:,:,2)=WhiteImN;
imageBrain(:,:,3)=GreyImN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotate Day5 images for manual segmentation only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isDay5 = contains(lower(filename),"day5");

if isDay5

    disp("Day5 detected -> rotating all spatial data")

    angle = 90 - 11.25;

    imageBrain = imrotate(imageBrain, angle, 'bilinear', 'crop');

    TumourImN = imrotate(TumourImN, angle, 'bilinear', 'crop');
    GreyImN   = imrotate(GreyImN, angle, 'bilinear', 'crop');
    WhiteImN  = imrotate(WhiteImN, angle, 'bilinear', 'crop');


    % rotate lipid cube
    hyperimageN_rot = zeros(size(hyperimageN));

    for k = 1:size(hyperimageN,3)

        hyperimageN_rot(:,:,k) = imrotate( ...
            hyperimageN(:,:,k), ...
            angle, ...
            'nearest', ...
            'crop');

    end

    hyperimageN = hyperimageN_rot;


    % rotate coordinates
    centre = [
        mean(unique_x)
        mean(unique_y)
    ]';

    coords_centered = coord - centre;

    theta = deg2rad(angle);

    R = [
        cos(theta) -sin(theta);
        sin(theta) cos(theta)
    ];

    rotated_coord = coords_centered * R';

    rotated_coord = rotated_coord + centre;

    end

% Create a figure for the message box
h = msgbox('Draw a polygon around the replica region exlcuding noise. Press ok to continue...');

% Wait for the user to press the OK button
uiwait(h);

binaryMask = roipoly(TumourImN);
binaryMask = logical(binaryMask);
if isempty(binaryMask)
    binaryMask = false(size(TumourImN)); % default empty mask matching image size
elseif ~isequal(size(binaryMask), size(TumourImN))
    binaryMask = imresize(binaryMask, size(TumourImN), 'nearest');
end

%define segmented regions
TumourThr=0.1;
GreyThr=0.1;
WhiteThr=0.1;
TumourSeg = (TumourImN > GreyImN) & (TumourImN > WhiteImN) & (TumourImN >TumourThr);
GreySeg = (GreyImN > TumourImN) & (GreyImN > TumourImN) & (GreyImN >GreyThr);
WhiteSeg = (WhiteImN > GreyImN) & (WhiteImN > TumourImN) & (WhiteImN >WhiteThr);

TumourSeg= TumourSeg.*binaryMask;
GreySeg = GreySeg.*binaryMask;
WhiteSeg = WhiteSeg.*binaryMask;


%erode regions
se = strel('disk', 1);
TumourSegE=imerode(TumourSeg,se);
GreySegE=imerode(GreySeg,se);
WhiteSegE=imerode(WhiteSeg,se);

%% Save segmentation masks

disp("========== REACHED SAVE SECTION ==========")
disp("Saving outputs to:")
disp(outDir)


imwrite(uint8(TumourSegE)*255,...
    fullfile(outDir,'tumour_mask.png'));

imwrite(uint8(GreySegE)*255,...
    fullfile(outDir,'grey_mask.png'));

imwrite(uint8(WhiteSegE)*255,...
    fullfile(outDir,'white_mask.png'));

% Region mask (uint8 labels)
regionMask = zeros(size(TumourSegE), 'uint8');
TumourSegE = logical(TumourSegE);
GreySegE   = logical(GreySegE);
WhiteSegE  = logical(WhiteSegE);
regionMask(TumourSegE) = 1;
regionMask(GreySegE)   = 2;
regionMask(WhiteSegE)  = 3;

% Save MAT, PNG and CSV to cwd
save(fullfile(outDir, 'regionMask.mat'), 'regionMask', '-v7.3');
imwrite(regionMask, fullfile(outDir, 'regionMask.png'));
writematrix(regionMask, fullfile(outDir, 'regionMask.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export pixel-by-pixel region labels and feature matrix rotated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Build grid coordinates that match regionMask ordering
[Xg,Yg] = meshgrid(unique_y, unique_x);
if ~isequal(size(regionMask), size(Xg))
    Xg = Xg';
    Yg = Yg';
end
coords_grid = [Xg(:), Yg(:)];   % Npix x 2

% If Day5, apply the same rotation used earlier to the grid coordinates
if isDay5
    centre = [mean(unique_x), mean(unique_y)];
    coords_centered = coords_grid - centre;
    theta = deg2rad(angle);                         % angle defined earlier
    R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    rotated_coords_grid = (coords_centered * R') + centre;
else
    rotated_coords_grid = coords_grid;
end

% pixel_regions table (same ordering as regionMask)
pixel_regions = table(...
    rotated_coords_grid(:,1), ...
    rotated_coords_grid(:,2), ...
    regionMask(:), ...
    'VariableNames', {'x','y','region'});

% Write pixel regions CSV
pixelFile = fullfile(outDir,'pixel_regions.csv');
writetable(pixel_regions, pixelFile);

% Build feature_matrix_rotated matching pixel order (no region column)
Nbands = size(hyperimage,3);
spectra_grid = reshape(hyperimage, [], Nbands);   % ordering matches Xg(:),Yg(:)

% Create output table with sample/day, x, y and spectra columns (no region)
Npix = size(rotated_coords_grid,1);
% Determine metadata columns from original table T (exclude x,y and spectra)
metaNames = T.Properties.VariableNames(1:4);        % adjust if metadata in different cols
metaNames = setdiff(metaNames, {'x','y'}, 'stable');

% Build table with metadata repeated per pixel
metaTable = table();
for m = 1:numel(metaNames)
    val = T{1, metaNames{m}};                        % assume file-level metadata same for all rows
    if iscell(val), val = val{1}; end
    metaTable.(metaNames{m}) = repmat(val, Npix, 1);
end

% Base coords table
coordsTable = table(rotated_coords_grid(:,1), rotated_coords_grid(:,2), ...
    'VariableNames', {'x','y'});

% Spectra columns
Nbands = size(hyperimage,3);
spectra_grid = reshape(hyperimage, [], Nbands);   % ordering matches Xg(:),Yg(:)

% Build final feature table: [meta | x y | spectra]
featureTable = [metaTable, coordsTable];

% Use lipidNames for spectra column names (truncate/pad to match Nbands)
for k = 1:min(Nbands, numel(lipidNames))
    colname = lipidNames{k};
    featureTable.(colname) = spectra_grid(:,k);
end

% If there are more bands than lipidNames, name them band1, band2, ...
if Nbands > numel(lipidNames)
    for k = numel(lipidNames)+1:Nbands
        colname = sprintf('band%d', k);
        featureTable.(colname) = spectra_grid(:,k);
    end
end

writetable(featureTable, fullfile(outDir,'feature_matrix_rotated.csv'));


disp("Saved pixel_regions.csv and feature_matrix_rotated.csv")

%for k = 1:Nbands
 %   featureTable.(featureVarNames{3+k}) = spectra_grid(:,k);
%end

writetable(featureTable, fullfile(outDir,'feature_matrix_rotated.csv'));

disp("Saved pixel_regions.csv and feature_matrix_rotated.csv")



%display eroded regions
figure()
subplot(2,2,1)
imshow(imageBrain);
title('overlay');
subplot(2,2,2)
imshow(TumourSegE);
title('Tumour Segment')
subplot(2,2,3)
imshow(GreySegE);
title('Grey Segment')
subplot(2,2,4)
imshow(WhiteSegE);
title('White Segment')

%calculate average spectra segmented regions

TumourSp= hyperimageN.*TumourSegE;
GreySp= hyperimageN.*GreySegE;
WhiteSp= hyperimageN.*WhiteSegE;

TumourSpN=sum(sum(TumourSegE,1),2);
GreySpN=sum(sum(GreySegE,1),2);
WhiteSpN=sum(sum(WhiteSegE,1),2);


TumourSpSum = sum(sum(TumourSp,1),2);
GreySpSum = sum(sum(GreySp,1),2);
WhiteSpSum = sum(sum(WhiteSp,1),2);

TumourSpAvg = squeeze(TumourSpSum./TumourSpN);
GreySpAvg = squeeze(GreySpSum./GreySpN);
WhiteSpAvg = squeeze(WhiteSpSum./WhiteSpN);

%calculate peak std for segmented regions

Tumour_sq_diff = (hyperimageN - reshape(TumourSpAvg,[1,1,length(TumourSpAvg)])).^2;
Tumour_sum_sq_diff = sum(sum(Tumour_sq_diff .* TumourSegE, 1), 2);
TumourSpStd = sqrt(squeeze(Tumour_sum_sq_diff ./ TumourSpN));

Grey_sq_diff = (hyperimageN - reshape(GreySpAvg,[1,1,length(GreySpAvg)])).^2;
Grey_sum_sq_diff = sum(sum(Grey_sq_diff .* GreySegE, 1), 2);
GreySpStd = sqrt(squeeze(Grey_sum_sq_diff ./ GreySpN));

White_sq_diff = (hyperimageN - reshape(WhiteSpAvg,[1,1,length(WhiteSpAvg)])).^2;
White_sum_sq_diff = sum(sum(White_sq_diff .* WhiteSegE, 1), 2);
WhiteSpStd = sqrt(squeeze(White_sum_sq_diff ./ WhiteSpN));

%compile the final matrix

outputMatrix(:,1)=lipids;
outputMatrix(:,2)=TumourSpAvg;
outputMatrix(:,3)=TumourSpStd;
outputMatrix(:,4)=TumourSpN;
outputMatrix(:,5)=GreySpAvg;
outputMatrix(:,6)=GreySpStd;
outputMatrix(:,7)=GreySpN;
outputMatrix(:,8)=WhiteSpAvg;
outputMatrix(:,9)=WhiteSpStd;
outputMatrix(:,10)=WhiteSpN;

figure()
subplot(3,1,1)
plot(outputMatrix(:,1),outputMatrix(:,2))
title('Tumour Avg')
subplot(3,1,2)
plot(outputMatrix(:,1),outputMatrix(:,5))
title('Grey Avg')
subplot(3,1,3)
plot(outputMatrix(:,1),outputMatrix(:,8))
title('White Avg')

% Temp fix for the p-values
WhiteSpN=50;
TumourSpN=50;
GreySpN=50;

% Calculate log2 fold change
log2TvG = log2(TumourSpAvg ./ GreySpAvg);
log2TvW = log2(TumourSpAvg ./ WhiteSpAvg);

% Calculate t-statistic and p-values for each variable
t_valuesTvG = (TumourSpAvg - GreySpAvg) ./ sqrt((TumourSpStd.^2 ./ TumourSpN) + (GreySpStd.^2 ./ GreySpN));
p_valuesTvG = 2 * (1 - tcdf(abs(t_valuesTvG), min(TumourSpN, GreySpN) - 1));

t_valuesTvW = (TumourSpAvg - WhiteSpAvg) ./ sqrt((TumourSpStd.^2 ./ TumourSpN) + (WhiteSpStd.^2 ./ WhiteSpN));
p_valuesTvW = 2 * (1 - tcdf(abs(t_valuesTvW), min(TumourSpN, WhiteSpN) - 1));


%temp p-value fix
p_valuesTvW(p_valuesTvW==0)=0.0000000000000001;
p_valuesTvG(p_valuesTvG==0)=0.0000000000000001;
% Generate Volcano Plots
figure;
scatter(log2TvG, -log10(p_valuesTvG), 'o', 'filled');
xlabel('Log2 Fold Change');
ylabel('-Log10 p-value');
title('GvT Volcano Plot');
grid on;

% Optional: Add a significance threshold line
hold on;
yline(-log10(p_val_thres), '--r'); % For a p-value threshold of 0.05
xline(log2(fold_thres), '--r'); % For a p-value threshold of 0.05
hold off;

figure;
scatter(log2TvW, -log10(p_valuesTvW), 'o', 'filled');
xlabel('Log2 Fold Change');
ylabel('-Log10 p-value');
title('WvT Volcano Plot');
grid on;

% Optional: Add a significance threshold line
hold on;
yline(-log10(p_val_thres), '--r'); % For a p-value threshold of 0.05
xline(log2(fold_thres), '--r'); % For a p-value threshold of 0.05
hold off;

VolcanoPlots(:,1)=lipids;
VolcanoPlots(:,2)=log2TvG;
VolcanoPlots(:,3)=-log10(p_valuesTvG);
VolcanoPlots(:,4)=log2TvW;
VolcanoPlots(:,5)=-log10(p_valuesTvW);

save(fullfile(outDir,'lipid_results.mat'), ...
    'outputMatrix',...
    'TumourSp',...
    'GreySp',...
    'WhiteSp',...
    'VolcanoPlots',...
    '-v7.3');

disp("Saved lipid_results.mat")