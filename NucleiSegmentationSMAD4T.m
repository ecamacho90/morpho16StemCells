%%%%%%%%% Analyse h5 data preprocessed in
%%%%%%%%% Ilastiks

% path to Ilastiks exported files:
datasetname='/exported_data';
% our home directory, find name it via command pwd
dataDir = '/Users/elenacamachoaguilar/Dropbox Warwick/Dropbox/StemCellProjectWarmflash/code';%
addpath(dataDir);

% READ FILES OF SEGMENTATION
% defines filenames of input data
filenameCells = fullfile(dataDir, 'ActivinBMP3rd_MIP_p0013_w0000_Simple Segmentation.h5');
filenameNuclei = fullfile(dataDir, 'ActivinBMP3rd_MIP_p0013_w0001_Simple Segmentation.h5');
% read segmentation data (squeeze removes excessive dimensions) and turn everything that was red into 1 in binary mode
foregroundLabel = 1;
% Segmentation definition, 1: nucleus, 2: background
H5 = hdf5read(filenameCells,datasetname);
cellsegmentation = squeeze(hdf5read(filenameCells,datasetname)) == foregroundLabel;
nucleisegmentation = squeeze(hdf5read(filenameNuclei,datasetname)) == foregroundLabel;

backgroundLabel = 2;
backgroundsegmentationG = squeeze(hdf5read(filenameCells,datasetname)) == backgroundLabel;
backgroundsegmentationR = squeeze(hdf5read(filenameNuclei,datasetname)) == backgroundLabel;

dimensions = size(cellsegmentation);
tmax = 1;%dimensions(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in raw data %%%%%%%%%%%%%%
% allocate empty array for data %

smad4raw = zeros([1024 1024 tmax],'uint16');
% read in raw data
for i = 1:tmax
    smad4raw(:,:,i) = imread('ActivinBMP3rd_MIP_p0013_w0000.tif', i)';
end

redchannelnuclei = zeros([1024 1024 tmax],'uint16');
% read in raw data
for i = 1:tmax
    redchannelnuclei(:,:,i) = imread('ActivinBMP3rd_MIP_p0013_w0001.tif', i)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Segmentation vizualisation 
% example time point to vizualise segmentation
t = 1;
% allocate empty array for data
R = nucleisegmentation(:,:,t);
G = cellsegmentation(:,:,t);
B = 0*G;
% show image for time t, cat concatenates nuclei and cell segmentation data
% to colour them differently
imshow(cat(3,R,G,B),[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% analyse raw tif data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% intensity_n_t will hold mean nuclear intensity of smad4 over time
intensity_n_t = zeros([1,tmax]);

% create a structure for cells for each time point
alltimes = {};

Cells = struct('x',{},'y',{},'R',{},'GN',{},'GC',{});

        %x,y are the coordinates of the centroid of each cell over time
        %R is the nuclear intensity in the red channel
        %GN is the nuclear intensity for smad4 (green channel)
        %GC is the cytoplasm intensity for smad4 (green channel)
        

for ti = 1:tmax
    
    alltimes{ti} = Cells;
    
end

    % Steps to improve image quality of nuclei red channel image
    % ---------------------------------------------
    nuclei_new = nucleisegmentation;
    
    % use bwareaopen to exclude noisy bright bits
    nuclei_new = bwareaopen(nuclei_new,5000);
    
    % use imaopen and imclose with the structuring element se to 
    % make edges smooth. Chose a disc of radius 10.
    se = strel('disk',10);
    nuclei_new = imopen(nuclei_new(:,:,1), se);
    
    opbackground = imdilate(nuclei_new, se);
    background = not(opbackground);
    Indback = find(background);
    backgroundintR = mean(redchannelnuclei(Indback));
    
    opbackground = imdilate(cellsegmentation(:,:,1), se);
    background = not(opbackground);
    Indback = find(background);
    backgroundintG = mean(smad4raw(Indback));

%% loop comment later
for t = 1:tmax
    
    redchannelnuclei(:,:,t) = redchannelnuclei(:,:,t) - backgroundintR;
    
    smad4raw(:,:,t) = smad4raw(:,:,t) - backgroundintG;
    
    % Option to use imclose or imclearboarders
    % se = strel('disk',1);
    % nuclei_new = imclose(nuclei_new, se);    
    % nuclei_new = imclearborder(nuclei_new,8);

    % get properties of transformed segmentation
    stats = regionprops(nuclei_new(:,:,t), 'PixelIdxList', 'Centroid');
    
    %stats(i): the ith nucleus identified on image. Then repeat for all
    %nculei identified to get average intensity accross sample
    intensity_ni = zeros(size(stats));
    intensity_n_t(t) = mean(intensity_ni);
    
    num = size(stats);
    num_nuc = num(1);
    
    cytmask = false(1024);
    
    for k = 1:num_nuc
        
        Cells(k).x = stats(k).Centroid(1);
        Cells(k).y = stats(k).Centroid(2);
        
        intensity_ni_R(k) = mean(redchannelnuclei(stats(k).PixelIdxList));
        Cells(k).R = intensity_ni_R(k);
        
        intensity_ni(k) = mean(smad4raw(stats(k).PixelIdxList));
        Cells(k).GN = intensity_ni(k);
        
        cytmask = false(1024);
        cytmask(stats(k).PixelIdxList) = true;
        donut1 = imdilate(cytmask,strel('disk',20))-cytmask;
        donut = donut1 & cellsegmentation(:,:,t);
    
        imshow(donut)
        intensity_cyt = mean(smad4raw(donut));
        Cells(k).GC = intensity_cyt;
    end
    alltimes{t} = Cells;  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cytoplasm segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%
cytoplasmsegmentation= cellsegmentation&(~nucleisegmentation); %Cells minus the nuclei

%% 

