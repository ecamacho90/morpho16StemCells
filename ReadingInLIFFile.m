%% Set the path
clear all;
dataDir = '/Users/Granton/Dropbox/StemCellProjectWarmflash/Data/08_02_16 _BMP_LDN1h/';
filename = fullfile(dataDir, '08_02_16_BMP_LDN1h.lif');


%% Read the LIF file
r = bfGetReader(filename);

% get number of series
numSeries = r.getSeriesCount();

%% Converting to TIFF

for s=9:numSeries % Set the projects in the lif file that will be converted to TIFF
    
%s = 12; 

r.setSeries(s - 1);

% allocate empty array for data
imclass = class(bfGetPlane(r, 1));
stackSize = [r.getSizeY(), r.getSizeX(), r.getSizeZ(), r.getSizeC(), r.getSizeT()];
data = zeros(stackSize, imclass);

% Save the data in the array

for i = 1:r.getImageCount()
    ZCTidx = r.getZCTCoords(i-1) + 1;
    fprintf('.');
    if rem(i,80) == 0
        fprintf('\n');
    end
    data(:,:, ZCTidx(1), ZCTidx(2), ZCTidx(3)) = bfGetPlane(r, i);
end

fprintf('\n');

% Allocate empty array for the TIFF data
stackSizeMIP = [r.getSizeY(), r.getSizeX(), r.getSizeT()];
sizex=r.getSizeX();
sizey=r.getSizeY();
sizet=r.getSizeT();
C1MIPdata = zeros(stackSizeMIP, imclass);
C2MIPdata = zeros(stackSizeMIP, imclass);


    for k = 1:sizet  % Loop over the time points
        for i = 1:sizey  % Loop over the y points
            for j = 1:sizex  % Loop over the x points

                C1MIPdata(i,j,k) = max(data(i,j,:,1,k)); %Maximum intensity of the mask1
                C2MIPdata(i,j,k) = max(data(i,j,:,2,k)); %Maximum intensity of the mask2

            end
        end

        %Show the time point to check that it's running
        disp(k)
        
        % Set name of the TIFF file for the mask for smad4 activity of project s
        outputFileName = sprintf(['S4_','08_02_16_BMP_LDN1h%i.tiff'],s);
        
        % Append the mask for smad 4 activity at time t to the TIFF file
        % with the specified name in outputFileName
        imwrite(C1MIPdata(:,:,k),outputFileName,'WriteMode','append');
        
        % Set name of the TIFF file for the mask for the nuclei of project s
        outputFileName2 = sprintf(['Nuclei_','08_02_16_BMP_LDN1h%i.tiff'],s);
        
        % Append the mask for nuclei at time t to the TIFF file
        % with the specified name in outputFileName
        imwrite(C2MIPdata(:,:,k),outputFileName2,'WriteMode','append');
    end

end

% Show the first image to check that it worked
tidx = 1; imshow(C1MIPdata(:,:,tidx),[]); %end of reading in Leica
