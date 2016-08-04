%% Set the path
clear all;
dataDir = '/Users/Granton/Dropbox/StemCellProjectWarmflash/Data/08_02_16 _BMP_LDN1h/';
filename = fullfile(dataDir, '08_02_16_BMP_LDN1h.lif');


%% Read the file
r = bfGetReader(filename);

% get number of series
numSeries = r.getSeriesCount();

for s=9:numSeries
    
%s = 12; 

r.setSeries(s - 1);
% allocate empty array for data
imclass = class(bfGetPlane(r, 1));
stackSize = [r.getSizeY(), r.getSizeX(), r.getSizeZ(), r.getSizeC(), r.getSizeT()];
data = zeros(stackSize, imclass);
for i = 1:r.getImageCount()
    ZCTidx = r.getZCTCoords(i-1) + 1;
    fprintf('.');
    if rem(i,80) == 0
        fprintf('\n');
    end
    data(:,:, ZCTidx(1), ZCTidx(2), ZCTidx(3)) = bfGetPlane(r, i);
end
fprintf('\n');
stackSizeMIP = [r.getSizeY(), r.getSizeX(), r.getSizeT()];
sizex=r.getSizeX();
sizey=r.getSizeY();
sizet=r.getSizeT();
C1MIPdata = zeros(stackSizeMIP, imclass);
C2MIPdata = zeros(stackSizeMIP, imclass);

for k = 1:sizet
for i = 1:sizey
    for j = 1:sizex
        C1MIPdata(i,j,k) = max(data(i,j,:,1,k));
        C2MIPdata(i,j,k) = max(data(i,j,:,2,k));
    end
end
disp(k)
outputFileName = sprintf(['S4_','08_02_16_BMP_LDN1h%i.tiff'],s);
imwrite(C1MIPdata(:,:,k),outputFileName,'WriteMode','append');
outputFileName2 = sprintf(['Nuclei_','08_02_16_BMP_LDN1h%i.tiff'],s);
imwrite(C2MIPdata(:,:,k),outputFileName2,'WriteMode','append');
end

end
%tiffobj = Tiff('C1MIP8-2Data.tif','w');
%write(tiffobj,C1MIPdata);
tidx = 1; imshow(C1MIPdata(:,:,tidx),[]); %end of reading in Leica
