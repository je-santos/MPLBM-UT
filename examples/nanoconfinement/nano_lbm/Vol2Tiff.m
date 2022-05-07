function [] = Vol2Tiff(imageVol, folder, outfoldername)

mkdir([folder '/' outfoldername]);
sizeK = size(imageVol,3);

for k = 1:sizeK
    imwrite(imageVol(:,:,k), ...
        [folder '/' outfoldername '/' outfoldername '_' num2str(k,'%0.4i')], ...
        'tiff');
end

end
