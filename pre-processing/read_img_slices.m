function [imageVol]=read_img_slices(address)

%read image slices sequentially

files = dir(address);
num_img = numel(files);
infile1=files(1).name;
foldername= extractBefore(address,'*.');

firstImage=uint8(imread([foldername infile1]));

imageVol=uint8(zeros(size(firstImage,1),size(firstImage,2),num_img));

imageVol(:,:,1)=firstImage;

for k = 2:num_img
    infile1=files(k).name;
    imageVol(:,:,k)=uint8(imread([foldername infile1]));
end
