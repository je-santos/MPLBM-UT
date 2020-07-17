% Script to create videos from .gif files


workingDir = 'tmp';
mkdir(workingDir)

imageNames = dir(['rho_f1_y_0','*.gif']);
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile(workingDir,'fluid1_flow.avi'));
outputVideo.FrameRate = 50;
open(outputVideo)

for ii = 1:length(imageNames)
   img = imread(fullfile(imageNames{ii}));
   writeVideo(outputVideo,img)
end

close(outputVideo)