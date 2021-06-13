clear all
% url = ‘https://ndownloader.figshare.com/files/26751203’;
% Go to this address to download this and save it where you can find it.
fName = ‘SupFig1c_BG_MAX_Cell01.tif’;
%% get size of the tiff
tiffInfo = imfinfo(fName); % return tiff structure, one element per image
nC = 3;
nT = length(tiffInfo)/nC
nY = tiffInfo(1).Height
nX = tiffInfo(1).Width
%% Load and make the image tensor
IMG = zeros(nT,nY,nX,nC,‘double’);
for iT = 1:nT
  for iC = 1:nC
    kTiff = (iT-1)*nC+iC;
    IMG(iT,:,:,iC) = imread(fName, kTiff);
  end
end
%% Rescale to have range [0:1]
maxRGB = squeeze(max(IMG,[],[1,2,3]));
minRGB = squeeze(min(IMG,[],[1,2,3]));
scaledIMG = 0*IMG;
for iC = 1:3
  scaledIMG(:,:,:,iC) = (IMG(:,:,:,iC)-minRGB(iC))/(maxRGB(iC)-minRGB(iC));
end
%% Show the 22 time point:
figure(1);
rangeX = [100 150];
rangeY = [100 150];
plotImage(scaledIMG)
%% Make a movie:
figure(2);
vidObj = VideoWriter(‘CellVideos.avi’);
open(vidObj);
for iT = 1:nT
  %% Crop around the brightest spot.
  [~,Jmax] = max(IMG(iT,:,:,3),[],[2,3],‘linear’);
  [JY,JX] = ind2sub([size(IMG,2),size(IMG,3)],Jmax);
  rangeX = JX+[-25 25];
  rangeY = JY+[-25 25];
  hold off
  plotImage(scaledIMG,iT,rangeX,rangeY)
  hold on
  plot(26,26,‘s’,‘MarkerSize’,20,‘MarkerEdgeColor’,[1 1 1])
  currFrame = getframe(gcf);
  writeVideo(vidObj,currFrame);
end
close(vidObj);
%% Function that makes all the plots for specified range of locations and specified time.
function plotImage(scaledIMG,iT,rangeX,rangeY)
arguments
  scaledIMG
  iT = 1;
  rangeX = [1,size(scaledIMG,2)];
  rangeY = [1,size(scaledIMG,3)];
end
for iC = 1:3
  subplot(1,4,iC);
  imshow(squeeze(scaledIMG(iT,rangeY(1):rangeY(2),rangeX(1):rangeX(2),iC)))
end
% Make color image with channels reordered
subplot(1,4,4);
colorOrder = [3,2,1];
imshow(squeeze(scaledIMG(iT,rangeY(1):rangeY(2),rangeX(1):rangeX(2),colorOrder)))
end