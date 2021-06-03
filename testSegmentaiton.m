%% test

data1 = a.polymer(:,:,:,1);
data2 = a.polymer(:,:,:,3);


%% Segmentation test1

%BW1 = imbinarize(data1);


[bg,sig] = getCornerStats(data1);

tHold = tholdSigBg(bg,sig);
BW1 = data1>2*tHold;

for i = 1:size(BW1,3)
    BW1(:,:,i) = bwareaopen(BW1(:,:,i),50);
    
end

se = strel('sphere',3);
BW1 = imclose(BW1,se);

figure
imagesc(BW1(:,:,69))


%%
[bg,sig] = getCornerStats(data2);

tHold = tholdSigBg(bg,sig);
BW2 = data2>2*tHold;

for i = 1:size(BW2,3)
    BW2(:,:,i) = bwareaopen(BW2(:,:,i),50);
    
end

se = strel('sphere',3);
BW2 = imclose(BW2,se);

figure
imagesc(BW2(:,:,69))



if sum(BW1(:)) > sum(BW2(:))
    
else
    error('oups');
end

%% Test2 

[bg,sig] = getCornerStats(data1);

tHold = tholdSigBg(bg,sig);

mask1 = stdfilt(data1, ones(15,15,3)) > tHold;
mask2 = data1 > 2*tHold;
% Then you'll undoubtedly have to do some clean up, like with bwareafilt() or imclose().
% Now create final mask.
mask = mask1 & mask2;
se = strel('sphere',5);
mask = imclose(mask,se);

mask = mask - test.cellMask(:,:,:,1);
figure
imagesc(mask(:,:,69))



%%
[bg,sig] = getCornerStats(data2);

tHold = tholdSigBg(bg,data2(:));

mask1 = stdfilt(data2, ones(15,15,3)) > tHold;
mask2 = data2 > 2*tHold;
% Then you'll undoubtedly have to do some clean up, like with bwareafilt() or imclose().
% Now create final mask.
mask2 = mask1 & mask2;
se = strel('sphere',5);
mask2 = imclose(mask2,se);

mask2 = mask2-test.cellMask(:,:,:,3);

figure
imagesc(mask2(:,:,69))




function [tHold] = tholdSigBg(bg,sig)
%THOLDSIGBG finds optimal treshhold that separates two distributions one
%comming from background and another from signal
%   Detailed explanation goes here

sigV = sig(:);
bgV  = bg(:);
assert(mean(bgV)<mean(sigV),'strange your background looks larger than your signal')

% find CDF of signal and CCDF of background
[sigCDF] = Plotting.getCDF(sigV);
[~, bgCCDF]  = Plotting.getCDF(bgV);

% find optimal tHold that separates two distributions
allX = [sigCDF.x; bgCCDF.x];
minX = min(allX);
maxX = max(allX);
globXq = linspace(double(minX),double(maxX),length(allX)*3);
sigVq = interp1(double(sigCDF.x),double(sigCDF.y),globXq);
bgVq  = interp1(double(bgCCDF.x),double(bgCCDF.y),globXq);
diffVq = abs(sigVq-bgVq);
%not sure what was this for:
% if mean(sigCDF.x(sigCDF.y>0.95))>2*max(bg)
%     tHold = max(bg);
% else
[~,idx] = nanmin(diffVq);
tHold = globXq(idx);
%end

end

 function [cornerMat,innerMat] = getCornerStats(currData)
   %get threshold based on corner
    height10 = round(size(currData,1)/10);
    width10  = round(size(currData,2)/10);



    %get corner
    cornerMat1 = currData(1:height10,1:end,:);
    cornerMat2 = currData(end-height10+1:end,1:end,:);
    cornerMat3 = currData(1:end,1:width10,:);
    cornerMat4 = currData(1:end,end-width10:end,:);
    cornerMat = [cornerMat1(:);cornerMat2(:);cornerMat3(:);cornerMat4(:)];
%             cornerMat(end-height10:end,1:width10,:) = currData(end-height10:end,1:width10,:);
%             cornerMat(1:height10,end-width10:end,:) = currData(1:height10,end-width10:end,:);
%             cornerMat(end-height10:end,end-width10:end,:) = currData(end-height10:end,end-width10:end,:);

    innerMat = currData(height10:end-height10,width10:end-width10,:);

end