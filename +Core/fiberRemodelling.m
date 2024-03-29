classdef fiberRemodelling < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        raw
        info
        channels
        results
    end
    
    methods
        function obj = fiberRemodelling(raw,info)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.raw = raw;
            obj.info = info;
        end
        
        function set.raw(obj,raw)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            assert(isstruct(raw),'raw is expected to be a structure');
            assert(and(isfield(raw,'path'),isfield(raw,'ext')),'raw is expected to be a structure with 2 fields: path and ext');
            assert(ischar(raw.path), 'Input path needs to be a char or string');
            assert(isfolder(raw.path),'Input path should be a folder containing one dataset to analyze');
            obj.raw = raw;
        end
        
        function loadData(obj,chan)
            path = obj.raw.path;
            ext  = obj.raw.ext;
            %let us check that there is no channel data existing
            if ~obj.existChannel
                disp('no channel data found, starting extraction ...')
                %get all file of appropriate extension in the file
                fileList = Core.fiberRemodelling.getFileInPath(path,ext);
                %get the different channel from the data
                channel  = obj.retrieveChannel(fileList,chan);
                                           
                %save the channel as matlab variable for the future
                filename = [obj.raw.path filesep 'channels.mat'];
                save(filename,'channel','-v7.3');
                disp('==========> DONE <==========')
            else
                disp('channel data found, loading from existing file ...')
                filename = [path filesep 'channels.mat'];
                tmp = load(filename);
                field = fieldnames(tmp);
                channel = tmp.(field{1});
                disp('==========> DONE <==========')
            end
         
            nFrames = size(channel.polymer,3);
            obj.raw.nFrames = nFrames;
            %store the data back into the object
            obj.channels = channel;
            
            
        end
        
        function channels = getChannel(obj)
            channels = obj.channels;
        end
        
        function results = getResults(obj)
           results = obj.results; 
        end
        
        function showChannel(obj)
            channel = obj.getChannel;
            field = fieldnames(channel);
            
            nField = length(field);
            sliceToShow = round(size(channel.cell,3)/2);
            frameToShow = round(size(channel.cell,4)/2);
            figure
            for i = 1:nField
                subplot(1,nField,i)
                currChan = channel.(field{i});
                try
                    imagesc(currChan(:,:,sliceToShow,frameToShow))
                catch
                    
                end
                axis image
                colormap('hot');
                title(field{i})
            end
            
        end
        
        function [allMask] = calc3DMask(obj)
            disp('Starting mask extraction...')
            %fileName = [obj.raw.path filesep 'SegmentedCell.tif'];
            
            disp('Performing segmentation on cell channel');
            chan = obj.getChannel;
            data1 = chan.cell;
            data2 = chan.nucleus;
            
            if isempty(data2)
                data = data1;
            else
                data = data1+data2;
            end
            allMask = zeros(size(data));
            for i = 1:size(data,4)
                disp(['Analyzing Frame ' num2str(i) ' / ' num2str(size(data,4))]);
                currData = data(:,:,:,i);
                % med filt
                
                currData = medfilt3(currData);
                
                % Gaussian filtering
                S = 2;
                % size of pixel in z vs x/y
                pixZ  = 4;
                zFactor = 2;
                sigma = [S,S,S*zFactor/pixZ];
                IMs = imgaussfilt3(currData, sigma);
                disp('DONE with filtering ------------')

                gBW = imbinarize(IMs,'adaptive','Sensitivity',0.2);
                
                gBW = bwareaopen(gBW,100);
                %Segmentation
              %  gBW = imbinarize(IMs);
              
%                 se = strel('disk',10);
%                 gBW = imopen(gBW,se);

                se = strel('disk',5);
                gBW = imclose(gBW,se);
                
                %find largest BW
                prop = regionprops3(gBW,'Volume','VoxelIDXList');
                
                voxelIdx = prop.VoxelIdxList{prop.Volume==max(prop.Volume)};
                gBW = zeros(size(gBW));
                gBW(voxelIdx)= 1;
                
                gBW = imfill(gBW);
                %storing temporary results             
            %    dataStorage.BinaryTiff(fileName,gBW);

                disp('Extracting Contour')
                %here we obtain the cell contour
                contour = obj.getCellContour(gBW);
                obj.results.cellContour{i} = contour;

                allMask(:,:,:,i) = gBW;
            
            end
            obj.results.cellMask = logical(allMask);       
            
        end
        
        function plotCellContour(obj,idx)
            assert(isfield(obj.results,'cellContour'),'No cell mask found run calc3DMask first');
            
            if nargin<2
                idx = 1;
            end
            
            contour = obj.results.cellContour{idx};
            figure(2)
            hold on
            for i=1:length(contour)
                if ~isempty(contour{i})
                    for j=1:length(contour{i})
                        plot3(contour{i}{j}(:,2),contour{i}{j}(:,1),repmat(i,1,length(contour{i}{j}(:,1))),'k')
                    end
                end
            end

        end
        
        function renderCell3D(obj,idx)
            assert(isfield(obj.results,'cellMask'),'No cell mask found run calc3DMask first');
        
            if nargin <2
                idx = 1;
            end
            
            data2Render = obj.results.cellMask(:,:,:,idx);
            iSurface = isosurface(data2Render,1/2);
            
            % smoothing using compiled c code
         %   smoothISurface = rendering3D.smoothpatch(iSurface,0,10);
            %comnvert to px
          %  smoothISurface.vertices(:,1) = (smoothISurface.vertices(:,1));
           % smoothISurface.vertices(:,2) = (smoothISurface.vertices(:,2));
            %smoothISurface.vertices(:,3) = (smoothISurface.vertices(:,3));

            %% Displaying network model
            %z-coloring
            colorModel = iSurface.vertices(:,3)/max(iSurface.vertices(:,3));
            zColor = true;
           
            %Plot the network with Z coloring or unique color depending on the user
            %input
            figure
            if zColor
                p = patch('Faces',iSurface.faces,'Vertices',iSurface.vertices,'FaceVertexCData',colorModel,'FaceColor','interp');
                colormap('jet')
                p.EdgeColor = 'none';
                daspect([2 2 1])
                view(3);
                axis tight
                camlight
                lighting gouraud
                title('Z-coloring')
            else
                p2 = patch(iSurface);
                p2.FaceColor = colorModel;
                p2.EdgeColor = 'none';
                view(3);
                axis tight
                camlight
                lighting gouraud
                title('unicolor');
            end

            
            
        end
        
        function getDensifiedNetwork(obj)
            assert(isfield(obj.channels,'polymer'));
            assert(isfield(obj.results,'cellMask'));
            method = 'tHold';% 'std' tHold
            disp('Extracting densified network...')
            
            disp('Performing segmentation on polymer channel');
            chan = obj.getChannel;
            data = chan.polymer;
            
            allMask = zeros(size(data));
            for i = 1:size(data,4)
           
                disp(['Analyzing Frame ' num2str(i) ' / ' num2str(size(data,4))]);
                currData = data(:,:,:,i);
                currData = medfilt3(currData);
                
                % Gaussian filtering
                S = 2;
                % size of pixel in z vs x/y
                pixZ  = 4;
                zFactor = 2;
                sigma = [S,S,S*zFactor/pixZ];
                IMs = imgaussfilt3(currData, sigma);
                disp('DONE with filtering ------------')

                %get threshold based on corner
                [stat,cornerMat] = obj.getCornerStats(IMs);
                
                tHold = Misc.tholdSigBg(cornerMat,IMs);
                
                if strcmp(method,'tHold')
                    gBW = IMs>2*tHold;

                    for j = 1:size(gBW,3)
                        gBW(:,:,j) = bwareaopen(gBW(:,:,j),50);

                    end

                    se = strel('sphere',3);
                    gBW = imclose(gBW,se);

                elseif strcmp(method,'std')
         
                    mask1 = stdfilt(IMs, ones(15,15,3)) > tHold;
                    mask2 = IMs > 2*tHold;
                    % Then you'll undoubtedly have to do some clean up, like with bwareafilt() or imclose().
                    % Now create final mask.
                    gBW = mask1 & mask2;
                    se = strel('sphere',5);
                    gBW = imclose(gBW,se);
                    
                end

                % binarize the image based on max value in corners
                %OLD METHOD SEE ISSUE IN THE PPT detection problem
              %  gBW = currData>stat.max;                              
                %get only biggest one
                prop = regionprops3(gBW,'Volume','Voxelidxlist');
                
                voxelIdx = prop.VoxelIdxList{prop.Volume==max(prop.Volume)};
                gBW = zeros(size(gBW));
                gBW(voxelIdx)= 1;
                
%                 se = strel('sphere',3);
%                 gBW = imclose(gBW,se);
                
             %   gBW = imfill(gBW);
                
                disp('Extracting Contour')
                %here we obtain the cell contour
                contour = obj.getCellContour(gBW);
                obj.results.polymerContour{i} = contour;

         
                %% building 3D mask
            
                currCellMask = obj.results.cellMask(:,:,:,i);
                finalMask = gBW-currCellMask;
                finalMask(finalMask<0)=0;
                allMask(:,:,:,i) = finalMask;

            end
            
            obj.results.polymerMask = logical(allMask);
            
        end
                
        function renderCellPolymer3D(obj,idx)
            assert(isfield(obj.results,'cellMask'),'No cell mask found run calc3DMask first');
            assert(isfield(obj.results,'polymerMask'),'No polymer mask found run calc3DMask first');
            
            if nargin <2
                idx = 1;
            end
            
            data2Render = obj.results.cellMask(:,:,:,idx);
            iSurface = isosurface(data2Render,1/2);
            
            % smoothing using compiled c code
%             cellsmoothISurface = rendering3D.smoothpatch(iSurface,0,10);
%             %comnvert to px
%             cellsmoothISurface.vertices(:,1) = (cellsmoothISurface.vertices(:,1));
%             cellsmoothISurface.vertices(:,2) = (cellsmoothISurface.vertices(:,2));
%             cellsmoothISurface.vertices(:,3) = (cellsmoothISurface.vertices(:,3));
            figure(idx)
            hold on
            p2 = patch(iSurface);
            p2.FaceColor = [0,1,0];
            p2.EdgeColor = 'none';
            
            clear cellsmoothISurface;
                      
            data2Render = obj.results.polymerMask(:,:,:,idx);
            iSurface = isosurface(data2Render,1/2);
            
            % smoothing using compiled c code
%             polsmoothISurface = rendering3D.smoothpatch(iSurface,0,1);
%             %comnvert to px
%             polsmoothISurface.vertices(:,1) = (polsmoothISurface.vertices(:,1));
%             polsmoothISurface.vertices(:,2) = (polsmoothISurface.vertices(:,2));
%             polsmoothISurface.vertices(:,3) = (polsmoothISurface.vertices(:,3));

            %% Displaying network model
           
            
            p3 = patch(iSurface);
            p3.FaceColor = [0.7 0.7 0.7];
            p3.FaceAlpha = 0.5;
            p3.EdgeColor = 'none';
           
            view(3);
            axis tight
            camlight
            lighting gouraud
            title('unicolor');
            
        end
        
        function [stats] = calcStats(obj)
            assert(isfield(obj.results,'cellMask'),'No cell mask found run calc3DMask first');
            assert(isfield(obj.results,'polymerMask'),'No polymer mask found run calc3DMask first');
            
            cellMask = obj.results.cellMask;
            polymerMask = obj.results.polymerMask;
            cellChan = obj.channels.cell;
            polChan = obj.channels.polymer;
           
            %memory preallocation
            stats.cellVol = zeros(1,size(cellMask,4));
            stats.polVol  = zeros(1,size(cellMask,4));
            stats.cellInt = zeros(1,size(cellMask,4));
            stats.polInt = zeros(1,size(cellMask,4));
            
            %remove the cell from the polymer mask
            %CELL IS NOW REMOVED FROM THE POLYMER MASK IN
            %GETDENSIFIEDNETWORK
%             polymerMask=polymerMask-cellMask;
%             polymerMask(polymerMask<0) = 0;
%             polymerMask = logical(polymerMask);
            
            pxSize = obj.info.pxSizeXY/1000;
            pxSizeZ = obj.info.pxSizeZ/1000;
            
            voxelSize = pxSize^2*pxSizeZ;
            for i = 1:size(cellMask,4)
                currentCellM = cellMask(:,:,:,i);
                currentPolM  = polymerMask(:,:,:,i);
                
                currentCell = cellChan(:,:,:,i);
                currentPol  = polChan(:,:,:,i); 
                cellCorner  = Core.fiberRemodelling.getCornerStats(currentCell);
                polCorner   = Core.fiberRemodelling.getCornerStats(currentPol);
                
                %normalize the mean intensity to the pixel volumes and the 
                %intensity in the corner
                stats.cellInt(i) = mean(currentCell(currentCellM));
                stats.cellCornerInt(i) = cellCorner.median;
                stats.polInt(i)  = mean(currentPol(currentPolM));
                stats.polCornerInt(i)  = polCorner.median;
                if isnan(stats.polInt(i))
                    stats.polInt(i) =0;
                end
                stats.cellVol(i) = sum(currentCellM(:))*voxelSize;
                stats.polVol(i)  = sum(currentPolM(:))*voxelSize;
                
                
            end
            stats.voxelSize = voxelSize;
            %store stats in results
            obj.results.stats = stats;
            
            %Volume plot
            figure
                subplot(1,2,1)
                    plot(stats.cellVol)
                    ylim([0 1+1.5*max(stats.cellVol)])
                    axis square
                    xlabel('Time')
                    ylabel('Volume (\mum^3)')
                    title('Cell')

                subplot(1,2,2)
                    plot(stats.polVol)
                    ylim([0 1+1.2*max(stats.polVol)])
                    axis square
                    xlabel('Time')
                    ylabel('Volume (\mum^3)')
                    title('Densified polymer')
                    
                    
             figure
                subplot(1,2,1)
                    plot(stats.cellInt)
                    ylim([0 1.5*max(stats.cellInt)])
                    axis square
                    xlabel('Time')
                    ylabel('Intensity(a.u.)')
                    title('Cell')

                subplot(1,2,2)
                    plot(stats.polInt)
                    ylim([0 0.001+1.2*max(stats.polInt)])
                    axis square
                    xlabel('Time')
                    ylabel('Intensity(a.u.)')
                    title('Densified polymer')

        end
        
        function [NOP]    = calcNOP(obj)
            assert(isfield(obj.results,'polymerMask'),'No polymer mask found run calc3DMask first');
            error('not Implemented yet')
            polymerMask = obj.results.polymerMask;
            data        = obj.channels.polymer;
            gradx = zeros(size(data));
            grady = gradx;
            gradz = grady;
            
            for i = 1:size(data,4)
            
                [Gx,Gy,Gz]  = imgradientxyz(data(:,:,:,i));
                 
                gradx(:,:,:,i) = Gx;
                grady(:,:,:,i) = Gy;
                gradz(:,:,:,i) = Gz;
                
                tmpGrad = sqrt(Gx.^2+Gy.^2+Gz.^2);
                
                NOP.polymer{i} = tmpGrad(polymerMask(:,:,:,i)==1);
                NOP.bkg{i}     = tmpGrad(polymerMask(:,:,:,i)==0);
                
            end
            
            
            
            
            
        end
        
        function [AllDist] = calcDistCellPol(obj)
            assert(~isempty(obj.channels),'Data not found, please loadData() before using this function')
            assert(isfield(obj.results,'cellMask'),'No cell mask found run calc3DMask first');
            
            cellMask = obj.results.cellMask;
            polMask  = obj.results.polymerMask;
            AllDist  = cell(size(cellMask,4),1);
            for i = 1:size(cellMask,4)
                currentCellM = cellMask(:,:,:,i);
                currentPolM  = ~polMask(:,:,:,i);
                
                dist = obj.getDistBetweenMasks(currentPolM,currentCellM);
                
                AllDist{i} = dist;
               
            end
            
            if size(cellMask,4)>1
                %here we assume that the first time frame is 24h before
                %drug so it is supposed to be the "big mask" and need to be
                %inverted
                %the last frame is the small mask
                
                polM1 = ~polMask(:,:,:,1);
                polM2 = polMask(:,:,:,end);
                dist = obj.getDistBetweenMasks(polM1,polM2);
                obj.results.polDist = dist;
            end
                
            
            obj.results.distances = AllDist;
            
        end
        
        
        function intensityDistrib(obj,weight,step)
            
            assert(~isempty(obj.channels),'Data not found, please loadData() before using this function')
            assert(isfield(obj.results,'cellMask'),'No cell mask found run calc3DMask first');
            data = obj.channels.polymer;
            mask = obj.results.cellMask;
          
            figure
            hold on
            leg = cell(size(mask,4),1);
            intRes = cell(size(mask,4),1);
            
            for j = 1:size(mask,4)
                currentMask = mask(:,:,:,j);
                stat = obj.getCornerStats(data(:,:,:,j));
                
                EDM = DistMap.calcWeightedDistMap(currentMask,weight);

                edgeMin = min(EDM(:));
                edgeMax = max(EDM(:));

                binEdges = edgeMin:step:edgeMax;
                tmpRes = table(zeros(length(binEdges)-1,1),zeros(length(binEdges)-1,1),...
                    'VariableNames',{'Distance','Intensity'});
               

                %Calculate the intensity vs distance from the cell
                for i = 1: length(binEdges)-1

                    idx = and(EDM>=binEdges(i), EDM<binEdges(i+1));

                    currDistance = (binEdges(i)-binEdges(1) +binEdges(i+1))/2;
                    currentIntensity = mean(data(idx));

                    tmpRes.Distance(i,:) = currDistance;
                    tmpRes.Intensity(i,:) = currentIntensity;

                end
                
                tmpRes.normInt = tmpRes.Intensity-stat.mean;
                tmpRes.normInt = tmpRes.normInt/max(tmpRes.normInt);
                
                intRes{j} = tmpRes;
                
                plot(intRes{j}.Distance,intRes{j}.normInt(j,:));
                xlabel('Distance (nm)')
                ylabel('Average intensity per pixel')
                axis square
                box on
                hold on
                leg{j} = ['T' num2str(j)];

                disp('=====> DONE <=====');
            end
            obj.results.intRes = intRes;
            %save the intensity curve
            filename = [obj.raw.path filesep 'IntensityResults.mat'];
            save(filename,'intRes');
            legend(leg)
        
        end
    end
    methods (Access = private)
        
        function [channel] = retrieveChannel(~,fileList,chan)
            
            fields = fieldnames(chan);
            %check movie size
            [movieInfo] = Load.Movie.tif.getinfo([fileList(1).folder filesep fileList(1).name]);
            for i = 1: length(fields)
                Misc.multiWaitbar('Channel Extraction',i/length(fields));
                index2Channel  = contains(lower({fileList.name}),fields{i},'IgnoreCase',true);
                currChanList   = fileList(index2Channel);
                endFile = fileList(end).name;
                
                %get max frame and max z layer
                idx2T  = strfind(lower(fileList(1).name),'t0');
                if isempty(idx2T)
                    maxT = 1;
                else
                    maxT   = str2double(endFile(idx2T+1:idx2T+3));
                end
                
                idx2Z  = strfind(lower(fileList(1).name),'z0');
                if isempty(idx2Z)
                    maxZ = 1;
                else
                    maxZ   = str2double(endFile(idx2Z+1:idx2Z+3));
                end
                
                
                currChan = zeros(movieInfo.Length,movieInfo.Width,maxZ,maxT);
                for j = 1:length(currChanList)
                    currPath = [currChanList(j).folder, filesep, currChanList(j).name];
                    
                    nZ = str2double(currChanList(j).name(idx2Z+1:idx2Z+3));
                    nT = str2double(currChanList(j).name(idx2T+1:idx2T+3));
                    if isnan(nT)
                        nT = 1;
                    end
                    currChan(:,:,nZ,nT) = Load.Movie.tif.getframes(currPath,1);
                    
                    Misc.multiWaitbar('Frames',j/length(currChanList));
                    
                end
                tmp.(fields{i}) = currChan;
                
            end
            Misc.multiWaitbar( 'CloseAll' );
            cellIDX = structfun(@(x) any(strcmp(x, 'cell')),chan);
            nucleusIDX = structfun(@(x) any(strcmp(x, 'nucleus')),chan);
            polIDX = structfun(@(x) any(strcmp(x, 'polymer')),chan);
            %we normalize to sum the channel and then multiply by mean of
            %the two max
            channel.cell    = tmp.(fields{cellIDX}); %normalizing factor
            channel.polymer = tmp.(fields{polIDX});
            
            if all(nucleusIDX==0)
                channel.nucleus = [];
            else
                channel.nucleus = tmp.(fields{nucleusIDX});
               
            end
            
            assert(all(size(channel.cell)==size(channel.polymer)),'Something is wrong with the size of the channels')
            
        end
        
        function check = existChannel(obj)
            path = obj.raw.path;
            filename = [path filesep 'channels.mat'];
            check = isfile(filename);
            
        end
        
        function check = existSegmentation(obj)
            path = obj.raw.path;
            filename = [path filesep 'SegmentedCell.tif'];
            check = isfile(filename);
            
        end
        
        function closestDist  = getDistBetweenMasks(obj,bigMask,smallMask)
            pxSizeXY = obj.info.pxSizeXY;
            pxSizeZ  = obj.info.pxSizeZ;
            %get big contour
            bigContour = imgradient3(bigMask);
            bigContour(bigContour>0) = 1;
            
            %get smallContour
            smallContour = imgradient3(smallMask);
            smallContour(smallContour>0) = 1;
            
            [y,x,z] = ind2sub(size(bigContour),find(bigContour));
            bigCoord = [x*pxSizeXY,y*pxSizeXY,z*pxSizeZ];
            
            [y,x,z] = ind2sub(size(smallContour),find(smallContour));
            smallCoord = [x*pxSizeXY,y*pxSizeXY,z*pxSizeZ];
            
            closestDist = zeros(size(bigContour,1),1);
            
            for i = 1:size(bigCoord,1)
               
                currentCoord = bigCoord(i,:);
                
                closestDist(i) = min(sqrt((currentCoord(1)-smallCoord(:,1)).^2+...
                    (currentCoord(2)-smallCoord(:,2)).^2+(currentCoord(3)-smallCoord(:,3)).^2));
                
                
                
                
                
            end
            
                
            
            
        end
        
        
    end
    
    methods (Static)
        
        function [stat,cornerMat] = getCornerStats(currData)
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

            stat.mean = mean(cornerMat(:));
            stat.median = median(cornerMat(:));
            %get maximum of corner without outliers
            stat.max    = prctile(cornerMat(:),99); 
            
        end
        
        function checkExtension(ext)
            ext = lower(ext);
            extensionList = {'.tif'};
            
            check = contains(extensionList,ext);
            
            assert(sum(check)==1,'Error, unknown extension provided');
            
        end
        
        function [file2Analyze] = getFileInPath(path, ext)
            %Small method to extract the file of a certain extension in a
            %given path
            assert(ischar(path),'The given path should be a char');
            assert(ischar(ext),'The given extension should be a char');
            assert(isfolder(path),'The path given is not a folder')
            
            folderContent = dir(path);
            index2Images  = contains(lower({folderContent.name}),ext,'IgnoreCase',true);
            file2Analyze  = folderContent(index2Images);
            
            assert(~isempty(file2Analyze),['No file of extension ', ext,' found in the input folder']')
        end
        
        function [contour] = getCellContour(gBW)
            contour = cell(1,size(gBW,3));
            for i = 1:size(gBW,3)
                currBW = gBW(:,:,i);
   
                %Get the largest area
                cBWarea = regionprops(currBW,'Area');
                [~,idx2BiggestArea] = max(cell2mat({cBWarea.Area}));
              
                if isempty(idx2BiggestArea)
                else
                    
                    [pContour] = bwboundaries(currBW);
                    for j = 1:length(pContour)
                        contour{i}{j} = pContour{j};
                    end
                end
            end
        end
    end
end

