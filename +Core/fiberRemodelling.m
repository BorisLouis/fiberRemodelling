classdef fiberRemodelling < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        raw
        info
        channels
        results
        validExt = {'.lif','.tif'};
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
            ext  = obj.validExt;
            %let us check that there is no channel data existing
            if ~obj.existChannel
                disp('no channel data found, starting extraction ...')
                
               
                for i = 1:length(ext)
                    currentExt = ext{i};
                    %get all file of appropriate extension in the file
                    fileList = Core.fiberRemodelling.getFileInPath(path,currentExt);
                    
                    if ~isempty(fileList)
                      
                        disp(['A file of extension ', currentExt,' found in the input folder'])
                        break;
                    end
                    
                end
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
         
            nFrames = cellfun(@size,{channel.polymer},'UniformOutput',0);
            nFrames = cellfun(@(x) x(3),nFrames);
            obj.raw.nFrames = nFrames;
            %store the data back into the object
            obj.channels = channel;
            
            disp(['Found ' num2str(length(nFrames)) ' cells in the file']);
            
            
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
            figure
            for j = 1:size(channel,2)
                
                sliceToShow = round(size(channel(j).cell,3)/2);
                frameToShow = round(size(channel(j).cell,4)/2);
                
                for i = 1:nField
                    subplot(size(channel,2),nField,(j-1)*nField+i)
                    currChan = channel(j).(field{i});
                    try
                        imagesc(currChan(:,:,sliceToShow,frameToShow))
                    catch

                    end
                    axis image
                    colormap('hot');
                    title([field{i} ' #' num2str(j)])
                end
            end
            
        end
        
        function [masks] = calc3DMask(obj)
            disp('Starting mask extraction...')
            disp('Performing segmentation on cell channel');
            
            chan = obj.getChannel;
           
            for j = 1:size(chan,2)
                disp(['Analyzing Cell ' num2str(j) ' / ' num2str(size(chan,2))]);
                data1 = chan(j).cell;
                try
                    data2 = chan.nucleus;
                catch 
                    data2 = []; 
                end

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
                    if median(IMs(:)) ==0
                        sig = IMs(IMs>0.1);
                        th = std(sig)/2; 
                        
                    else
                        th = adaptthresh(IMs,0.2,'neigh',[101 101 51],'Fore','bright');
                    end

                    %before
                    %gBW = imbinarize(IMs,'adaptive','Sensitivity',0.2);
                    gBW = imbinarize(IMs,th);
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
                    obj.results.cellContour{j,i} = contour;

                    allMask(:,:,:,i) = gBW;

                end
                
               masks{j} = logical(allMask);
            end
            obj.results.cellMask = masks;
            
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
            for i = 1:size(obj.results.cellMask,2)
                disp(['Rendering Cell ' num2str(i) ' / ' num2str(size(obj.results.cellMask,2))]) 
                data2Render = obj.results.cellMask{i}(:,:,:,idx);
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
                   
                else
                    p2 = patch(iSurface);
                    p2.FaceColor = colorModel;
                    p2.EdgeColor = 'none';
                    view(3);
                    axis tight
                    camlight
                    lighting gouraud
                   
                end
                axis image
                
                title(['Data ' num2str(i)]);
                set(gcf,'Color','w')
                filename = [obj.raw.path filesep 'Figures' filesep 'Cell_0' num2str(i) '.fig'];
                saveas(gcf,filename)
                filename = [obj.raw.path filesep 'Figures' filesep 'Cell_0' num2str(i) '.png'];
                saveas(gcf,filename)
               
                
            end
            
            
        end
        
        function getDensifiedNetwork(obj)
            assert(isfield(obj.channels,'polymer'));
            assert(isfield(obj.results,'cellMask'));
            method = 'tHold';% 'std' tHold
            disp('Extracting densified network...')
            
            disp('Performing segmentation on polymer channel');
            chan = obj.getChannel;
            
            for j = 1:size(chan,2)
                disp(['Analyzing polymer #' num2str(j) ' / ' num2str(size(chan,2))]);
                data = chan(j).polymer;

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

                        for k = 1:size(gBW,3)
                            gBW(:,:,k) = bwareaopen(gBW(:,:,k),50);

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
                    obj.results.polymerContour{j,i} = contour;
         
                %% building 3D mask
            
                currCellMask = obj.results.cellMask{j}(:,:,:,i);
                finalMask = gBW-currCellMask;
                finalMask(finalMask<0)=0;
                allMask(:,:,:,i) = finalMask;
                masks{j} = logical(allMask);
                end
            end
            
            obj.results.polymerMask = masks;
            
        end
                
        function renderCellPolymer3D(obj,idx)
            assert(isfield(obj.results,'cellMask'),'No cell mask found run calc3DMask first');
            assert(isfield(obj.results,'polymerMask'),'No polymer mask found run calc3DMask first');
            
            if nargin <2
                idx = 1;
            end
            for i = 1:size(obj.results.cellMask,2)
                disp(['Rendering Cell ' num2str(i) ' / ' num2str(size(obj.results.cellMask,2))]) 
                data2Render = obj.results.cellMask{i}(:,:,:,idx);
                iSurface = isosurface(data2Render,1/2);

                % smoothing using compiled c code
    %             cellsmoothISurface = rendering3D.smoothpatch(iSurface,0,10);
    %             %comnvert to px
    %             cellsmoothISurface.vertices(:,1) = (cellsmoothISurface.vertices(:,1));
    %             cellsmoothISurface.vertices(:,2) = (cellsmoothISurface.vertices(:,2));
    %             cellsmoothISurface.vertices(:,3) = (cellsmoothISurface.vertices(:,3));
                figure
                hold on
                p2 = patch(iSurface);
                p2.FaceColor = [0,1,0];
                p2.EdgeColor = 'none';

                clear cellsmoothISurface;

                data2Render = obj.results.polymerMask{i}(:,:,:,idx);
                iSurface = isosurface(data2Render,1/2);

                % smoothing using compiled c code
    %             polsmoothISurface = rendering3D.smoothpatch(iSurface,0,1);
    %             %comnvert to px
    %             polsmoothISurface.vertices(:,1) = (polsmoothISurface.vertices(:,1));
    %             polsmoothISurface.vertices(:,2) = (polsmoothISurface.vertices(:,2));
    %             polsmoothISurface.vertices(:,3) = (polsmoothISurface.vertices(:,3));

                %% Displaying network model


                p3 = patch(iSurface);
                p3.FaceColor = [0.7 0 0];
                p3.FaceAlpha = 0.3;
                p3.EdgeColor = 'none';

                view(3);
                axis tight
                camlight
                lighting gouraud
                
                title(['Data ' num2str(i)]);
                
                
                set(gcf,'Color','w')
                filename = [obj.raw.path filesep 'Figures' filesep 'Cell_Pol_0' num2str(i) '.fig'];
                saveas(gcf,filename)
                filename = [obj.raw.path filesep 'Figures' filesep 'Cell_Pol_0' num2str(i) '.png'];
                saveas(gcf,filename)
                
                
            end
        end
        
        function getCellGap(obj)
            assert(isfield(obj.channels,'polymer'));
            assert(isfield(obj.results,'cellMask'));
            
            disp('Extracting polymer gap left by cell...')
            
            disp('Performing segmentation on polymer channel');
            chan = obj.getChannel;
            
            for j = 1:size(chan,2)
                disp(['Analyzing dataset #' num2str(j) ' / ' num2str(size(chan,2))]);
                data = chan(j).polymer;
                cellMask = obj.results.cellMask{j};
                allMask = zeros(size(data));
                for i = 1:size(data,4)
                    currCellMask = cellMask(:,:,:,i);
                    currData = data(:,:,:,i);
                    %filtering
                    pixZ  = 4;
                    zFactor = 2;
                    sigma = [2,2,2*zFactor/pixZ];
                    IMs = imgaussfilt3(currData, sigma);
                    filt = imgradient3(IMs);
                    %threwsholding
                    th = adaptthresh(filt./max(filt(:)),0.8,'neigh',[101 101 51],'Fore','bright');
                    bwfilt = imbinarize(filt./max(filt(:)),th);

                    se = strel('disk',5);
                    bwfilt = imclose(bwfilt,se);

                    bwfilt = ~bwfilt;
                    
                    %remove small object in 2D
                    for k=1:size(bwfilt,3)
                        bwfilt(:,:,k) = bwareaopen(bwfilt(:,:,k),200);
                        
                    end
                    
                  
                    
                    %extract largest
                    stats = regionprops3(bwfilt,'Volume','VoxelIdxList','Extent');
                    [~,idx] = max(stats.Volume.*stats.Extent.^2);
                    idxList = stats.VoxelIdxList{idx,:};

                    bwMask = zeros(size(bwfilt));

                    bwMask(idxList) = 1;

                    se = strel('disk',3);
                    bwMask = imdilate(bwMask,se);
                    finalM = bwMask - currCellMask;
                    finalM(finalM<0)= 0;
                    allMask(:,:,:,i) = finalM;
                    
                    masks{j} = logical(allMask);
%                 
%                 figure
%                 imagesc(bwMask(:,:,100));
%                 
               end
                
            end
            
            obj.results.holeMask = masks;
            
            
        end
        
        function renderAll(obj)
             assert(isfield(obj.results,'polymerMask'),'No polymer mask found run calc3DMask first');
            
            if nargin <2
                idx = 1;
            end
            for i = 1:size(obj.results.cellMask,2)
                disp(['Rendering Cell ' num2str(i) ' / ' num2str(size(obj.results.cellMask,2))]) 
                data2Render = obj.results.cellMask{i}(:,:,:,idx);
                iSurface = isosurface(data2Render,1/2);

                figure
                hold on
                p2 = patch(iSurface);
                p2.FaceColor = [0,1,0];
                p2.EdgeColor = 'none';
                p3.FaceAlpha = 0.7;
              
                data2Render = obj.results.holeMask{i}(:,:,:,idx);
                iSurface = isosurface(data2Render,1/2);
                %% Displaying network model


                p3 = patch(iSurface);
                p3.FaceColor = [0.5 0.5 0.5];
                p3.FaceAlpha = 0.3;
                p3.EdgeColor = 'none';

                
                %% displaying gap
                data2Render = obj.results.polymerMask{i}(:,:,:,idx);
                iSurface = isosurface(data2Render,1/2);

                p4 = patch(iSurface);
                p4.FaceColor = [0.7 0 0];
                p4.FaceAlpha = 0.5;
                p4.EdgeColor = 'none';
                
                
                
                
                
                
                
                view(3);
                axis tight
                camlight
                lighting gouraud
                title(['Data ' num2str(i)]);
                
                set(gcf,'Color','w')
                filename = [obj.raw.path filesep 'Figures' filesep 'Cell_Pol_0' num2str(i) '.fig'];
                saveas(gcf,filename)
                filename = [obj.raw.path filesep 'Figures' filesep 'Cell_Pol_0' num2str(i) '.png'];
                saveas(gcf,filename)
                
            end
        end
        
        function [stats] = calcStats(obj)
            assert(isfield(obj.results,'cellMask'),'No cell mask found run calc3DMask first');
            assert(isfield(obj.results,'polymerMask'),'No polymer mask found run calc3DMask first');
            assert(isfield(obj.results,'holeMask'),'No hole mask found run getCellGap first');
            
            cellMask = obj.results.cellMask;
            polymerMask = obj.results.polymerMask;
            holeMask    = obj.results.holeMask;
            cellChan = obj.channels.cell;
            polChan = obj.channels.polymer;
            
            %memory preallocation
            stats.cellVol = zeros(size(cellMask,2),size(cellMask{1},4));
            stats.polVol  = zeros(size(cellMask,2),size(cellMask{1},4));
            stats.holeVol = zeros(size(cellMask,2),size(cellMask{1},4));
            stats.cellInt = zeros(size(cellMask,2),size(cellMask{1},4));
            stats.polInt = zeros(size(cellMask,2),size(cellMask{1},4));
            stats.holeInt = zeros(size(cellMask,2),size(cellMask{1},4));
            stats.cellCornerInt = zeros(size(cellMask,2),size(cellMask{1},4));
            stats.polCornerInt  = zeros(size(cellMask,2),size(cellMask{1},4));
            
            %remove the cell from the polymer mask
            %CELL IS NOW REMOVED FROM THE POLYMER MASK IN
            %GETDENSIFIEDNETWORK
%             polymerMask=polymerMask-cellMask;
%             polymerMask(polymerMask<0) = 0;
%             polymerMask = logical(polymerMask);
            
            pxSize = obj.info.pxSizeXY/1000;
            pxSizeZ = obj.info.pxSizeZ/1000;
            
            voxelSize = pxSize^2*pxSizeZ;
            for j = 1:size(cellMask,2)
                for i = 1:size(cellMask{1},4)
                    currentCellM = cellMask{j}(:,:,:,i);
                    currentPolM  = polymerMask{j}(:,:,:,i);
                    currentHoleM = holeMask{j}(:,:,:,i);
                    
                    currentCell = cellChan(:,:,:,i);
                    currentPol  = polChan(:,:,:,i); 
                    cellCorner  = Core.fiberRemodelling.getCornerStats(currentCell);
                    polCorner   = Core.fiberRemodelling.getCornerStats(currentPol);

                    %normalize the mean intensity to the pixel volumes and the 
                    %intensity in the corner
                    stats.cellInt(j,i) = mean(currentCell(currentCellM));
                    stats.cellCornerInt(j,i) = cellCorner.median;
                    stats.polInt(j,i)  = mean(currentPol(currentPolM));
                    stats.holeInt(j,i) = mean(currentPol(currentHoleM));
                    stats.polCornerInt(j,i)  = polCorner.median;
                    if isnan(stats.polInt(j,i))
                        stats.polInt(j,i) =0;
                    end
                    stats.cellVol(j,i) = sum(currentCellM(:))*voxelSize;
                    stats.polVol(j,i)  = sum(currentPolM(:))*voxelSize;
                    stats.holeVol(j,i) = sum(currentHoleM(:))*voxelSize;
                end
            end
            stats.voxelSize = voxelSize;
            %store stats in results
            obj.results.stats = stats;
%             hVol = figure(100);
%             hInt = figure(101);
%             for j = 1:size(cellMask,2)
%                  subplot(1,2,1)
%                     plot(stats.cellVol)
%                     ylim([0 1+1.5*max(stats.cellVol)])
%                     axis square
%                     xlabel('Time')
%                     ylabel('Volume (\mum^3)')
%                     title('Cell')
% 
%                 subplot(1,2,2)
%                     plot(stats.polVol)
%                     ylim([0 1+1.2*max(stats.polVol)])
%                     axis square
%                     xlabel('Time')
%                     ylabel('Volume (\mum^3)')
%                     title('Densified polymer')
% 
%             end
%             %Volume plot
%             figure
%                 subplot(1,2,1)
%                     plot(stats.cellVol)
%                     ylim([0 1+1.5*max(stats.cellVol)])
%                     axis square
%                     xlabel('Time')
%                     ylabel('Volume (\mum^3)')
%                     title('Cell')
% 
%                 subplot(1,2,2)
%                     plot(stats.polVol)
%                     ylim([0 1+1.2*max(stats.polVol)])
%                     axis square
%                     xlabel('Time')
%                     ylabel('Volume (\mum^3)')
%                     title('Densified polymer')
% 
%              figure
%                 subplot(1,2,1)
%                     plot(stats.cellInt)
%                     ylim([0 1.5*max(stats.cellInt)])
%                     axis square
%                     xlabel('Time')
%                     ylabel('Intensity(a.u.)')
%                     title('Cell')
% 
%                 subplot(1,2,2)
%                     plot(stats.polInt)
%                     ylim([0 0.001+1.2*max(stats.polInt)])
%                     axis square
%                     xlabel('Time')
%                     ylabel('Intensity(a.u.)')
%                     title('Densified polymer')

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
            AllDist  = cell(size(cellMask,2),size(cellMask{1},4),1);
            
            for j = 1:size(cellMask,2)
                for i = 1:size(cellMask{1},4)
                    currentCellM = cellMask{j}(:,:,:,i);
                    currentPolM  = ~polMask{j}(:,:,:,i);

                    dist = obj.getDistBetweenMasks(currentPolM,currentCellM);

                    AllDist{j,i} = dist;

                end
            end
            
            if size(cellMask{1},4)>1
                %here we assume that the first time frame is 24h before
                %drug so it is supposed to be the "big mask" and need to be
                %inverted
                %the last frame is the small mask
                for j=1:size(polMask,2)
                    polM1 = ~polMask{j}(:,:,:,1);
                    polM2 = polMask{j}(:,:,:,end);
                    dist = obj.getDistBetweenMasks(polM1,polM2);
                    obj.results.polDist{j} = dist;
                end
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
        
        function [channel] = retrieveChannel(obj,fileList,chan)
            
            [folder,file,ext] = fileparts([fileList(1).folder filesep fileList(1).name]);
            fields = fieldnames(chan);
            
            switch (ext)
                case '.tif'
                    [movInfo] = Load.Movie.tif.getInfo([fileList(1).folder filesep fileList(1).name]);
                    [channel]   = obj.extractChannelTIF(movInfo, fields,fileList);
                
                case '.lif'
                    channel  = obj.extractChannelLIF([fileList(1).folder filesep fileList(1).name],chan);
                otherwise
                    error('Unkown extension')
            end
             
            
           
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
        
        function channel = extractChannelTIF(~,movieInfo, fields,fileList)
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
        
        function channel = extractChannelLIF(~,path2File,fields)
            [path,file,ext] = fileparts(path2File);
            
            %extract cell and polymer idx based on user input
            nf = fieldnames(fields);
            for i = 1:length(nf)
                f = fields.(nf{i});
                
                if (strcmpi(f,'cell'))
                    idCell = i;
                elseif (strcmpi(f,'polymer'))
                    idPol = i;
                end
   
            end
            
            data = bfopen(path2File);

            nSeries = size(data,1);
            for i = 1: nSeries
                Misc.multiWaitbar('Channel Extraction',i/nSeries);
                currData = data{i}(:,1);
                nameData = data{i}(1,2);
                
                idx = strfind(nameData{1}, 'C=1/');
                nChan = str2double(nameData{1}(idx+4:end));

                %assign cell and polymer channel based on user input 
                idx2Cell = idCell:nChan:size(currData,1);
                idx2Pol  = idPol:nChan:size(currData,1);
               
                    
                currCellChan = currData(idx2Cell);
                currPolChan = currData(idx2Pol);
                cellChan = zeros(size(currData{1},1),size(currData{1},2),length(currCellChan));
                polChan  = cellChan;
                for j = 1:length(currCellChan)
                    
                    cellChan(:,:,j) = currCellChan{j};
                    polChan(:,:,j)  = currPolChan{j};
                    
                    Misc.multiWaitbar('Frames',j/length(currCellChan));
                end
                   
                channel(i).cell    = cellChan; %normalizing factor
                channel(i).polymer = polChan;
                channel(i).nucleus = [];
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
            
            if isempty(file2Analyze)
                disp(['No file of extension ', ext,' found in the input folder'])
        
            end
        end
        function [contour] = getCellContour(gBW)
            contour = cell(1,size(gBW,3));
            for i = 1:size(gBW,3)
                currBW = gBW(:,:,i);
   
                %Get the largest areaf
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

