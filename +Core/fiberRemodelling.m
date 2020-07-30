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
                
                pxSizeXY = obj.info.pxSizeXY;
                pxSizeZ  = obj.info.pxSizeZ;
                                
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
        
        function showChannel(obj)
            channel = obj.getChannel;
            field = fieldnames(channel);
            
            nField = length(field);
            frameToShow = round(size(channel.cell,3)/2);
            figure
            for i = 1:nField
                subplot(1,nField,i)
                currChan = channel.(field{i});
                try
                    imagesc(currChan(:,:,frameToShow))
                catch
                    
                end
                axis image
                colormap('hot');
                title(field{i})
            end
            
        end
        
        function [mask3D] = calc3DMask(obj)
            disp('Starting mask extraction...')
            fileName = [obj.raw.path filesep 'SegmentedCell.tif'];
            
            disp('Performing segmentation on cell channel');
            chan = obj.getChannel;
            data1 = chan.cell;
            data2 = chan.nucleus;
            
            if isempty(data2)
                data = data1;
            else
                data = data1+data2;
            end
            
            % Gaussian filtering
            S = 2;
            % size of pixel in z vs x/y
            pixZ  = 4;
            zFactor = 2;
            sigma = [S,S,S*zFactor/pixZ];
            IMs = imgaussfilt3(data, sigma);
            disp('DONE with filtering ------------')

            gBW = imbinarize(IMs,'adaptive','Sensitivity',0.2);
            %Segmentation
          %  gBW = imbinarize(IMs);
            se = strel('disk',10);
            gBW = imopen(gBW,se);
            se = strel('disk',5);
            gBW = imclose(gBW,se);

            %storing temporary results             
            dataStorage.BinaryTiff(fileName,gBW);
          
            disp('Extracting Contour')
            %here we obtain the cell contour
            contour = obj.getCellContour(gBW);
            obj.results.contour = contour;
            
            %% Extract contour
            fContour = [];
            for z = 1 :size(IMs,3)
                if ~isempty(contour{z})
                    for i = 1:length(contour{z})
                    
                        fContour = [fContour ; contour{z}{i}(:,1) contour{z}{i}(:,2) ones(length(contour{z}{i}(:,1)),1)*z];
                    
                    end
                end
            end   
            %% building 3D mask
            mask3D = zeros(size(IMs));
            for i = 1:size(IMs,3)
                idx = find(fContour(:,3)==i);
                mask3D(:,:,i) = poly2mask(fContour(idx,2),fContour(idx,1),size(IMs,1),size(IMs,2));
            end
            
            %% Cleaning mask
            
            test = bwlabeln(mask3D);
            data = regionprops3(test,'Volume','VoxelList','VoxelIdxList');
            %Delete mask region that are extended only on 1 z slice
            for i = 1:height(data)
               currentIdx = data.VoxelList{i,:};
               if length(unique(currentIdx(:,3)))>1
                   
               else
                   idx2Delete = data.VoxelIdxList{i,:};
                   mask3D(idx2Delete) = 0;
               end
                
            end
            
            %Keep only the main mask region (delete all smaller regions)
            maskBW = bwlabeln(mask3D);
            
            maskProps = regionprops3(test,'vol','VoxelIDXList');
            
            if height(maskProps)>1
                [val,idx] = max(maskProps.Volume);
                maskProps(idx,:) = [];
                
                for i = 1 : height(maskProps)
                    mask3D(maskProps.VoxelIdxList{i}) = 0;
                end
                
            end
            maskBW = bwlabeln(mask3D);
            assert(length(unique(maskBW))==2,'More than one region in the mask, something is strange')
                        
            obj.results.mask = logical(mask3D);       
            
        end
        
        function plotCellContour(obj)
            
            contour = obj.results.contour;
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
        
        function renderCell3D(obj)
            %compile c code for smoothing
            mex +rendering3D\smoothpatch_curvature_double.c -v
            mex +rendering3D\smoothpatch_inversedistance_double.c -v
            mex +rendering3D\vertex_neighbours_double.c -v
            
            
            data2Render = obj.results.mask;
            iSurface = isosurface(data2Render,1/2);
            
            % smoothing using compiled c code
            smoothISurface = rendering3D.smoothpatch(iSurface,0,10);
            %comnvert to px
            smoothISurface.vertices(:,1) = (smoothISurface.vertices(:,1));
            smoothISurface.vertices(:,2) = (smoothISurface.vertices(:,2));
            smoothISurface.vertices(:,3) = (smoothISurface.vertices(:,3));

            %% Displaying network model
            %z-coloring
            colorModel = smoothISurface.vertices(:,3)/max(smoothISurface.vertices(:,3));
            zColor = true;
           
            %Plot the network with Z coloring or unique color depending on the user
            %input
            figure
            if zColor
                p = patch('Faces',smoothISurface.faces,'Vertices',smoothISurface.vertices,'FaceVertexCData',colorModel,'FaceColor','interp');
                colormap('jet')
                p.EdgeColor = 'none';
                daspect([2 2 1])
                view(3);
                axis tight
                camlight
                lighting gouraud
                title('Z-coloring')
            else
                p2 = patch(smoothISurface);
                p2.FaceColor = colorModel;
                p2.EdgeColor = 'none';
                view(3);
                axis tight
                camlight
                lighting gouraud
                title('unicolor');
            end

            
            
        end
        
        function intensityDistrib(obj,weight,step)
            
            assert(~isempty(obj.channels),'Data not found, please loadData() before using this function')
            assert(~isempty(obj.results.mask),'mask not found, please calc3DMask before using this function')
            data = obj.channels.polymer;
            mask = obj.results.mask;
            
            EDM = DistMap.calcWeightedDistMap(mask,weight);
            
            edgeMin = min(EDM(:));
            edgeMax = max(EDM(:));
            
            binEdges = edgeMin:step:edgeMax;
            
            intRes = table(zeros(length(binEdges)-1,1),zeros(length(binEdges)-1,1),...
                'VariableNames',{'Distance','Intensity'});
            
            %Calculate the intensity vs distance from the cell
            for i = 1: length(binEdges)-1
               
                idx = and(EDM>=binEdges(i), EDM<binEdges(i+1));
                
                currDistance = (binEdges(i)-binEdges(1) +binEdges(i+1))/2;
                currentIntensity = mean(data(idx));
                
                intRes.Distance(i) = currDistance;
                intRes.Intensity(i) = currentIntensity;
                
            end
            
            intRes.Mean = ones(length(intRes.Distance),1)*mean(data(EDM>0));
            intRes.Median = ones(length(intRes.Distance),1)*median(data(EDM>0));
            
            obj.results.intRes = intRes;
            %save the intensity curve
            filename = [obj.raw.path filesep 'IntensityResults.mat'];
            save(filename,'intRes');
            
            figure
            plot(intRes.Distance,intRes.Intensity);
            xlabel('Distance (nm)')
            ylabel('Average intensity per pixel')
            axis square
            box on
            hold on
            
            plot(intRes.Distance,ones(1,length(intRes.Distance))*mean(data(EDM>0)))
            plot(intRes.Distance,ones(1,length(intRes.Distance))*median(data(EDM>0)))
            legend({'Experiment', 'Mean', 'Median'})
            
            disp('=====> DONE <=====');
        end
        
        
    end
    methods (Access = private)
        
        function [channel] = retrieveChannel(~,fileList,chan)
            
            fields = fieldnames(chan);
            %check movie size
            [movieInfo] = Load.Movie.tif.getinfo([fileList(1).folder filesep fileList(1).name]);
            for i = 1: length(fields)
                Misc.multiWaitbar('Channel Extraction',i/length(fields))
                index2Channel  = contains(lower({fileList.name}),fields{i},'IgnoreCase',true);
                currChanList  = fileList(index2Channel);
                
                currChan = zeros(movieInfo.Length,movieInfo.Width,length(currChanList));
                for j = 1:length(currChanList)
                    currPath = [currChanList(j).folder, filesep, currChanList(j).name];
                    currChan(:,:,j) = Load.Movie.tif.getframes(currPath,1);
                    
                    Misc.multiWaitbar('Frames',j/length(currChanList));
                    
                end
                tmp.(fields{i}) = currChan;
                
            end
            multiWaitbar( 'CloseAll' );
            cellIDX = structfun(@(x) any(strcmp(x, 'cell')),chan);
            nucleusIDX = structfun(@(x) any(strcmp(x, 'nucleus')),chan);
            polIDX = structfun(@(x) any(strcmp(x, 'polymer')),chan);
            %we normalize to sum the channel and then multiply by mean of
            %the two max
            channel.cell    = tmp.(fields{cellIDX}); %normalizing factor
            channel.polymer = tmp.(fields{polIDX});
            if ~isempty(tmp.(fields{nucleusIDX}))
                channel.nucleus = tmp.(fields{nucleusIDX});
            else
                channel.nucleus = [];
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
        
        
    end
    
    methods (Static)
        
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
                    %kill all the other area found
                    [pContour] = bwboundaries(currBW);
                    contour{i} = {pContour{idx2BiggestArea}};

                    idx = find(cell2mat({cBWarea.Area})>5000);
                    idx(idx==idx2BiggestArea) = [];

                    for j = 1:length(idx)
                        contour{i}{j+1} = pContour{idx(j)};
                    end
                end
            end
        end
    end
end

