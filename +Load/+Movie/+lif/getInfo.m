function [movInfo] = getInfo(path2File)
    [path,file,ext] = fileparts(path2File);
    
    data = bfopen(path2File);
    
    nSeries = size(data,1);
    %get the number of frame
    maxFrame = cellfun(@length,data(:,1));
    %get the exposure time
    expT   = 0;
    %store a few thing for the rest of the code
    isMultiImage = true;
    
    isZStack = true;
    
    Cam  = 0;
    for i = 1:nSeries
        X(i) = size(data{i,1}{1,1},1);
        Y(i) = size(data{i,1}{1,1},2);
    end
        
    %Store info for output
   
    
    movInfo.Width  = X;
    movInfo.Length = Y;
    movInfo.Path   = path;
    movInfo.fullPath = path2File;
    movInfo.isMultiImage = isMultiImage;
    movInfo.isZStack = isZStack;
    movInfo.Cam = Cam;
    movInfo.expT = expT;
    movInfo.maxFrame = maxFrame;
    
    
    
    
    
    


end