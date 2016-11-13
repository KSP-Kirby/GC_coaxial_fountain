function [ rayOut ] = XY2radial(uv, numQuads)
    step = 1;
    rows = 240;
    cols = 320;

    % this assumes the motion is horizontal
    img = uv(:,:,1);
    [m,n] = size(img);
    
    minDim = min(m,n);
    maxDim = max(m,n);

    dimDiff = maxDim-minDim;

    pad = zeros(dimDiff/2,n);

    img = [pad;img;pad];

    %make sure image has odd number of pixels by choping off a 
    resize = 0;
    if mod(m,2)== 0
        m = m - 1;
        resize = 1;
    end

    if mod(2,2)== 0
        n = n - 1;
        resize = 1;
    end

    if resize == 1
        img = imresize(img,[n,n]);
    end

    indexOut = 1;
    rayOut = [];
    h = waitbar(0,'Extracting rays')
    for quadrant = 1:numQuads
        imgIn = imrotate(img,-90*(quadrant-1));
        ray = [];
        index = 1;
        for i = 0:step:90
            ray(index,:) = extractRay(imgIn,i);
            if ~(i == 0 && indexOut ~= 0)
                rayOut(indexOut,:) = ray(index,:);
                indexOut = indexOut + 1;
            end
            index = index + 1;
            waitbar((i + (quadrant - 1)*90)/360)
        end
    end
    close(h)
    rayOut(361,:) = ray(1,:);   %this wraps the rays around 360 degrees
end

