function [ imgOut ] = radial2XY(rayIn, numQuads)
    rows = 240;
    cols = 320;
    for quadrant = 1:numQuads
        
        quadStart = 1+(quadrant-1)*90;
        quadEnd = quadStart + 90;
        ray = rayIn(quadStart:quadEnd,:);

        % now reconstruct it
        % do first row
        imgOut = [];
        imgOut(1,:) = ray(1,:);

        % do lower triangle
        % below 45 degrees the ray is indexed on columns
        if quadrant > 1
            temp = rows;
            rows = cols;
            cols = temp;
        end
        for i = 2:rows
            for j = i:cols
                if j == i
                    imgOut(i,j) = ray(46,j);        % this is 45 degrees because the index is one above the degrees
                else
                    imgOut(i,j) = radial2XYpoint(ray,i,j);
                end
            end
        end

        %do upper triangle, this does columns from the 45 degree line to the top
        %edge
        % from 45 to 90 degrees the ray is indexed on the row
        for i = 2:cols
            for j = i+1:rows                        % start at i + 1 to avoid redoing the 45 degree line which was done above
                imgOut(j,i) = radial2XYpoint(ray,j,i);
            end
        end

        imgOut = mirrorVert(imgOut);
        if quadrant == 1
            imgOut1 = imgOut;
            %imtool(imgOut1/256)
        elseif quadrant == 2
            imgOut2 = imrotate(imgOut, 90);
            imgOut2 = imgOut2(81:end,:);
            %imtool(imgOut2/256);
        elseif quadrant == 3
            imgOut3 = imrotate(imgOut, 180);
            %imtool(imgOut3/256);    
        elseif quadrant == 4
            imgOut4 = imrotate(imgOut, 270);
            imgOut4 = imgOut4(1:end-80,:);
            %imtool(imgOut4/256);  
        end
    end

    if numQuads == 1
        imgOut = imgOut1(:, 2:end);
    elseif numQuads == 2
        imgOut = imgOut2;
    elseif numQuads == 3
        imgOut = imgOut3;
    elseif numQuads == 4
        imgTop = [imgOut2, imgOut1(:, 2:end)];
        imgBottom = [imgOut3(:, 1:end-1), imgOut4];
        imgOut = [imgTop(1:end-1,:);imgBottom(2:end,:)];
    end
    %imtool(imgOut/max(max(imgOut)))
    imshow(imresize(imgOut, 2/6.5)/max(max(imgOut)))
    %imshow(imgOut/max(max(imgOut)))
    imageSizer
end

