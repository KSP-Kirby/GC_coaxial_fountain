function [ ray ] = extractRay(img, rayAngle)
% extractRay creates a 1d array of image values (gray or flow) from the input 
% image img at rayAngle in degrees

displayOn = 0;
    
    [m, n, p] = size(img);
    center = [floor(m/2)+1,floor(n/2)+1]; 
    
    if rayAngle >= 0 && rayAngle <= 90
        img = img(1:center(1), center(2):end);
    elseif rayAngle > 90 && rayAngle <= 180
        img = img(1:center(1), 1:center(2));
        img = img';
        img = mirrorHorz(img);
        rayAngle = rayAngle - 90;
    elseif rayAngle > 180 && rayAngle <= 270
        img = img(center(1):end, 1:center(2));
        %imtool(img)  
        mirrorHorz(img);
        mirrorVert(img);
    else
        img = img(center(1):end, center(2):end);
        %imtool(img)          
    end
        
    if rayAngle == 0  
        % [1,1] is upper left
        % [1,m] is lower left
        % [n,1] is upper right
        % [n,m] is lower right
        % disp('Horizontal line to right')
        ray = img(end,1:end,1);
        
    elseif rayAngle < 45
        % disp('Lower triangle')
        
        tanRatio = tan(rayAngle*pi/180);
        
        if displayOn
            imshow(img)
            hold on
        end
        
        rows = size(img,1); 
        ray(1) = img(end,1);
        for col = 2:size(img,2)
            r = rows - col*tanRatio;
            a = ceil(r);
            b = floor(r);
            if b < 1
                b = 1;
            end
            if a < 1
                a = 1;
            end
            if a > rows
                a = rows;
            end
            u_b = img(b,col);
            u_a = img(a,col);
            ray(col) = u_b+(u_a-u_b)*(r-b);
            if displayOn == 1
                plot(col,r,'*r')
            end
        end
    elseif rayAngle == 45
        % disp('Quadrant 1: 45 degree angle')
        
        rows = size(img,1);
        if displayOn == 1
            imshow(img(:,:,1))
            hold on
        end
        
        
        minDim = min(size(img,1),size(img,2));
        for r = 1:minDim
            r_i = rows - r + 1;
            ray(r) = img(r_i,r);
            if displayOn == 1
                plot(r,r_i,'*r')
            end
        end      
    elseif rayAngle < 90
        % disp('Quadrant 1 upper triangle')
        
        tanRatio = tan(rayAngle*pi/180);
        rows = size(img,1); 
        
        if displayOn == 1
            imshow(img(:,:,1))
            hold on
        end

        ray(1) = img(end,1);
        for r = 2:rows
            c = r/tanRatio;
            r_i = rows - r + 1;
            a = ceil(c);
            b = floor(c);
            if b < 1
                b = 1;
                a = 2;
            end 
            u_b = img(r_i,b);
            u_a = img(r_i,a);
            ray(r) = u_b+(u_a-u_b)*(c-a);
            if displayOn == 1
                plot(c+1,r_i,'*r')
            end
        end
        
    elseif rayAngle == 90
        %disp('90 degree straight up')
        if displayOn == 1
            imshow(img(:,:,1))
            hold on
        end
        rows = size(img,1);
        for r = 1:rows
            r_i = rows - r + 1;
            ray(r) = img(r_i,1); 
            if displayOn == 1
                plot(1,r_i,'*r')
            end
        end     
    end
    
%     if saveOn
%         figure
%         plot(frontFlow,'LineWidth',5)
%         hold all
%         plot(backFlow,'LineWidth',5)
%         legend('Front Flow','Back Flow')
%         title(strcat('Flow in u direction along radial line with angle:', num2str(rayAngle),' degrees'))
%         xlabel('pixels from center')
%         ylabel('flow in u direction in pixels')
%         axis([1,size(uvF,2),6,20])
%         openFigures = findall(0,'type','figure');
%         filename = strcat('flowPlot',num2str(uint16(rayAngle)),'.jpg');
%         saveas(openFigures(1),filename)
%         close all
%     else
%         filename = '';
%     end


end

