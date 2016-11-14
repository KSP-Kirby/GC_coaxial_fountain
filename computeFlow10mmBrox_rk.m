% computes flow using brox
% 29 Sep 2015 data, 10mm increments, images: [1, 8, 11, 13, 15]
% method = {'classic+nl-fast', 'classic+nl', 'classic+nl-full', 'hs-brightness', 'hs'};
% 'hs-brightness' and 'hs' don't work, may need to make image square

addpath 'C:\Users\Richard\Documents\MATLAB\pami2010Matlab';
addpath 'C:\Users\Richard\Documents\MATLAB\flow-code-matlab';

h = waitbar(0,'Processing Data')
images = [1, 2, 3, 4, 5];
%images = [4, 5, 6];

%b = 55.84;
fl_f = [ 1285.77473   1284.27991 ] * .006;  % for front camera use y
fl_b = [ 973.17464   972.30769  ] * .006;   % for back camera use x

%for jj = 1:length(method)

framesBetweenFlowComp = 4;

    for j = 1:length(images)-framesBetweenFlowComp
        i = images(j);
        i1 = images(j+framesBetweenFlowComp);
        if i < 10
            fc=imread(strcat('fc0',num2str(i),'_rect.tif'));
            bc=imread(strcat('bc0',num2str(i),'_rect.tif'));
            titlef1 = strcat('fc0',num2str(i));
            titleb1 = strcat('bc0',num2str(i));
        else
            fc=imread(strcat('fc',num2str(i),'_rect.tif'));
            bc=imread(strcat('bc',num2str(i),'_rect.tif'));
            titlef1 = strcat('fc',num2str(i));
            titleb1 = strcat('bc',num2str(i));
        end

        if i1 < 10
            fc2=imread(strcat('fc0',num2str(i1),'_rect.tif'));
            bc2=imread(strcat('bc0',num2str(i1),'_rect.tif'));
            titlef2 = strcat('fc0',num2str(i1));
            titleb2 = strcat('bc0',num2str(i1));
        else
            fc2=imread(strcat('fc',num2str(i1),'_rect.tif'));
            bc2=imread(strcat('bc',num2str(i1),'_rect.tif'));
            titlef2 = strcat('fc',num2str(i1));
            titleb2 = strcat('bc',num2str(i1));
        end

        plotTitle = {strcat('From 2 July 2015 Images:', titlef1,'-',titlef2,'<<->>', titleb1, '-', titleb2);...
            strcat('method: Brox')};
        
        [m, n, p] = size(fc);
        if p == 1
            fc = convert1Dto3D(fc);
            bc = convert1Dto3D(bc);
            fc2 = convert1Dto3D(fc2);
            bc2 = convert1Dto3D(bc2);
        end
        
        trimX = 80;
        trimY = 80;
        horzLine = (480 - 2*trimX)/2;

        fc = mirrorHorz(fc);
        fc2 = mirrorHorz(fc2);
        
        %uv_fm{i} = estimate_flow_interface(fc(100:end-100,100:end-100,:), fc2(100:end-100,100:end-100,:),method{jj});
        %uv_fm{j} = mex_LDOF(double(fc(trimX:end-trimX,trimY:end-trimY,:)), double(fc2(trimX:end-trimX,trimY:end-trimY,:)));
        uv_fm{j} = mex_LDOF(double(fc), double(fc2));
        waitbar((j-.5)/length(images))

        uvi_f{j} = uint8(flowToColor(uv_fm{j}));

        %uv_b{i} = estimate_flow_interface(bc(100:end-100,100:end-100,:), bc2(100:end-100,100:end-100,:),method{jj});
        %uv_b{j} = mex_LDOF(double(bc(trimX:end-trimX,trimY:end-trimY,:)), double(bc2(trimX:end-trimX,trimY:end-trimY,:)));
        uv_b{j} = mex_LDOF(double(bc), double(bc2));
        uvi_b{j} = uint8(flowToColor(uv_b{j}));

        curFigure = figure;
        plot(uv_fm{j}(horzLine,:,1),'LineWidth',2)
        hold all
        plot(uv_b{j}(horzLine,:,1),'LineWidth',2)
        title(plotTitle);
        legend('Front Camera', 'Back Camera')
        ylabel('pixels')
        xlabel('Flow in pixels along center horizontal line')
        %axis([0,450,0,20])
        hold off
        filename = strcat('brox',num2str(images(j)),'.jpg');
        saveas(curFigure,filename)
        filename = strcat('brox',num2str(images(j)),'.fig');
        saveas(curFigure,filename)

        waitbar((j)/length(images))
    end
    
%end

close(h)

save('flow','uv_b','uv_fm','uvi_b','uvi_f')
clear 'uv_b'
clear 'uv_fm'
clear 'uvi_b'
clear 'uvi_f'
close all
