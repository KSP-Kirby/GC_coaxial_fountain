function [ rmsErr ] = computErrorGC( frontFlow, depthMap, minLabel, labels )
% Call: computErrorGC(rayOut1_f, depthMap, minLabel, labels )
% Run depthFromFlowDriver_Coaxial_GC first
% August 25 2016 for Coaxaial fountain images

left = 315;
right = 325;
top = 235;
bottom = 245;

depthMapXY = mirrorHorz(radial2XY(depthMap, 4));
frontFlow = mirrorHorz(radial2XY(frontFlow, 4));

imshow(depthMapXY/length(labels));


hold on
plot(left:right,top:top,'r')
plot(right:right,top:bottom,'r')
plot(left:right,bottom:bottom,'r')
plot(left:left,top:bottom,'r')
hold off;

stepSize = 10;

depth = [];
flow = [];
Vel = [];
if stepSize == 10
    index = 1;
    for i = left:right
        for j = top:bottom
            depth(index) = ((depthMapXY(j,i) - 1) + minLabel)*10;
            flow(index) = frontFlow(j,i);
            Vel(index) = flow(index).*depth(index)/1320.1734;
            index = index + 1;
        end
    end
else
    index = 1;
    for i = left:right
        for j = top:bottom
            depth(index) = ((depthMapXY(j,i) - 1)*10 + minLabel)*10;
            flow(index) = frontFlow(j,i);
            Vel(index) = flow(index).*depth(index)/1320.1734;
            index = index + 1;
        end
    end
end
    

vMean = 20;
rmsErr = sqrt(mean((Vel-vMean).*conj(Vel-vMean)));
disp(strcat('RMS Error:',num2str((rmsErr/20)*100),'%'))

vMean = mean(Vel);
rmsErr = sqrt(mean((Vel-vMean).*conj(Vel-vMean)));
disp(strcat('RMS Error:',num2str((rmsErr/20)*100),'%'))

end

