IMT = 1;
for smoothingFactor = [6]
    for dataFactor = [80, 85, 90, 95, 100, 105, 110]
        for horzNeighborMaskWeight = [3, 4, 5, 6]
            for verticalNeighborMaskWeight = [1, 2, 3, 4]
                %smoothingFactor = 10;    % higher is smoother, zero = no smoothing
                
                filename=strcat('depthMap_',num2str(smoothingFactor),'_',num2str(dataFactor),'_',num2str(horzNeighborMaskWeight),'_',num2str(verticalNeighborMaskWeight))
                load(filename); 
                [ imgOut ] = radial2XY(depthMap, 4);
                disp(strcat('IMT#:',num2str(IMT),'; Scale: 0.5;SF:',num2str(smoothingFactor),'; DF:',num2str(dataFactor),'; HNMW:',num2str(horzNeighborMaskWeight),'; VNMW:',num2str(verticalNeighborMaskWeight)));
                IMT = IMT+1;
                
            end
        end
    end
end