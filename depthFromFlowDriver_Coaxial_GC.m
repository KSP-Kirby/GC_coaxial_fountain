% version 25 Aug 2016 for coaxial-fountain images
% to get ray_out1_f and ray_out1_b run radial2XYdriver on flow_10m.mat
tic

radial2XYdriver

for smoothingFactor = [2]
    for dataFactor = [90]
        for horzNeighborMaskWeight = [5]
            for verticalNeighborMaskWeight = [5]
                %smoothingFactor = 10;    % higher is smoother, zero = no smoothing
                %dataFactor = 90;
                %horzNeighborMaskWeight = 5;
                %verticalNeighborMaskWeight = 2;


                params.columns = 400;
                params.pixelDim = .006;
                deltaX = [3500/30];
                deltaZ = [0];
                params.b = 300;
                results = [];

                minLabel = 60;
                maxLabel = 300;
                f_b = 989.31982*.006;       % from calibration done on 2016/7/16
                f_f = 1307.53356*.006;       % from calibration done on 2016/7/16
                b = 143.251/10;             % this is from my previous coaxial camera, it needs to be measured again, but this should be close
                %b = 115.56/10;              % computed using center flow, 20 mm displacement.
                b = 119.159/10;                 % computed using computeB for 10X10 pixel area around center
                b = 127.76/10;
                f_b = 1005.52555*.006;       % from calibration done on 2016/8/25
                f_f = 1320.1734*.006;       % from calibration done on 2016/8/25
                
                %error.rms_Z should be multiplied by 100 to get percent depth error
                %error.rms_h is the pixel error, don't multiple by 100 and don't report as
                %a percent
                %error.dispErrPercent is the number of pixels that have an error greater
                %than 1

                % possible disparities range from 11 - 20, use labels 10 to 21.
                % 440 pixels (sites), 12 disparities (labels) one epipolar line
                % include the directory gco-v3.0/matlab
                % run GCO_UnitTest() to make sure everything works, otherwise you might
                % need to specify compiler per instructions in ??
                % run GCO_Delete(h) to get rid of object

                %load('flow.mat')
                %load('rayOut1')
                addpath('C:\Users\Richard\Documents\MATLAB\flow-code-matlab')   %this is flow to color

                imageSet = 3;

                usePreviousDataCost = 0;

                startAngle = 1;
                endAngle = 361;
                wf = rayOut1_f(startAngle:endAngle,:);
                wb = rayOut1_b(startAngle:endAngle,:);

                [numRows,numCols] = size(wf);
                optimalLabellingOut = [];

                % Add path where gco-v3.0\matlab is located
                addpath('C:\Users\Richard\Documents\MATLAB\gco-v3.0\matlab')

                numSites = numRows*numCols;

                % labels are in cm
                % min label for the coaxial camera setup is about
                % max label is about
                labels = (minLabel:1:maxLabel);
                numLabels = length(labels);
                h = GCO_Create(numSites,numLabels);

                if usePreviousDataCost == 1
                    load('dataCost.mat')
                else
                    dataCost = zeros(numLabels, numSites);

                    site = 1;

                    h1 = waitbar(0, 'Constructing data cost matrix');
                    for k = 1:numRows
                        for i = 1:numCols
                            for j = 1: numLabels
                                Z = labels(j);
                                %if i + disparity <= numCols
                                    m = (f_b/f_f)*(Z/(Z+b));
                                    scaledPixelLocation = i * m;                                                                % this is a fractional pixel
                                    %wb_intrp = wb(k, floor(scaledPixelLocation)) + ((wb(k,ceil(scaledPixelLocation))-wb(k,floor(scaledPixelLocation))) * (scaledPixelLocation-floor(scaledPixelLocation)));                      % this interpolates between the two interger wb values
                                    sp = round(scaledPixelLocation);
                                    if sp < 1
                                        sp = 1;
                                    end
                                    wb_intrp = wb(k,sp);
                                    dataCost(j, site) = abs(m*wf(k,i) - wb_intrp);
                                %else
                                %    dataCost(j,site) = NaN;
                                %end   
                            end
                            site = site + 1;
                        end
                        waitbar(k/numRows)
                    end
                    close(h1)


                    minDataCost = min(min(dataCost));
                    maxDataCost = max(max(dataCost));

                    %scale data cost between 1 and 1000

                    scaleFactor =50/maxDataCost;
                    dataCost = dataFactor*scaleFactor*dataCost+1;
                end


                GCO_SetDataCost(h,cast(dataCost,'int32'));

                % Smooth cost is number of labels by number of labels
                smoothCost = zeros(numLabels, numLabels);
                for i = 1:numLabels
                    smoothCost(i,:) = abs((1:1:numLabels) - i);
                end

                smoothCost = smoothCost*smoothingFactor;

                GCO_SetSmoothCost(h,cast(smoothCost, 'int32'));

                neighbors = sparse(numSites,numSites);

                colCount = 1;
                for i = 1:numSites-1
                    if i ~=numCols
                        neighbors(i,i+1) = horzNeighborMaskWeight;
                        colCount = colCount + 1;
                    else
                        colCount = 1;
                    end
                end

                % Radial line to radial line neighbors
                for i = 1:numSites-numCols - 1
                    neighbors(i,i+numCols) = verticalNeighborMaskWeight;
                end

                % Tie together first and last ray
                for i = 1:numCols
                    neighbors(i,numSites-numCols + i) = verticalNeighborMaskWeight;
                end

                GCO_SetNeighbors(h,neighbors);

                GCO_Expansion(h);

                [E, D, S] = GCO_ComputeEnergy(h); 

                optimalLabelling = cast(GCO_GetLabeling(h),'double');

                %figure
                %plot(optimalLabelling);

                GCO_Delete(h);

                depthMap = reshape(optimalLabelling,numCols,length(optimalLabelling)/numCols);

                depthMap = depthMap';
                filename=strcat('depthMap_',num2str(smoothingFactor),'_',num2str(dataFactor),'_',num2str(horzNeighborMaskWeight),'_',num2str(verticalNeighborMaskWeight));
                save(filename,'depthMap');     
                
                 
                [ imgOut ] = mirrorHorz(radial2XY(depthMap, 4));
                imtool(imgOut)
            end
        end
    end
    toc
end
