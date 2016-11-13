% version 25 Aug 2016
% to get ray_out1_f and ray_out1_b run radial2XYdriver on flow_10m.mat


for smoothingFactor = [1]
    for dataFactor = [90]
        for horzNeighborMaskWeight = [1]
            for verticalNeighborMaskWeight = [1]
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

                minLabel = 40;
                maxLabel = 300;
                f_b = 989.31982*.006;       % from calibration done on 2016/7/16
                f_f = 1307.53356*.006;       % from calibration done on 2016/7/16
                b = 143.251/10;             % this is from my previous coaxial camera, it needs to be measured again, but this should be close

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

                imageSet = 1;

                usePreviousDataCost = 0;

                startAngle = 4;
                endAngle = 4;
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
                                    if scaledPixelLocation < 1
                                        wb_intrp = wb(k,1);
                                    else 
                                        wb_intrp = wb(k, floor(scaledPixelLocation)) + ((wb(k,ceil(scaledPixelLocation))-wb(k,floor(scaledPixelLocation))) * (scaledPixelLocation-floor(scaledPixelLocation)));                      % this interpolates between the two interger wb values
                                    end
                                    %wb_intrp = wb(k,round(scaledPixelLocation));
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
                % this has to be commented out for a single ray
%                 for i = 1:numCols
%                     neighbors(i,numSites-numCols + i) = verticalNeighborMaskWeight;
%                 end

                GCO_SetNeighbors(h,neighbors);

                GCO_Expansion(h);

                [E, D, S] = GCO_ComputeEnergy(h); 

                optimalLabelling = cast(GCO_GetLabeling(h),'double');

                %figure
                %plot(optimalLabelling);

                GCO_Delete(h);
                
                % figure
% plot(optimalLabelling);


% title('Optimal Labeling vs. actual disparity')
% xlabel('Pixels')
% ylabel('Disparity')
% legend('Optimal Labeling', 'Actual Disparity')

                shiftWb = [];
                for i = 1:numSites
                    Z = optimalLabelling(i)+minLabel-1;
                    m = (f_b/f_f)*(Z/(Z+b));
                    shiftWb(1,i) = i/m;
                    if round(i*m) <= numSites
                        shiftWb(2,i) = wb(round(i*m))/m;
                    end
                end

                for i = 1:length(shiftWb)
                    if shiftWb(2,i) == 0
                        shiftWb(2,i) = nan;
                    end
                end

                figure
                plot(wf)
                hold all
                plot(wb)
                plot(shiftWb(2,:))
                
                dispError = wf-shiftWb(2,:);
                rmsError = sqrt(mean((dispError(3:320)).*conj(dispError(3:320))))
                
                title('Flow along epipolar line adjusted by estimated depth')
                legend('Flow1','Flow2', '(1/m)*Flow1(x+h)')
                xlabel('Pixel')
                ylabel('Flow (pixels/frame)')
                
                dispError = shiftWb(2,3:end)-wf(3:end);
                rmsError = sqrt(mean((dispError).*conj(dispError)))
                
                figure
                plot(wf)
                hold all
                plot(wb) 
                
                % evaluate at a pixel
                frontPixel = 149;
                plot(frontPixel,wf(frontPixel),'g*')
                x_b = frontPixel*.006;
                for i = 1:length(labels)
                    z = labels(i);
                    m = (f_b/f_f)*(z/(z+b));
                    scaledPixelLocation = frontPixel * m;                                                                % this is a fractional pixel
                    %wb_intrp = wb(k, floor(scaledPixelLocation)) + ((wb(k,ceil(scaledPixelLocation))-wb(k,floor(scaledPixelLocation))) * (scaledPixelLocation-floor(scaledPixelLocation)));                      % this interpolates between the two interger wb values
                    wb_intrp = wb(round(scaledPixelLocation));
                    back_flow(i) = (m)*wf(frontPixel);
                    pix(i) = scaledPixelLocation;
                end
                plot(pix,back_flow)

                %axis([180,405,5,35])
                title('Energy minimization search path')
                xlabel('Pixels')
                ylabel('Flow (pixels/frame)')
                legend('Flow 1', 'Flow 2',strcat('Eval Point:',num2str(frontPixel)), 'Eval Energy')
                
                figure
                plot(wf)
                hold all
                plot(wb) 
                plot(frontPixel,wf(frontPixel),'g*')
                z = optimalLabelling(frontPixel)+minLabel-1;
                m = (f_b/f_f)*(z/(z+b));
                scaledPixelLocation = frontPixel * m;                                                                % this is a fractional pixel
                wb_intrp = wb(k, floor(scaledPixelLocation)) + ((wb(k,ceil(scaledPixelLocation))-wb(k,floor(scaledPixelLocation))) * (scaledPixelLocation-floor(scaledPixelLocation)));                      % this interpolates between the two interger wb values
                %wb_intrp = wb(round(scaledPixelLocation));
                plot(scaledPixelLocation,wb_intrp,'r*')
                
                title('Flow along a radial epipolar line')
                xlabel('Pixels')
                ylabel('Flow (pixels/frame)')
                legend('Front Flow', 'Back Flow')
                
                figure
                plot(labels, dataCost(:,frontPixel))
                title('Typical Energy Response to change in Z')
                xlabel('Z Estimate in cm')
                ylabel('Energy')
                
                

% % this filles in gaps where optimal labelling jumps by using the previous
% % optimal label.  
% for i = 2:numSites
%     if shiftW0Left(1,i) == 0;
%         shiftW0Left(1,i) = shiftW0Left(1,i-1);
%     end
% end

%figure
%     plot(shiftW0Left)
%     title('Flow matching with Graph Cuts')
%     xlabel('Pixels')
%     ylabel('Flow in pixels')
%     legend('Flow from IR camera', 'Flow from RGB Camera', 'IR flow shifted by GC disparity')
%axis([0,640,11,15])

%rmsError = sqrt(mean((dispError).*conj(dispError)))
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                

                depthMap = reshape(optimalLabelling,numCols,length(optimalLabelling)/numCols);

                depthMap = depthMap';
                filename=strcat('depthMap_',num2str(smoothingFactor),'_',num2str(dataFactor),'_',num2str(horzNeighborMaskWeight),'_',num2str(verticalNeighborMaskWeight));
                save(filename,'depthMap');   
                
            end
        end
    end
end
