function [ b ] = computeB( uv_b, uv_fm )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here



    frontFlow = uv_fm{1}(:,:,1);
    backFlow = uv_b{1}(:,:,1);
    
    delta = 4;
    deltaBack = round(989.31982/1307.53356*delta);
    
    centerFrontFlow = mean(mean(frontFlow(240-delta:240+delta,320-delta:320+delta)));
    depthFront = 20 * 1320.1734 / centerFrontFlow;

    centerBackFlow = mean(mean(backFlow(240-deltaBack:240+deltaBack,320-deltaBack:320+deltaBack)));
    depthBack = 20 * 1005.52555 / centerBackFlow;

    b = depthBack - depthFront;


end

