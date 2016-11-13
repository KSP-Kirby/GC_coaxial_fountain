% 25 August 2016
% creates an array of radial lines as an input to
% depthFromFlowDriver_coaxial_GC

load('flow_rect20.mat');

imageSet = 1;
numQuads = 4;

rayOut1_f = XY2radial(uv_fm{imageSet},numQuads);
imgOut_f = radial2XY(rayOut1_f,numQuads);

rayOut1_b = XY2radial(uv_b{imageSet},numQuads);
imgOut_b = radial2XY(rayOut1_b,numQuads);