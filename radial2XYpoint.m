function [ val ] = radial2XYpoint( rayArray, row, col )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % lower left corner is (0,0)
    % point is (row, column)
    angleOfRayToPoint = (180/pi) * atan2(row, col);
    % find ray below
    r_low = floor(angleOfRayToPoint)+1;
    r_high = ceil(angleOfRayToPoint)+1;
    
    delta = angleOfRayToPoint + 1 - r_low; 
    
    % find lower value above
    if angleOfRayToPoint <= 45
        u_b = rayArray(r_low, col);
        u_a = rayArray(r_high, col);  
    elseif angleOfRayToPoint < 90
        u_b = rayArray(r_low, row);
        u_a = rayArray(r_high, row); 
    end
    
    val = u_b+(u_a-u_b)* delta;

end

