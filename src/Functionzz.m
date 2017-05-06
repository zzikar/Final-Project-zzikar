function [ F ] = Functionzz( x,y )
%This function evaluates the right hand side at the values that will be required 
%   The right hand side values are generated in a mesh for all the points
%   that are used within the interior. 
%   inputs: x = x values of points within boundaries
%           y = y values of points within boundaries

ax = 0;
ay = 0;
bx = 2*pi;
by = 2*pi;

[X,Y] = meshgrid(x,y);
F =  cos( (0.5*pi)*(( 2 * (X-ax) / (bx -ax) ) + 1)).*sin( pi* (Y-ay) / (by - ay) );
F = -F;
end
