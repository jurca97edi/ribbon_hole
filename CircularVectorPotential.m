% Function describing circular vector potential  -  Based on EQuUs v4.8
%    Copyright (C) 2018 Borbala Farkas
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.
%
%> @brief Function describing circular vector potential
%> @param x X coordinates of the sites.
%> @param y Y coordinates of the sites.
%> @param phi The piercing flux
%> @param cCircle An instance of sructure circle describing a circle shaped area positioned in a two-dimensional space.
%> @return Returns with N x 2 array containing the x and y components of the vector potential. N is the number of sites.
function ret = CircularVectorPotential( x , y, phi, cCircle)

    % determine the x,y coordinates relative to the center of the circle
    x = x - cCircle.center.x;
    y = y - cCircle.center.y;

    % preallocating array
    ret = zeros( length(x),2);
    
    % calculating the radial distance
    r_array = sqrt(x.^2+y.^2);
    
     
       
    %add a homogene magnetic field
    ret(:,1) = -phi/(2*pi)*y./(r_array.^2);
    ret(:,2) = phi/(2*pi)*x./(r_array.^2);
    
    
end