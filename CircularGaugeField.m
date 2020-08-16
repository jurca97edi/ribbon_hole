% Function describing circular gauge field  -  Based on EQuUs v4.8
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
%> @brief Function describing circular gauge field 
%> @param x X coordinates of the sites.
%> @param y Y coordinates of the sites.
%> @param eta_B Dimensionless strength of the magnetic field.
%> @param cCircle An instance of sructure circle describing a circle shaped area positioned in a two-dimensional space.
%> @return Returns with N x 1 array containing the scalar potential field. N is the number of sites.
function ret = CircularGaugeField(x,y, flux, cCircle)


% determine the x,y coordinates relative to the center of the circle
x = x - cCircle.center.x;
y = y - cCircle.center.y;

% determine the angle of the [x;y] vector
z = x + 1i*y;
phi = angle(z);
   
   %this term transform the homogenious B from the circular gauge to the Landau gauge
   ret = -flux/(2*pi)*phi;
   
end
            
   