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
%> @brief Function describing constant vector potential
%> @param x X coordinates of the sites.
%> @return Returns with N x 2 array containing the x and y components of the vector potential. N is the number of sites.
function ret = ConstantVectorPotential( x , y, flux )

    % preallocating array
    ret = flux/(2*pi);
    
end