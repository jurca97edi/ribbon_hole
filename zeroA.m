% Function describing zero vector potential  -  Based on EQuUs v4.8
%    Copyright (C) 2016 Peter Rakyta, Ph.D.
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
%% zeroA
%> @brief Vector potential in the Landau gauge parallel to the y direction.
%> @param x X coordinates of the sites.
%> @param y Y coordinates of the sites.
%> @return Returns with N x 2 array containing the x and y components of the vector potential. N is the number of sites.
    function ret = zeroA( x,y )
        ret = zeros(size(x,1),2);
    end