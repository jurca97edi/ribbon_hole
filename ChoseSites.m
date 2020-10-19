%% ChoseSites
%> @brief function to pick the central sites in the scattering region
%> @param coords An instance of class #Coordinates containing the geometry data
%> @return Returns with an array of logical values to label the indexes to be kept.
    function center_sites = ChoseSites( coords )
        %ribbon_length = (max(coords.y) - min(coords.y));
        center_sites = abs(coords.y - mean(coords.y)) < norm(coords.a)/4 +0.001;
        %unique_x = unique(coords.x)
        %end_sites = unique_x(1:2);
        %uneven_site = (mod(coords.x - end_sites(1),3) == 0) | (mod(coords.x - end_sites(2),3) == 0);
        %plot(coords.x,coords.y,'x')
        %hold on
        %plot(coords.x(center_sites),coords.y(center_sites),'x')
        
        right_min = min(coords.x(coords.x > mean(coords.x)));
        right_max = max(coords.x(coords.x > mean(coords.x)));
        center_right = (right_min + right_max)/2;
        
        left_min = min(coords.x(coords.x < mean(coords.x)));
        left_max = max(coords.x(coords.x < mean(coords.x)));
        center_left = (left_min + left_max)/2;
        
        center_sites = center_sites & ( abs(coords.x - center_right) < 30 | abs(coords.x -center_left) < 30 );
        
        %center_sites = center_sites & uneven_site;
        %center_sites2 = coords.y == mean(coords.y);
        %plot(coords.x(center_sites),coords.y(center_sites),'x')
        
        %plot(coords.x(uneven_site),coords.y(uneven_site),'x')
    end