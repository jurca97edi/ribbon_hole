%% ChoseSites
%> @brief function to pick the central sites in the scattering region
%> @param coords An instance of class #Coordinates containing the geometry data
%> @return Returns with an array of logical values to label the indexes to be kept.
    function center_sites = ChoseSitesOuter( coords )
        %ribbon_length = (max(coords.y) - min(coords.y));
      
        center_sites = abs(coords.y - mean(coords.y)) < norm(coords.a)/4;
    
        figure1 = figure('rend','painters','pos',[10 10 1200 800]);
        
        plot(coords.x,coords.y,'x')
        hold on

        right_min = min(coords.x(coords.x > mean(coords.x) & center_sites));
        right_max = max(coords.x(coords.x > mean(coords.x) & center_sites));
        center_right = (right_min + right_max)/2;
       
        left_min = min(coords.x(coords.x < mean(coords.x) & center_sites));
        left_max = max(coords.x(coords.x < mean(coords.x) & center_sites));
        center_left = (left_min + left_max)/2;

        branch_width = left_max-left_min;

        %original
        %center_sites = abs(coords.y - mean(coords.y)) < 0.01*ribbon_length;
        
        %choosing sites in at the outer edge
        center_sites = center_sites & ( abs(coords.x - right_max) < branch_width/3 | abs(coords.x - left_min) < branch_width/3 );
        
        % choosing sites in the center
        %center_sites = center_sites & ( abs(coords.x - center_right) < branch_width/6 | abs(coords.x - center_left) < branch_width/6 );
        
        %choosing sites in the inner edges
        %center_sites = center_sites & ( abs(coords.x - right_min) < branch_width/3 | abs(coords.x - left_max) < branch_width/3 );
        
%{a
        plot(coords.x(center_sites),coords.y(center_sites),'x')
        
        height = ( max(coords.y) - min(coords.y) )/sqrt(3)+0.5;
                
        print('-dpng', ['ChoseSites_H',num2str(height),'_nu',num2str(sum(center_sites)),'_outer.png'])
        close(figure1);
%}
    end
