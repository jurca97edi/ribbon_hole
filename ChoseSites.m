%% ChoseSites
%> @brief function to pick the central sites in the scattering region
%> @param coords An instance of class #Coordinates containing the geometry data
%> @return Returns with an array of logical values to label the indexes to be kept.
    function center_sites = ChoseSites( coords )
        ribbon_length = (max(coords.y) - min(coords.y));
        center_sites = abs(coords.y - mean(coords.y)) < 0.1*ribbon_length;
    end