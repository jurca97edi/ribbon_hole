%%    Eotvos Quantum Transport Utilities - Ribbon_Keldysh
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
%
%> @addtogroup utilities Utilities
%> @{
%> @file Ribbon_Keldysh.m
%> @brief  A class representing a two-terminal structure defined on a preprogrammed lattices for steady state non-equilibrium calculations.
%> @image html Ribbon_structure.jpg
%> @image latex Ribbon_structure.jpg
%> @} 
%> @brief A class representing a two-terminal structure defined on a preprogrammed lattices for steady state non-equilibrium calculations.
%> @Available
%> EQuUs v4.9 or later
%> <tr class="heading"><td colspan="2"><h2 class="groupheader"><a name="avail"></a>Structure described by the class</h2></td></tr>
%> @image html Ribbon_structure.jpg
%> @image latex Ribbon_structure.jpg
%> The drawing represents a two-terminal structure made of two leads, of a scattering region and of two interface regions between the leads and scattering center.
%> Each rectangle describes a unit cell including singular and non-singular sites. 
%> The scattering center is also described by a set of identical unit cells, but arbitrary potential can be used.
%> Arrows indicate the hopping direction stored by the attributes in the corresponding classes (see attributes #CreateLeadHamiltonians.H1 and #InterfaceRegion.Hcoupling for details).
%> The orientation of the lead is +1 if the lead is terminated by the interface region in the positive direction, and -1 if the lead is terminated by the interface region in the negative direction.
%> (see attribute #CreateLeadHamiltonians.Lead_Orientation for details)
%%
classdef Ribbon_hole < Ribbon

    properties (Access = protected)
        %> an instance of structure circle describing the inner hole
        cCircle_in
        %> an instance of structure circle describing the outher boundaries of the scattering region
        cCircle_out
        %> The width of the ribbon in the middle
        middle_width
        %> The width of the lead
        lead_width
    end       
    
    properties (Access = public)    
    end
    


methods ( Access = public )
    
    
%% Contructor of the class
%> @brief Constructor of the class.
%> @param varargin Cell array of optional parameters. For details see #InputParsing.
%> @return An instance of the class
function obj = Ribbon_hole( varargin )
    obj = obj@Ribbon();
    
    if strcmpi( class(obj), 'Ribbon_hole') 
        % processing the optional parameters
        obj.InputParsing(varargin{:})                 

        
        obj.display(['EQuUs:Utils:', class(obj), ':Ribbon_Keldysh: Creating a Ribbon_Keldysh object'])
        
        % create the shape of the scattering region
        obj.createShape();
        
        % set the Fermi level
        obj.setFermiEnergy();
        
        % exporting calculation parameters into an XML format
        createOutput( obj.filenameOut, obj.Opt, obj.param );
        
        % create class intances and initializing class attributes
        obj.CreateHandles();
        
        % create Hamiltonians and coordinates of the unit cells 
        obj.CreateRibbon();
    end        
    


end


%% CreateScatter
%> @brief Initializes class #CreateHamiltonians for storing and manipulate the Hamiltonian of the the scattering region. The created object is stored in attribute #CreateH.
    function Scatter_UC = CreateScatter( obj ) 
        Scatter_UC = obj.Scatter_UC.CreateClone('empty', true);
        
        % create Hamiltonian of one unit cell of the scattering region
        Scatter_UC.CreateHamiltonians( 'toSave', 0);
        Scatter_UC.ShiftCoordinates( obj.shift );
        
        % shifting the coordinates
        coordinates = obj.Scatter_UC.Read('coordinates');
        coordinates = coordinates.Shift( -floor((obj.middle_width - obj.lead_width )/2*[1.5; 0]/3)*3);
        Scatter_UC.Write('coordinates', coordinates );
               
        %applying transverse potential
        obj.ApplyTransversePotential( Scatter_UC )
       
        
       	% apply magnetic field in the unit cell of the scattering region
        % can be applied if the vector potential is identical in each unit cells
        if ~isempty( obj.PeierlsTransform_Scatter ) && obj.Opt.magnetic_field_trans_invariant %for finite q the vector potential must be parallel to q, and perpendicular to the unit cell vector
           	obj.display(['EQuUs:Utils:',class(obj),':CreateScatter: Applying magnetic field in the unit cell of the scattering region']);
           	obj.PeierlsTransform_Scatter.PeierlsTransformLeads( Scatter_UC );
        end

        
        % pin the hole and remove sites beyond the outher boundary
        
        % Create the Hamiltonian of the scattering region
        CreateH = CreateHamiltonians(obj.Opt, obj.param, 'q', obj.q);
        CreateH.CreateScatterH( 'Scatter_UC', Scatter_UC );    
        
        % obtaining coordinates
        coordinates_scatter = CreateH.Read('coordinates');
%        plot( coordinates_scatter.x, coordinates_scatter.y, 'bx')
%        hold on

        % determine sites to be removed
        indexes_hole = (coordinates_scatter.x - obj.cCircle_in.center.x).^2 + (coordinates_scatter.y - obj.cCircle_in.center.y).^2 <= obj.cCircle_in.radius^2;

        % determine sites beyound the circular boundary
        indexes_out =  floor(abs(coordinates_scatter.x - obj.cCircle_out.center.x)/coordinates.b*obj.width*2)*[1.5; 0] >= 1.5*obj.lead_width/2 & ...
            (coordinates_scatter.x - obj.cCircle_out.center.x).^2 + (coordinates_scatter.y - obj.cCircle_out.center.y).^2 > obj.cCircle_out.radius^2;
            
        %removing the sites of the inner hole        
        CreateH.RemoveSites( indexes_hole | indexes_out );
%{
        figure
        plot( coordinates_scatter.x, coordinates_scatter.y, 'bx')
        hold on      
        
        plot( coordinates_scatter.x(~( indexes_out | remove_left | remove_right)), coordinates_scatter.y(~( indexes_out | remove_left | remove_right)), 'kx');
         
        plot( coordinates_scatter.x( keep_left_lead | keep_right_lead), coordinates_scatter.y( keep_left_lead | keep_right_lead), 'rx');      
        
        %removing the sites of the inner hole account for leads        
        %CreateH.RemoveSites( indexes_hole | (indexes_out | remove_left | remove_right ) & not_lead );
%}      
                
        % obtaining the modified coordinates
        coordinates_scatter = CreateH.Read('coordinates');
%        plot( coordinates_scatter.x, coordinates_scatter.y, 'rx')

%         % determine the sheet distance
%         sheet_distance = max(coordinates_scatter.z) - min(coordinates_scatter.z);
%         
%         if obj.num_of_layers > 1
%         % doubling the Hamiltonians and the coordinates
%             Hscatter = [Hscatter, sparse([],[],[],size(Hscatter,1), size(Hscatter,2));
%                     sparse([],[],[],size(Hscatter,1), size(Hscatter,2)), Hscatter];
%                 
%             Hscatter_transverse = [Hscatter_transverse, sparse([],[],[],size(Hscatter_transverse,1), size(Hscatter_transverse,2));
%                     sparse([],[],[],size(Hscatter_transverse,1), size(Hscatter_transverse,2)), Hscatter_transverse];
%                 
%             coordinates_scatter_upper_sheet   = coordinates_scatter;
%             coordinates_scatter_upper_sheet.z = coordinates_scatter_upper_sheet.z + 2*sheet_distance*(obj.num_of_layers-1);
%             coordinates_scatter = coordinates_scatter.Combine( coordinates_scatter_upper_sheet );
%         end
%         
%         CreateH.Write('Hscatter', Hscatter);
%         CreateH.Write('Hscatter_transverse', Hscatter_transverse);
%         CreateH.Write('coordinates', coordinates_scatter);        
        
       	% apply magnetic field in the whole Hamiltonian of the scattering region
        % can be applied for non-translational invariant vector potentials
        if ~isempty( obj.PeierlsTransform_Scatter ) && ~obj.Opt.magnetic_field_trans_invariant 
           	obj.display(['EQuUs:Utils:',class(obj),':CreateScatter: Applying magnetic field in the whole Hamiltonian of the scattering region']);
           	obj.PeierlsTransform_Scatter.PeierlsTransform( CreateH );
        end          
        
        
        obj.CreateH = CreateH;    
        
        % Determine the sites that are coupled to the leads
        y_min = min( coordinates_scatter.y );
        y_max = max( coordinates_scatter.y );
        x_min = min( coordinates_scatter.x );
        x_max = max( coordinates_scatter.x );
        
        non_singular_sites_logical = abs(y_min - coordinates_scatter.y ) < 1e-6 | abs(y_max - coordinates_scatter.y ) < 1e-6 | abs(x_min - coordinates_scatter.x ) < 1e-6 | abs(x_max - coordinates_scatter.x ) < 1e-6;
        
%         figure
%         plot( coordinates_scatter.x, coordinates_scatter.y, 'bx')
%         hold on 
%         plot( coordinates_scatter.x(non_singular_sites_logical), coordinates_scatter.y(non_singular_sites_logical), 'rx');

        non_singular_sites = 1:length(non_singular_sites_logical);
        non_singular_sites = non_singular_sites(non_singular_sites_logical);
        CreateH.Write('kulso_szabfokok', non_singular_sites );

        
    end


%% CreateInterface
%> @brief Creates the Hamiltonians for the interface regions between the leads and scattering center.
%> @param idx Identification number of the interface region. 
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'UseHamiltonian' Logical value. Set true if the interface region should be created to match to the whole Hamiltonian of the scattering center, false (default) if only the surface Green operator of the scattering center is used in the calculations.
	function CreateInterface( obj, idx, varargin )
        
        p = inputParser;
        p.addParameter('UseHamiltonian', false); %true if the interface region is matched to the whole Hamiltonian of the scattering center, false if the surface Green operator of the scattering center is used in the calculations.
        p.parse(varargin{:});
        UseHamiltonian     = p.Results.UseHamiltonian;
        
        %> Hamiltoninans of the interface region		
        Interface_Region = obj.Interface_Regions{idx};
        
        % The regularization of the interface is performed according to the Leads
        Leads = obj.FL_handles.Read( 'Leads' );
        Lead = Leads{idx};
        
        % first retrive the coordinates of the scattering region
        coordinates_scatter = obj.CreateH.Read('coordinates');
        non_singluar_sites_scatter = obj.CreateH.Read('kulso_szabfokok');
        non_singluar_sites_scatter_logical = false(size(non_singluar_sites_scatter));
        non_singluar_sites_scatter_logical(non_singluar_sites_scatter) = true;
        coordinates_scatter_extarnal = coordinates_scatter.KeepSites( non_singluar_sites_scatter_logical );
        
        Interface_Region.Write( 'coordinates2', coordinates_scatter_extarnal );
        Interface_Region.Write( 'K0', Lead.Read('K0'));
        Interface_Region.Write( 'K1', Lead.Read('K1'));
        Interface_Region.Write( 'K1adj', Lead.Read('K1adj'));
        Interface_Region.Write( 'K1_transverse', Lead.Read('K1_transverse'));
        Interface_Region.Write( 'K1_skew_left', Lead.Read('K1_skew_left'));
        Interface_Region.Write( 'K1_skew_right', Lead.Read('K1_skew_right'));
        Interface_Region.Write( 'coordinates', Lead.Read('coordinates'));
        Interface_Region.Write( 'kulso_szabfokok', Lead.Read('kulso_szabfokok'));
        Interface_Region.Write( 'OverlapApplied', true);

        if length(Leads) == 2
            coordinates_shift = [1, -1 ]; %relative to the leads
        else
            coordinates_shift = [0, 0, 1, -1 ]; %relative to the leads
        end
        Interface_Region.ShiftCoordinates( coordinates_shift(idx) );

        % determine the coupling between the interface and the scattering region
        coordinates_interface = Interface_Region.Read('coordinates');

        edge_regions_scatter = sparse(coordinates_scatter.y == min(coordinates_scatter.y) | coordinates_scatter.y == max(coordinates_scatter.y) | coordinates_scatter.x == min(coordinates_scatter.x) | coordinates_scatter.x == max(coordinates_scatter.x) );
        
        edge_regions_interface = sparse(ones(length(coordinates_interface.y),1));% == min(coordinates_interface.y) | coordinates_interface.y == max(coordinates_interface.y) ) | sparse(coordinates_interface.x == min(coordinates_interface.x) | coordinates_interface.x == max(coordinates_interface.x) );
       
        distance_x = ( sparse(coordinates_interface.x)*( edge_regions_scatter' ) - ...
                       edge_regions_interface*( coordinates_scatter.x.*edge_regions_scatter )' ).^2; 
        distance_y = ( sparse(coordinates_interface.y)*( edge_regions_scatter' ) - ...
                       edge_regions_interface*( coordinates_scatter.y.*edge_regions_scatter )' ).^2;
        
        indexes = distance_x + distance_y <= 1.01^2 & distance_x + distance_y > 0;   

%{
         figure1 = figure('rend','painters','pos',[10 10 900 400]);
         plot( coordinates_scatter.x, coordinates_scatter.y, 'bx')
         hold on
         plot( coordinates_interface.x, coordinates_interface.y, 'rx')

             for jdx = 1:size(indexes,1)
                scatter_x = coordinates_scatter.x(indexes(jdx,:));
                scatter_y = coordinates_scatter.y(indexes(jdx,:));
                interface_x = coordinates_interface.x(jdx)*ones(size(scatter_x));
                interface_y = coordinates_interface.y(jdx)*ones(size(scatter_y));
                plot( [scatter_x'; interface_x'], [scatter_y'; interface_y'], 'k' )
             end
        %close(figure1);
%}
        % first determine the coupling constant
        params = Interface_Region.Read('params');
        coupling_constant = params.vargamma;       

        % now construct the coupling matrices
        Hcoupling = sparse(coupling_constant*indexes);
        if ~isempty( coordinates_interface.BdG_u )
            Hcoupling( ~coordinates_interface.BdG_u, ~coordinates_scatter.BdG_u ) = -Hcoupling( ~coordinates_interface.BdG_u, ~coordinates_scatter.BdG_u );
            Hcoupling( ~coordinates_interface.BdG_u, coordinates_scatter.BdG_u ) = 0;
            Hcoupling( coordinates_interface.BdG_u, ~coordinates_scatter.BdG_u ) = 0;
        end
        
        Kcoupling = Hcoupling;
        Kcouplingadj = Kcoupling';
        
        Lead_Orientation = params.Lead_Orientation;
        
        if Lead_Orientation == 1
            
            if UseHamiltonian             
            else
                Kcoupling = Kcoupling(:, non_singluar_sites_scatter_logical);
                Kcouplingadj = Kcouplingadj(non_singluar_sites_scatter_logical, :);                                 
            end
            

            
        elseif Lead_Orientation == -1
            
            if UseHamiltonian  
                non_singular_sites_lead = Lead.Read('kulso_szabfokok');
                Kcoupling = Kcoupling(non_singular_sites_lead, :);
                Kcouplingadj = Kcouplingadj(:, non_singular_sites_lead);
            else
                non_singular_sites_lead = Lead.Read('kulso_szabfokok');
                Kcoupling = Kcoupling(non_singular_sites_lead, non_singluar_sites_scatter);
                Kcouplingadj = Kcouplingadj(non_singluar_sites_scatter, non_singular_sites_lead);
            end
            

        else
            error('EQuUs:Utils:Ribbon:CreateInterface', 'Unknown lead orientation');            
        end
        
               
        Interface_Region.Write('Kcoupling', Kcoupling);
        Interface_Region.Write('Kcouplingadj', Kcouplingadj);  
        
        % method to adjust the interface region and coupling to the scattering region by an external function.
        if ~isempty( obj.interfacemodel )
            obj.interfacemodel( Interface_Region ); 
        end
        
        Interface_Region.Calc_Effective_Hamiltonians( 0, 'Lead', Lead );


    end


end % end methods public



methods (Access=protected)  

    
    
%% createShape
%> @brief Creates the geometry data of the ribbon shaped scattering region.
    function createShape( obj )
        
        if ~isempty( obj.middle_width ) && ~isempty( obj.height )
            obj.param.scatter.shape.width = obj.middle_width;
            obj.param.scatter.shape.height = obj.height;
        end
        
        
        obj.calculate_lead_attach_points();
       
    end
    
    
    
%% calculate_lead_attach_points
%> @brief Determines the site indexes at which the leads are connected to the scattering center.
    function calculate_lead_attach_points( obj )
        for idx = 1:length(obj.param.Leads)
            obj.param.Leads{idx}.M = obj.width;
        end
    end    
    
    %% InputParsing
%> @brief Parses the optional parameters for the class constructor.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'filenameIn' The input filename containing the computational parameters. (Use parameters 'Op' and 'param' instead)
%> @param 'filenameOut' The output filename to export the computational parameters.
%> @param 'WorkingDir' The absolute path to the working directoy.
%> @param 'CustomHamiltoniansHandle' function handle for the custom Hamiltonians. Has the same inputs as #Custom_Hamiltonians.LoadHamiltonians and output values defined by the example #Hamiltonians.
%> @param 'E' The energy value used in the calculations (in the same units as the Hamiltonian).
%> @param 'EF' The Fermi energy in the same units as the Hamiltonian. Attribute #E is measured from this value. (Use for equilibrium calculations in the zero temperature limit. Overrides the one comming from the external source)
%> @param 'silent' Set true to suppress output messages.
%> @param 'leadmodel' A function handle #Lead=f( idx, E, varargin ) of the alternative lead model with equivalent inputs and return values as #Transport_Interface.SurfaceGreenFunctionCalculator and with E standing for the energy.
%> @param 'interfacemodel' A function handle f( #InterfaceRegion ) to manually adjus the interface regions. (Usefull when 'leadmodel' is also given. For example see @InterfaceModel)
%> @param 'Opt' An instance of the structure #Opt.
%> @param 'param' An instance of the structure #param.
%> @param 'q' The transverse momentum quantum number.
%> @param 'mu' An array containing the chemical potentials of the leads.
    function InputParsing( obj, varargin)
    
        p = inputParser;
        p.addParameter('filenameIn', obj.filenameIn, @ischar);
        p.addParameter('filenameOut', obj.filenameOut, @ischar);
        p.addParameter('WorkingDir', obj.WorkingDir, @ischar);
        p.addParameter('E', obj.E, @isscalar);
        p.addParameter('silent', obj.silent);   
        p.addParameter('leadmodel', obj.leadmodel); %individual physical model for the contacts
        p.addParameter('interfacemodel', obj.interfacemodel); %individual physical model for the interface regions
        p.addParameter('Opt', obj.Opt);
        p.addParameter('param', obj.param);
        p.addParameter('q', obj.q);
        p.addParameter('EF', obj.EF);
        p.addParameter('cCircle_in', obj.cCircle_in);
        p.addParameter('cCircle_out', obj.cCircle_out);
        p.addParameter('middle_width', obj.middle_width);
        p.addParameter('lead_width', obj.lead_width);
        
        
        p.addParameter('width', obj.width);
        p.addParameter('height', obj.height); 
        p.addParameter('transversepotential', obj.transversepotential);      
        
        p.parse(varargin{:});
        
        
        InputParsing@NTerminal( obj, 'filenameIn', p.Results.filenameIn, ...
            'filenameOut', p.Results.filenameOut, ...
            'WorkingDir', p.Results.WorkingDir, ...
            'E', p.Results.E, ...
            'silent', p.Results.silent, ...
            'leadmodel', p.Results.leadmodel, ...
            'interfacemodel', p.Results.interfacemodel, ...
            'q', p.Results.q, ...
            'Opt', p.Results.Opt, ...
            'param', p.Results.param, ...
            'EF', p.Results.EF);
        
        InputParsing@Ribbon( obj, 'width', p.Results.width, ...
            'height', p.Results.height, ...
            'transversepotential', p.Results.transversepotential);
        
        obj.middle_width = p.Results.middle_width;
        obj.cCircle_in   = p.Results.cCircle_in;
        obj.cCircle_out  = p.Results.cCircle_out;
        obj.lead_width  = p.Results.lead_width;
        
    end
    

end % methods private

end

