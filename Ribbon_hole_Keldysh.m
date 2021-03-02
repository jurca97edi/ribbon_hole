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
classdef Ribbon_hole_Keldysh < Ribbon_Keldysh & Ribbon_hole

    properties (Access = protected)
    end       
    
    properties (Access = public)    
    end
    


methods ( Access = public )
    
    
%% Contructor of the class
%> @brief Constructor of the class.
%> @param varargin Cell array of optional parameters. For details see #InputParsing.
%> @return An instance of the class
function obj = Ribbon_hole_Keldysh( varargin )
    obj = obj@Ribbon_Keldysh();
    obj = obj@Ribbon_hole();
    
    if strcmpi( class(obj), 'Ribbon_hole_Keldysh') 
        % processing the optional parameters
        obj.InputParsing(varargin{:})                 

        obj.display(['EQuUs:Utils:', class(obj), ':Ribbon_hole_Keldysh: Creating a Ribbon_twoalayer_Keldysh object'])
        
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

%% CustomDysonFunc
%> @brief Custom Dyson function for a two terminal arrangement on a two dimensional lattice.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'gfininv' The inverse of the Greens function of the scattering region. For default the inverse of the attribute #G is used.
%> @param 'constant_channels' Logical value. Set true (default) to keep constant the number of the open channels in the leads for each energy value, or false otherwise.
%> @param 'onlyGinverz' Logical value. Set true to calculate only the inverse of the total Green operator, or false (default) to calculate #G as well.
%> @param 'recalculateSurface' A vector of the identification numbers of the lead surfaces to be recalculated.
%> @param 'decimate' Logical value. Set true (default) to eliminate all inner sites in the Greens function and keep only the selected sites. Set false to omit the decimation procedure.
%> @param 'kulso_szabfokok' Array of sites to be kept after the decimation procedure. (Use parameter 'keep_sites' instead)
%> @param 'selfEnergy' Logical value. Set true to use the self energies of the leads in the Dyson equation, or false (default) to use the surface Green function instead.
%> @param 'keep_sites' Name of sites to be kept in the resulted Green function (Possible values are: 'scatter', 'interface', 'lead').
%> @param 'UseHamiltonian' Set true if the interface region is matched to the whole Hamiltonian of the scattering center, or false (default) if the surface Green operator of the scattering center is used in the calculations.
%> @return [1] The calculated Greens function. 
%> @return [2] The inverse of the Green operator. 
%> @return [3] An instance of structure #junction_sites describing the sites in the calculated Green operator.
    function [Gret, Ginverz, junction_sites] = CustomDysonFunc( obj, varargin ) %NEW output
        
    p = inputParser;
    p.addParameter('gfininv', []);
    p.addParameter('constant_channels', true);
    p.addParameter('onlyGinverz', false );
    p.addParameter('recalculateSurface', [1 2] );
    p.addParameter('decimate', true );
    p.addParameter('kulso_szabfokok', []); %The list of sites to be left after the decimation procedure
    p.addParameter('SelfEnergy', false); %set true to calculate the Dyson equation with the self energy
    p.addParameter('keep_sites', 'lead'); %Name of sites to be kept (scatter, interface, lead)
    p.addParameter('UseHamiltonian', false); %true if the interface region is matched to the whole Hamiltonian of the scattering center, false if the surface Green operator of the scattering center is used in the calculations.
    p.parse(varargin{:});
    gfininv     = p.Results.gfininv;
    constant_channels = p.Results.constant_channels;
    onlyGinverz        = p.Results.onlyGinverz;
    recalculateSurface = p.Results.recalculateSurface;  
    decimate           = p.Results.decimate; 
    kulso_szabfokok = p.Results.kulso_szabfokok;
    useSelfEnergy         = p.Results.SelfEnergy;
    keep_sites        = p.Results.keep_sites;
    UseHamiltonian    = p.Results.UseHamiltonian;
    

    if ~isempty(recalculateSurface)
    
        % creating interfaces for the Leads
        if constant_channels
            shiftLeads = ones(length(obj.param.Leads),1)*obj.E;
        else
            shiftLeads = ones(length(obj.param.Leads),1)*0;
        end
    
        % creating Lead instaces and calculating the retarded surface Green operator/self-energy
        coordinates_shift = [round(obj.height/2)-2, round(obj.height/2)-2, -2, obj.height+1] + obj.shift;             
        obj.FL_handles.LeadCalc('coordinates_shift', coordinates_shift, 'shiftLeads', shiftLeads, 'transversepotential', obj.transversepotential, ...
            'SelfEnergy', useSelfEnergy, 'SurfaceGreensFunction', ~useSelfEnergy, 'gauge_field', obj.gauge_field, 'leads', recalculateSurface, 'q', obj.q, ...
            'leadmodel', obj.leadmodel, 'bias_leads', obj.bias_leads, 'T', obj.T);
        
        for idx = 1:length(recalculateSurface)
            obj.CreateInterface( recalculateSurface(idx), 'UseHamiltonian', UseHamiltonian );
        end
        
    end
    
    [Gret, Ginverz, junction_sites] = CustomDysonFunc@Ribbon_Keldysh(obj, 'gfininv', gfininv, ...
                                    'onlyGinverz', onlyGinverz, ...
                                    'recalculateSurface', [], ...
                                    'decimate', decimate, ...
                                    'kulso_szabfokok', kulso_szabfokok, ...
                                    'SelfEnergy', useSelfEnergy, ...
                                    'keep_sites', keep_sites);                                 

    end
    
%% setEnergy
%> @brief Sets the energy for the calculations
%> @param Energy The value of the energy in the same units as the Hamiltonian.
    function setEnergy( obj, Energy )
        
        setEnergy@Ribbon_Keldysh( obj, Energy );
        
        
    end  
    
    

end % end methods public



methods (Access=protected)    
    
%% CreateHandles
%> @brief Initializes the attributes of the class.
    function CreateHandles( obj )  
        
        CreateHandles@Ribbon_Keldysh( obj )
        
        obj.Scatter_UC = obj.FL_handles.SurfaceGreenFunctionCalculator([], 'createCore', 1, 'q', obj.q, 'T', obj.T, 'mu', obj.mu);
 
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
        p.addParameter('T', obj.T);
        p.addParameter('EF', obj.EF);
        p.addParameter('bias_leads', obj.bias_leads);
        p.addParameter('cCircle_in', obj.cCircle_in);
        p.addParameter('cCircle_out', obj.cCircle_out);
        p.addParameter('middle_width', obj.middle_width);
        p.addParameter('lead_width', obj.lead_width);
        
        p.addParameter('width', obj.width);
        p.addParameter('height', obj.height); 
        p.addParameter('transversepotential', obj.transversepotential);      
        
        p.parse(varargin{:});
        
        
        InputParsing@Ribbon_Keldysh( obj, 'filenameIn', p.Results.filenameIn, ...
            'filenameOut', p.Results.filenameOut, ...
            'WorkingDir', p.Results.WorkingDir, ...
            'E', p.Results.E, ...
            'silent', p.Results.silent, ...
            'leadmodel', p.Results.leadmodel, ...
            'interfacemodel', p.Results.interfacemodel, ...
            'q', p.Results.q, ...
            'Opt', p.Results.Opt, ...
            'param', p.Results.param, ...
            'EF', p.Results.EF, ...
            'bias_leads', p.Results.bias_leads, ...
            'T', p.Results.T,...
            'width', p.Results.width, ...
            'height', p.Results.height, ...
            'transversepotential', p.Results.transversepotential);
        
        InputParsing@Ribbon_hole( obj, 'cCircle_in', p.Results.cCircle_in, 'cCircle_out', p.Results.cCircle_out, 'middle_width', p.Results.middle_width, 'lead_width', p.Results.lead_width );

    end

end % methods private

end

