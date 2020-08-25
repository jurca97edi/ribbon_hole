%%    Eotvos Quantum Transport Utilities - Peierls
%    Copyright (C) 2009-2015 Peter Rakyta, Ph.D.
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
%> @addtogroup basic Basic Functionalities
%> @{
%> @file Peierls.m
%> @brief A class responsible for the Peierls and gauge transformations.
%> @} 
%> @brief A class responsible for the Peierls and gauge transformations.
%%
classdef Peierls < handle & Messages

    properties ( Access = protected )
        %> String containing the name of the built-in gauge ('LandauX', 'LandauY').
        gauge
        %> The strength of the magnetic field for the built-in vector potentials.
        eta
        %> Function handle A = f( x,y) of the vector potential to be used in the calculations (A is a N x 2 vector, where N is the number of the points given by the x and y coordinates.)
        Vectorpotential
    end



methods ( Access = public )

%% Contructor of the class
%> @brief Constructor of the class.
%> @param Opt An instance of the structure #Opt.
%> @param varargin Cell array of optional parameters. For details see #InputParsing.
%> @return An instance of the class
    function obj = Peierls(Opt, varargin)
        obj = obj@Messages( Opt );
        
        obj.gauge     = []; % 'LandauX', 'LandauY' gauges or [] for load from input files
        obj.eta       = 0;  % eta_B = 2*pi/phi0*(rCC)^2*B;
        obj.Vectorpotential = [];

        obj.InputParsing( varargin{:} );
        
        if isfield( obj.Opt, 'eta' )
            obj.eta   = obj.Opt.eta;
        else
            obj.eta = 0;
        end    
    
    
        if isempty( obj.Vectorpotential )
            if isfield( obj.Opt, 'gauge' )
                obj.gauge     = obj.Opt.gauge;
            else
                obj.gauge = [];
            end
            obj.Vectorpotential = obj.CreateHandleForVectorpotential();
        end        
        
    end

%% PeierlsTransformLeads
%> @brief Algorithm to perform Peierls transformation in the Hamiltonians of the leads.
%> @param Lead An instance of class #CreateLeadHamiltonians (or a derived class) 
function PeierlsTransformLeads( obj, Lead )
    
    if ~strcmpi( class(Lead), 'CreateLeadHamiltonians' )        
        supClasses = superclasses(Lead);
        if sum( strcmp( supClasses, 'CreateLeadHamiltonians' ) ) == 0
            error(['EQuUs:', class(obj.junction), ':PeierlsTransformLeads'], 'Input Lead is not valid.');
        end
    end

    if Lead.Read( 'MagneticFieldApplied' )
		obj.display(['EQuUs:',class(obj),':PeierlsTransformLeads: Magnetic field already applied in the Hamiltonians']);
		return
    end
    
    if isempty(obj.Vectorpotential)
        obj.display(['EQuUs:',class(obj),':PeierlsTransformLeads: Vectorpotential is not set.']);
        return
    end
		

    if ~Lead.Read( 'HamiltoniansCreated' ) || Lead.Read( 'HamiltoniansDecimated' ) 
       err = MException(['EQuUs:',class(obj),':PeierlsTransformLeads'], 'Hamiltonians in Lead are not created, or they are decimated. Thus, unable to perform Peierls transformation.');
       save('Error_Peierls_PeierlsTransformLeads.mat');
       throw(err); 
    end       
    
    % obtaining Hamiltonians from the leads
    H0 = Lead.Read( 'H0' );
    Lead.Clear( 'H0' );
    H1 = Lead.Read( 'H1' );
    Lead.Clear( 'H1' );
    H1_transverse = Lead.Read( 'H1_transverse' );
    Lead.Clear( 'H1_transverse' );
    H1_skew_right = Lead.Read( 'H1_skew_right' );
    Lead.Clear( 'H1_skew_right' );
    H1_skew_left = Lead.Read( 'H1_skew_left' );
    Lead.Clear( 'H1_skew_left' );
    
    if obj.Opt.BdG
        H0_size = size(H0,1)/2;
        mtx_pair_potential = H0( 1:H0_size, H0_size+1:2*H0_size);
        H0                 = H0( 1:H0_size, 1:H0_size);
        H1                 = H1( 1:H0_size, 1:H0_size);
        if ~isempty(H1_transverse)
            H1_transverse      = H1_transverse( 1:H0_size, 1:H0_size);
        end
        
        if ~isempty(H1_skew_right)
            H1_skew_right      = H1_skew_right( 1:H0_size, 1:H0_size);
        end
        
        if ~isempty(H1_skew_left)
            H1_skew_left      = H1_skew_left( 1:H0_size, 1:H0_size);
        end
        
        if norm(mtx_pair_potential,1) >= 1e-6
            error( ['EQuUs:',class(obj),':PeierlsTransformLeads'], 'The peierls transformation with nonzero B is forbidden in the superconducting phase. Set the pair potential to zero, or use gauge transformation instead.')
            mtx_pair_potential = mtx_pair_potential*0;
        end
    end
    
    
    coordinates = Lead.Read( 'coordinates' );
       
    
    % Peierls in the unit cell
    [sor_idx,oszlop_idx,elements] = find(H0);
        
    startpoint = [coordinates.x( oszlop_idx ), coordinates.y( oszlop_idx )];
    endpoint     = [coordinates.x( sor_idx ), coordinates.y( sor_idx )];
    
    fazis        = obj.Peierls_phase(startpoint,endpoint);        
    fazis_mtx_H0 = sparse(sor_idx,oszlop_idx,fazis,size(H0,1), size(H0,2) );
    elements     = exp(1i*fazis).*elements;
    H0           = sparse(sor_idx,oszlop_idx,elements,size(H0,1), size(H0,2) );
        
    Lead.Write( 'fazis_mtx_H0', fazis_mtx_H0);
    
    
    % Peierls in the coupling of unit cells
    [sor_idx,oszlop_idx,elements] = find(H1);
    orientation = Lead.Read('Lead_Orientation');
    Hanyadik_Lead = Lead.Read('Hanyadik_Lead' );
    if orientation == 1 && ~isempty( Hanyadik_Lead )
        startpoint = [coordinates.x( oszlop_idx ) , coordinates.y( oszlop_idx )];
        endpoint   = [coordinates.x( sor_idx ) - coordinates.a(1), coordinates.y( sor_idx ) - coordinates.a(2)];
    else
        startpoint = [coordinates.x( oszlop_idx ) + coordinates.a(1), coordinates.y( oszlop_idx ) + coordinates.a(2)];
        endpoint   = [coordinates.x( sor_idx ), coordinates.y( sor_idx )];
    end
    
    fazis        = obj.Peierls_phase(startpoint,endpoint);        
    fazis_mtx_H1 = sparse(sor_idx,oszlop_idx,fazis,size(H1,1), size(H1,2) );
    elements     = exp(1i*fazis).*elements;
    H1           = sparse(sor_idx,oszlop_idx,elements,size(H1,1), size(H1,2) );
        
    Lead.Write( 'fazis_mtx_H1', fazis_mtx_H1); 
        
        
    % Peierls in the transverse coupling of unit cells
    if ~isempty(H1_transverse)
        [sor_idx,oszlop_idx,elements] = find(H1_transverse);
        startpoint = [coordinates.x( oszlop_idx ) + coordinates.b(1), coordinates.y( oszlop_idx ) + coordinates.b(2)];
        endpoint     = [coordinates.x( sor_idx ), coordinates.y( sor_idx )];
        
        fazis         = obj.Peierls_phase(startpoint,endpoint);        
        fazis_mtx_H1t = sparse(sor_idx,oszlop_idx,fazis,size(H1_transverse,1), size(H1_transverse,2) );
        elements      = exp(1i*fazis).*elements;
        H1_transverse = sparse(sor_idx,oszlop_idx,elements,size(H1_transverse,1), size(H1_transverse,2) );       
        
        Lead.Write( 'fazis_mtx_H1t', fazis_mtx_H1t); 
    end
    
    % Peierls in the skew couling in the right direction
    if ~isempty(H1_skew_right)
        [sor_idx,oszlop_idx,elements] = find(H1_skew_right);
        orientation = Lead.Read('Lead_Orientation');
        Hanyadik_Lead = Lead.Read('Hanyadik_Lead' );
        if orientation == 1 && ~isempty( Hanyadik_Lead )
            startpoint = [coordinates.x( oszlop_idx ) + coordinates.b(1) , coordinates.y( oszlop_idx ) + coordinates.b(2)];
            endpoint   = [coordinates.x( sor_idx ) + coordinates.a(1), coordinates.y( sor_idx ) + coordinates.a(2)];
        else
            startpoint = [coordinates.x( oszlop_idx ) - coordinates.a(1) + coordinates.b(1), coordinates.y( oszlop_idx ) - coordinates.a(2) + coordinates.b(2)];
            endpoint     = [coordinates.x( sor_idx ), coordinates.y( sor_idx )];
        end        
        
        fazis                   = obj.Peierls_phase(startpoint,endpoint);        
        fazis_mtx_H1_skew_right = sparse(sor_idx,oszlop_idx,fazis,size(H1_transverse,1), size(H1_transverse,2) );
        elements                = exp(1i*fazis).*elements;
        H1_skew_right           = sparse(sor_idx,oszlop_idx,elements,size(H1_transverse,1), size(H1_transverse,2) );       
        
        Lead.Write( 'fazis_mtx_H1_skew_right', fazis_mtx_H1_skew_right); 
    end
    
    % Peierls in the transverse coupling of unit cells
    if ~isempty(H1_skew_left)
        [sor_idx,oszlop_idx,elements] = find(H1_skew_left);
        orientation = Lead.Read('Lead_Orientation');
        Hanyadik_Lead = Lead.Read('Hanyadik_Lead' );
        if orientation == 1 && ~isempty( Hanyadik_Lead )
            startpoint = [coordinates.x( oszlop_idx ) + coordinates.b(1) , coordinates.y( oszlop_idx ) + coordinates.b(2)];
            endpoint     = [coordinates.x( sor_idx ) - coordinates.a(1), coordinates.y( sor_idx ) - coordinates.a(2)];
        else
            startpoint = [coordinates.x( oszlop_idx ) + coordinates.a(1) + coordinates.b(1), coordinates.y( oszlop_idx ) + coordinates.a(2) + coordinates.b(2)];
            endpoint     = [coordinates.x( sor_idx ), coordinates.y( sor_idx )];
        end        
        
        fazis                  = obj.Peierls_phase(startpoint,endpoint);        
        fazis_mtx_H1_skew_left = sparse(sor_idx,oszlop_idx,fazis,size(H1_transverse,1), size(H1_transverse,2) );
        elements               = exp(1i*fazis).*elements;
        H1_skew_left           = sparse(sor_idx,oszlop_idx,elements,size(H1_transverse,1), size(H1_transverse,2) );       
        
        Lead.Write( 'fazis_mtx_H1_skew_left', fazis_mtx_H1_skew_left); 
    end
        
   
    
    if obj.Opt.BdG
        H0 = [H0, mtx_pair_potential; mtx_pair_potential', -conj(H0)];
        H1 = [H1, zeros(size(H1)); zeros(size(H1')), -conj(H1)];
        H1_transverse = [H1_transverse, zeros(size(H1_transverse)); zeros(size(H1_transverse')), -conj(H1_transverse)];
    end   
    
    
     Lead.Write( 'H0', H0 );
     Lead.Write( 'H1', H1 );
     Lead.Write( 'H1adj', H1' );
     Lead.Write( 'H1_transverse', H1_transverse );
     Lead.Write( 'H1_skew_right', H1_skew_right );
     Lead.Write( 'H1_skew_left', H1_skew_left );
     Lead.Write( 'MagneticFieldApplied', 1 );
     %Surface_tmp.Write( 'q', q);
    
    
end

%% PeierlsTransform
%> @brief Algorithm to perform Peierls transformation in Hamiltonian stored in the #CreateHamiltonians object. 
%> @param CreateH An instance of class #CreateHamiltonians (or a derived class).
function PeierlsTransform( obj, CreateH )
    
    if ~strcmpi( class(CreateH), 'CreateHamiltonians' )        
        supClasses = superclasses(CreateH);
        if sum( strcmp( supClasses, 'CreateHamiltonians' ) ) == 0
            error(['EQuUs:', class(obj.junction), ':PeierlsTransform'], 'Input CreateH is not valid.');
        end
    end
    
	if CreateH.Read( 'MagneticFieldApplied' )
		obj.display(['EQuUs:',class(obj),':PeierlsTransform: Magnetic field already applied in the Hamiltonians'], 1);
		return
    end
    
    if isempty(obj.Vectorpotential)
        obj.display(['EQuUs:',class(obj),':PeierlsTransform: Vectorpotential is not set.']);
        return
    end
    
    if ~CreateH.Read( 'HamiltoniansCreated' ) || CreateH.Read( 'HamiltoniansDecimated' )
       err = MException(['EQuUs:',class(obj),':PeierlsTransform'], 'Hamiltonians in are not created, or they are decimated. Thus, unable to perform Peierls transformation.');
       save('Error_Peierls_PeierlsTransform.mat');
       throw(err); 
    end    
    
    NameOfH = 'Hscatter';
    NameOfH_transverse = 'Hscatter_transverse';
    Hamiltoni = CreateH.Read( NameOfH );
    Hamiltoni_transverse = CreateH.Read( NameOfH_transverse );
    CreateH.Clear( NameOfH );
    CreateH.Clear( NameOfH_transverse );
    coordinates =  CreateH.Read( 'coordinates' );
    
    if obj.Opt.BdG
        mtx_pair_potential = Hamiltoni( ~coordinates.BdG_u, coordinates.BdG_u );
        Hamiltoni                 = Hamiltoni( coordinates.BdG_u, coordinates.BdG_u );
        if ~isempty(Hamiltoni_transverse)
            Hamiltoni_transverse      = Hamiltoni_transverse( coordinates.BdG_u, coordinates.BdG_u );
        end
        coordinates.x             = coordinates.x( coordinates.BdG_u );
        coordinates.y             = coordinates.y( coordinates.BdG_u );
        
        if norm(mtx_pair_potential,1) >= 1e-6
            error( 'The peierls transformation with nonzero B is forbidden in the superconducting phase. Set the pair potential to zero, or use gauge transformation instead.')
            mtx_pair_potential = mtx_pair_potential*0;
        end
    end

               
    %Peierls in the scattering region
    [sor_idx,oszlop_idx,elements] = find(Hamiltoni);
    startpoint = [coordinates.x( oszlop_idx ), coordinates.y( oszlop_idx )];
    endpoint     = [coordinates.x( sor_idx ), coordinates.y( sor_idx )];
    fazis          = obj.Peierls_phase(startpoint,endpoint);
        
    fazis_mtx = sparse(sor_idx,oszlop_idx,fazis,size(Hamiltoni,1), size(Hamiltoni,2) );
    elements        = exp(1i*fazis).*elements;
    Hamiltoni = sparse(sor_idx,oszlop_idx,elements,size(Hamiltoni,1), size(Hamiltoni,2) );               
        
    CreateH.Write( 'fazis_mtx_scatter', fazis_mtx);
        
        
    %Peierls in the transverse coupling of the scattering region
    if ~isempty(Hamiltoni_transverse)
        [sor_idx,oszlop_idx,elements] = find(Hamiltoni_transverse);
        startpoint = [coordinates.x( oszlop_idx ) + coordinates.b(1), coordinates.y( oszlop_idx ) + coordinates.b(2)];
        endpoint     = [coordinates.x( sor_idx ), coordinates.y( sor_idx )];
        fazis          = obj.Peierls_phase(startpoint,endpoint);
        
        fazis_mtx_t = sparse(sor_idx,oszlop_idx,fazis,size(Hamiltoni_transverse,1), size(Hamiltoni_transverse,2) );
        elements        = exp(1i*fazis).*elements;
        Hamiltoni_transverse = sparse(sor_idx,oszlop_idx,elements,size(Hamiltoni_transverse,1), size(Hamiltoni_transverse,2) );  
        
        CreateH.Write( 'fazis_mtx_scatter_t', fazis_mtx_t);
    end
            
    
    if obj.Opt.BdG
        Hamiltoni_tmp = sparse([],[],[], size(Hamiltoni,1)*2, size(Hamiltoni,2)*2 );
        Hamiltoni_transverse_tmp = sparse([],[],[], size(Hamiltoni_transverse,1)*2, size(Hamiltoni_transverse,2)*2 );
        Hamiltoni_tmp( coordinates.BdG_u, coordinates.BdG_u ) = Hamiltoni;
        Hamiltoni_tmp( ~coordinates.BdG_u, ~coordinates.BdG_u ) = -conj(Hamiltoni);
        Hamiltoni_tmp( ~coordinates.BdG_u, coordinates.BdG_u ) = mtx_pair_potential;
        Hamiltoni_tmp( coordinates.BdG_u, ~coordinates.BdG_u ) = mtx_pair_potential';
        
        if ~isempty(Hamiltoni_transverse)
            Hamiltoni_transverse_tmp( coordinates.BdG_u, coordinates.BdG_u ) = Hamiltoni_transverse;
            Hamiltoni_transverse_tmp( ~coordinates.BdG_u, ~coordinates.BdG_u ) = -conj(Hamiltoni_transverse);
        end
        
        Hamiltoni            = Hamiltoni_tmp;
        Hamiltoni_transverse = Hamiltoni_transverse_tmp;
    end   
    
    CreateH.Write( NameOfH, Hamiltoni );
    CreateH.Write( NameOfH_transverse, Hamiltoni_transverse );
    CreateH.Write( 'MagneticFieldApplied', 1 );
    %Way2Hamiltonian.Write( 'q', q);
end

%% gaugeTransformation
%> @brief Algorithm to perform gauge transformation on Green's function and/or on Hamiltonian.
%> @param mtx The matrix of the Green's function or the Hamiltonian.
%> @param coordinates An instance of structure #coordinates containing the site coordinates
%> @param gauge_field Function handle S = f( x,y) of the gauge transformation. (S is a N x 1 vector, where N is the number of the points given by the x and y coordinates.)
%> @return Returns with the transformaed matrix.
    function ret = gaugeTransformation( obj, mtx, coordinates, gauge_field ) 
        if isempty( gauge_field )
            ret = mtx;
            return
        end    
        % for Landau_y = Landau_x + grad( gauge_field )
        %gauge_field = @(x,y)( eta*x.*y);
        
        mtx_size = size(mtx,1);
        if length( coordinates.x ) ~= mtx_size
            warning( ['EQuUs:', class(obj.junction), ':gaugeTransformation'], 'The coordinates and matrix elements do not correspond to each other.' )
            return
        end
        
        % "u" and "v" like components must be transformed differently in BdG theory
        if isprop( coordinates, 'BdG_u' ) && ~isempty(coordinates.BdG_u)
            fact = -(-1).^coordinates.BdG_u;
        else
            fact = ones( size(coordinates.x) );
        end
        
                
        if issparse( mtx )
            Umtx = sparse( 1:mtx_size, 1:mtx_size, exp(1i*fact.*gauge_field( coordinates.x, coordinates.y )), mtx_size, mtx_size);
        else
            Umtx = diag( exp(1i*fact.*gauge_field( coordinates.x, coordinates.y )) );
        end
        
        ret = Umtx * mtx * Umtx';        
        
    end

%% gaugeTransformationOnLead
%> @brief Algorithm to perform gauge transformation in the Hamiltonians of a lead. The transformed matrices are stored within the input class.
%> @param Lead An instance of class #CreateLeadHamiltonians (or a derived class) 
%> @param gauge_field Function handle S = f( x,y) of the gauge transformation. (S is a N x 1 vector, where N is the number of the points given by the x and y coordinates.)
    function gaugeTransformationOnLead( obj, Lead, gauge_field ) 
        
        if ~strcmpi( class(Lead), 'CreateLeadHamiltonians' )        
            supClasses = superclasses(Lead);
            if sum( strcmp( supClasses, 'CreateLeadHamiltonians' ) ) == 0
                error(['EQuUs:', class(obj.junction), ':gaugeTransformationOnLead'], 'Input Lead is not valid.');
            end
        end
        
        if isempty( gauge_field )
            return
        end         
        
        if Lead.Read( 'GaugeTransformationApplied' )
            obj.display('Gauge transformation already applied in the Hamiltonians');
            return
        end
        
        H0 = Lead.Read('H0');
        H1 = Lead.Read('H1');
        H1adj = Lead.Read('H1adj');
        orientation = Lead.Read('Lead_Orientation');
        Hanyadik_Lead = Lead.Read('Hanyadik_Lead' );
        coordinates = Lead.Read( 'coordinates' );
        coordinates_tmp = structures( 'coordinates');
        if orientation == 1 && ~isempty( Hanyadik_Lead )
            coordinates_tmp.x = [coordinates.x-coordinates.a(1); coordinates.x];
            coordinates_tmp.y = [coordinates.y-coordinates.a(2); coordinates.y];
        else
            coordinates_tmp.x = [coordinates.x; coordinates.x+coordinates.a(1)];
            coordinates_tmp.y = [coordinates.y; coordinates.y+coordinates.a(2)];
        end
            
            
        
        H0size = size(H0, 1);
        
        if isprop( coordinates, 'BdG_u' )
            coordinates_tmp.BdG_u = [coordinates.BdG_u; coordinates.BdG_u];
        end
        
        H_tmp = [H0,H1;H1adj,H0];%          
        H_tmp = obj.gaugeTransformation( H_tmp, coordinates_tmp, gauge_field );
        
        
        H0 = H_tmp(1:H0size,1:H0size);
        H1 = H_tmp(1:H0size,H0size+1:end);  
        H1adj = H_tmp(H0size+1:end,1:H0size);
        
        Lead.Write( 'H0', H0);
        Lead.Write( 'H1', H1);    
        Lead.Write( 'H1adj', H1adj);
        
        % gauge transformation of the transverse H1
        H1_transverse = Lead.Read('H1_transverse');
        if ~isempty(H1_transverse) && ~Lead.Read('HamiltoniansDecimated')
            coordinates_tmp.x = [coordinates.x; coordinates.x+coordinates.b(1)];
            coordinates_tmp.y = [coordinates.y; coordinates.y+coordinates.b(2)];
            
            H_tmp = [sparse([],[],[],size(H1_transverse,2), size(H1_transverse,1)), H1_transverse; ...
                H1_transverse', sparse([],[],[],size(H1_transverse,1), size(H1_transverse,2))];     
            
            H_tmp = obj.gaugeTransformation( H_tmp, coordinates_tmp, gauge_field );
            H1_transverse = H_tmp( 1:size(H1_transverse,2), size(H1_transverse,1)+1:end );
            Lead.Write( 'H1_transverse', H1_transverse);
        end
        
        
        %gauge transformation of the transverse momentum
%         q = Surface_tmp.Read('q');
%         identity_mtx = sparse(1:size(H1_transverse,1),1:size(H1_transverse,1),1, size(H1_transverse,1), size(H1_transverse,2) );
%         if ~isempty(q) && ~Surface_tmp.Read('HamiltoniansDecimated')
%             coordinates_tmp.x = [coordinates.x+coordinates.b(1); coordinates.x];
%             coordinates_tmp.y = [coordinates.y+coordinates.b(2); coordinates.y];
%             
%             H_tmp = [sparse([],[],[],size(identity_mtx,2), size(identity_mtx,1)), identity_mtx; ...
%                 identity_mtx', sparse([],[],[],size(identity_mtx,1), size(identity_mtx,2))];     
%             
%             H_tmp = gaugeTransformation( H_tmp, coordinates_tmp, gauge_field );
%             identity_mtx = H_tmp( 1:size(H1_transverse,2), size(H1_transverse,1)+1:end );
%             q = q - angle( diag(identity_mtx) );
%             %Surface_tmp.Write( 'q', q);
%         end
        
        
        Lead.Write( 'GaugeTransformationApplied', true);
        
        
    end

%% CreateClone
%> @brief Creates a clone of the present object.
%> @return Returns with the cloned object.
    function ret = CreateClone( obj )        
        ret = Peierls(obj.Opt, ...
            'Vectorpotential', obj.Vectorpotential);             
        
    end   
    
%% setVectorPotential
%> @brief Use to set the handle for the vector potential.
%> @param vectorpot Function handle A = f( x,y) of the vector potential to be used in the calculations (A is a N x 2 vector, where N is the number of the points given by the x and y coordinates.)
function setVectorPotential( obj, vectorpot )
	obj.Vectorpotential = vectorpot;
end  


end % methods public


methods ( Access = protected )

%% Peierls_phase
%> @brief Calculates the Peierls phase along bonds.
%> @param startpoint Coordinates ([x,y]) of the starting points. 
%> @param endpoint Coordinates ([x,y]) of the end points. 
%> @return Returns with the peierls phases
function fazis = Peierls_phase(obj, startpoint, endpoint)
    
    % assume a dirac-delta like magnetic field
    % the vector potential of such a magnetic field can be written as the
    % spatial derivative of a scalar potential A = grad(Phi), 
    % where Phi = C * angle
    
    consts = obj.Vectorpotential(1,1);
    center_x = consts(2);
    center_y = consts(3);

    startpoint_complex  = ( startpoint(:,1) - center_x ) + 1i* ( startpoint(:,2) - center_y );
    endpoint_complex    = ( endpoint(:,1)   - center_x ) + 1i* ( endpoint(:,2) - center_y );
    
    consts = obj.Vectorpotential(1,1);

    fazis = consts(1) * phase(startpoint_complex./ endpoint_complex);

    return

end

%% CreateHandleForVectorpotential
%> @brief Creates function handles for in-built Vector potentials in a gauge given by attribute #Opt.gauge.
%> @return Returns A function handle
function HandleForA = CreateHandleForVectorpotential( obj )    
        
    
    if ~isempty( obj.gauge )  
        if strcmpi( obj.gauge, 'LandauX' )
            HandleForA = @(x,y)([-obj.eta*y;0]);
        elseif strcmpi( obj.gauge, 'LandauY' )
            HandleForA = @(x,y)([0;obj.eta*x]);
        else
            HandleForA = [];
        end
    else
        HandleForA = [];
    end

end


end % methods protected



methods ( Access = protected )


%% InputParsing
%> @brief Parses the optional parameters for the class constructor.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'Vectorpotential' Function handle A = f( x,y) of the vector potential to be used in the calculations (A is a N x 2 vector, where N is the number of the points given by the x and y coordinates.)
function InputParsing( obj, varargin)
    p_validating = inputParser;
    p_validating.addParameter('Vectorpotential', []);
    p_validating.parse(varargin{:});
        
    obj.Vectorpotential = p_validating.Results.Vectorpotential;
    

end

end % methods private

end %Global end


    


