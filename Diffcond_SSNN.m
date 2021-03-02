%%    Eotvos Quantum Transport Utilities - Josephson_NEGF
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
%> @addtogroup utilities Utilities
%> @{
%> @file Josephson_NEGF.m
%> @brief A class to calculate the DC Josephson current in a non-equilibrium multiterminal setup.
%> @} 
%> @brief A class to calculate the DC Josephson current in a non-equilibrium multiterminal setup.
%%
classdef Diffcond_SSNN < SNSJosephson

    properties    
    end



methods (Access=public)
    
%% Contructor of the class
%> @brief Constructor of the class.
%> @param Opt An instance of the structure #Opt.
%> @param varargin Cell array of optional parameters. For details see #InputParsing.
%> @return An instance of the class
    function obj = Diffcond_SSNN( Opt, varargin )     
        
        obj = obj@SNSJosephson( Opt, varargin{:} ); 
                
        if strcmpi(class(obj), 'Diffcond_SSNN') 
            obj.InputParsing( varargin{:} );
        end
        
        
    end

    
%% Differential conductance
%> @brief Calculates the nonlocal differential conductance
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'mu_db' The number of the energy points in the energy array. (default value is 50)
%> @param 'mu_vec' An array containing the energy points.
%> @param 'mu_min' A minimum of the energy array on the real axis.
%> @param 'mu_max' A maximum of the energy array on the real axis.
%> @return An array of the calculated differential conductance.
    function diffcond_vec = diffCond(obj, varargin ) 
         
        p = inputParser;      
        p.addParameter('mu_db', 50 ); % number of energy points on the contour 
        p.addParameter('mu_vec', [] );
        p.addParameter('mu_min', [] );
        p.addParameter('mu_max', [] );
        
        p.parse(varargin{:});       
        mu_db      = p.Results.mu_db;
        mu_vec     = p.Results.mu_vec;
        mu_min     = p.Results.mu_min;
        mu_max     = p.Results.mu_max;
        
        if ~isempty(mu_vec)
            mu_min = min(real(mu_vec));
            mu_max = max(real(mu_vec)); 
            mu_db = numel(mu_vec);
        elseif (~isempty(mu_max) && ~isempty(mu_min))
        else
            err = MException(['EQuUs:Utils:', class(obj), ':diffCond', 'Unrecognized mu range']);
            save('Error_Josephson_NEGF_diffCond.mat');
            throw(err); 
        end
        
        if mu_min >= mu_max
            diffcond_vec = zeros(3, size(mu_vec));
            return
        end    
               
        % creating the energy array 
        %mu_vec = mu_min:(mu_max-mu_min)/(mu_db-1):mu_max;
        
        
        obj.display(['EQuUs:Utils:', class(obj), ':diffCond: Calculating the differential conductance between mu_min = ', num2str(mu_min),' and mu_max = ',num2str(mu_max),],1);
        diffcond_vec = zeros(3, length(mu_vec));

        % open the parallel pool
        %poolmanager = Parallel( obj.junction.FL_handles.Read('Opt') );
        %poolmanager.openPool(); 

        % setting function handles for the parallel calculations
        hDiffcond_e = @(junction_loc, recalculateSurface)diffCond(junction_loc, recalculateSurface, 'electron');
        %hDiffcond_h = @(junction_loc, recalculateSurface)diffCond(junction_loc, recalculateSurface, 'hole');
        hSetEnergy = @obj.SetEnergy;
        
        junction_loc = obj.junction;
        
        parfor iidx = 1:length(mu_vec)
        %for iidx = 1:length(mu_vec)

            % the self energies of all terminals should be recalculated
            recalculateSurface = 1:length(junction_loc.Interface_Regions);
        
            %%%%%%%%%%% Calculate the electron part %%%%%%%%%%%%%
            % set the energy value in the current iteration
            feval(hSetEnergy, mu_vec(iidx), 'junction', junction_loc );
            
            % calculating the Green operator of the scattering center
            obj.create_scatter_GreensFunction('junction', junction_loc);
                               
            % evaluate the current at a complex energy
            diffcond_vec(:, iidx) = 2*feval(hDiffcond_e, junction_loc, recalculateSurface);
            
%             %%%%%%%%%%% Calculate the hole part %%%%%%%%%%%%%
%             % set the energy value in the current iteration
%             feval(hSetEnergy, -mu_vec(iidx), 'junction', junction_loc );
%             
%             % calculating the Green operator of the scattering center
%             obj.create_scatter_GreensFunction('junction', junction_loc);
%                                
%             % evaluate the current at a complex energy
%             diffcond_vec(iidx) = diffcond_vec(iidx) + feval(hDiffcond_h, junction_loc, recalculateSurface);

        end
        
        % closing the parallel pool
        %poolmanager.closePool();

        %--------------------------------------------------
        % nested functions
		% calculates the differential conductance
		function diffcond = diffCond( junction_loc, recalculateSurface, chargeType )  
            
            Dysonfunc = @()(junction_loc.CustomDysonFunc( ...
                'decimate', false, 'constant_channels', false, ...
                'recalculateSurface', recalculateSurface, ...
                'SelfEnergy', true));
            
            %Calculate retarded and advanced Green-functions of the whole junction at the surface points
            % and get the junction surface sites
            [Gret, ~, junction_sites] = junction_loc.FL_handles.DysonEq( 'CustomDyson', Dysonfunc );
            
            % The advanced Greens function of the whole junction at the surface points
            Gav=Gret';
            
            % The retarded Greens function of "normal lead 2 ---> scattering region"
            Gret_sc_lead2 = Gret(junction_sites.Scatter.site_indexes, junction_sites.Interface{2}.site_indexes);
            
            % The advanced Greens function of "normal lead 1 ---> normal lead 2"
            Gav_lead_2_lead_1 = Gav(junction_sites.Interface{2}.site_indexes,junction_sites.Interface{1}.site_indexes);
            
            % Getting retarded self energy of normal lead 2
            Leads=junction_loc.FL_handles.Read('Leads');
            sigma_NormalLead2=Leads{2}.Read('Sigma');
            
            % Get the number of nonsingular sites in the lead2
            Neffs = junction_loc.FL_handles.Get_Neff(); 
            Neff = Neffs(2);
                        
            %Create retarded and advanced Green-functions of the electron-like particles in the normal lead 2
            coordinates_NormalLead_2 = Leads{2}.Read('coordinates');
            
            % getting the nonsingular sites of the second lead
            non_singular_sites_lead2 = Leads{2}.Read('kulso_szabfokok');                        
            
            % Turn the hole/eelctron part of the self energy to zero
            if strcmp( chargeType, 'electron' )
                sigma_NormalLead2(~coordinates_NormalLead_2.BdG_u(non_singular_sites_lead2), ~coordinates_NormalLead_2.BdG_u(non_singular_sites_lead2)) = 0; 
                fact = -1;
            elseif strcmp( chargeType, 'hole' ) 
                sigma_NormalLead2(coordinates_NormalLead_2.BdG_u(non_singular_sites_lead2), coordinates_NormalLead_2.BdG_u(non_singular_sites_lead2)) = 0;     
                fact = 1;
            else
                error('Unknown charge carrier type')
            end
            
            % The electron part of the lesser self energy of the normal lead 2
            %Sigma_lesser = zeros( length(coordinates_NormalLead_2.BdG_u) );
            %Sigma_lesser(1:Neff, 1:Neff) = sigma_NormalLead2' - sigma_NormalLead2;
            Sigma_lesser = sigma_NormalLead2' - sigma_NormalLead2;

         
            %Create current operator for normal lead 1
            %[~, ~, ~, K1_sc, ~] = junction_loc.Interface_Regions{1}.Get_Effective_Hamiltonians(); 
            [~, ~, ~, K1_sc, ~, coordinates_interface_1] = junction_loc.Interface_Regions{1}.Get_Effective_Hamiltonians(); 
            
            % apply tau_3 to the current operator
            K1_sc(~coordinates_interface_1.BdG_u, :) = -K1_sc(~coordinates_interface_1.BdG_u, :);
            
            %Calc diff.cond.
            %diffcond_total = fact*(2/pi)*real(trace(K1_sc*Gret_sc_lead2*Sigma_lesser*Gav_lead_2_lead_1));
            
            %Calc electron part and hole part of the diff. cond.
            K1_sc_e = K1_sc;
            K1_sc_e(~coordinates_interface_1.BdG_u, :) = 0;
            
            K1_sc_h = K1_sc;
            K1_sc_h(coordinates_interface_1.BdG_u, :) = 0;
            
            diffcond_e = fact*(2/pi)*real(trace(K1_sc_e*Gret_sc_lead2*Sigma_lesser*Gav_lead_2_lead_1));
            diffcond_h = fact*(2/pi)*real(trace(K1_sc_h*Gret_sc_lead2*Sigma_lesser*Gav_lead_2_lead_1));
            
            diffcond = [diffcond_e+diffcond_h; diffcond_e; diffcond_h];

        end

    end
 

    
    end % public methods
    
    
    methods (Access=protected)    


    

%% InputParsing
%> @brief Parses the optional parameters for the class constructor.
%> @param varargin Cell array of optional parameters (https://www.mathworks.com/help/matlab/ref/varargin.html):
%> @param 'T' The temperature in Kelvin (scalar or an array) 
%> @param 'mu' The Chemical potential in the same unit as other energy scales in the Hamiltonians.
%> @param 'scatterPotential' A function handle pot=f( #coordinates ) or pot=f( #CreateHamiltonians, Energy) for the potential to be applied in the Hamiltonian (used when FiniteGreensFunctionFromHamiltonian=true).
%> @param 'gfininvfromHamiltonian' Logical value. Set true calculate the surface Greens function of the scattering region from the Hamiltonaian of the scattering region, or false (default) to calculate it by the fast way (see Phys. Rev. B 90, 125428 (2014) for details).
%> @param 'junction' An instance of class #NTerminal (or its derived class) representing the junction.
    function InputParsing( obj, varargin )
                
        p = inputParser;
        p.addParameter('T', obj.T);
        p.addParameter('mu', obj.mu, @isscalar);
        p.addParameter('scatterPotential', obj.scatterPotential);
        p.addParameter('gfininvfromHamiltonian', obj.gfininvfromHamiltonian);
        p.addParameter('junction', obj.junction);     
        
        p.parse(obj.varargin{:});
            
        
        InputParsing@SNSJosephson( obj, 'T', p.Results.T, ...
                                'mu', p.Results.mu, ...
                                'scatterPotential', p.Results.scatterPotential, ...
                                'gfininvfromHamiltonian', p.Results.gfininvfromHamiltonian, ...
                                'useSelfEnergy', true, ...
                                'junction', p.Results.junction);
               

    end 
    
    
end % protected methods
end
