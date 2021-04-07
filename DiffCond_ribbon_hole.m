%    Three-terminal equilibrium Josephson junction - based on EQuUs v4.9.36
%    Copyright (C) 2015 Peter Rakyta, Ph.D.
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
function diffcond_summed = DiffCond_ribbon_hole( DeltaPhi, height, width , Circ_in , Circ_in2 , EF , resolution, PotentialStrength_lower, PotentialStrength_upper, mu_vec , outputdir)

    if ~exist('filenum', 'var')
        filenum = 1;
    end

    filename = mfilename('fullpath');
    [directory, fncname] = fileparts( filename );

    % The outfilename
    outfilename = [fncname, '_',num2str( filenum )];
    
    % The input and output XML files
    inputXML = 'Graphene_Input_4leads.xml';
    workingdir = [];
    
    % Loading the input parameters
    [Opt, param] = parseInput( inputXML );
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    %dope the leads, usally 1
    dope = 1;%.5;
    mu_N = 0.5;

    % the width in the middle in l.c.
    middle_width = 2*width;
    lead_width = 1.5*width;
    
    % setting the Fermi energy in the leads    
    if(length(param.Leads) == 2)
        param.Leads{1}.epsilon = param.Leads{1}.epsilon - dope;
        param.Leads{2}.epsilon = param.Leads{2}.epsilon - dope;
        param.Leads{1}.M = lead_width;
        param.Leads{2}.M = lead_width;
    else    
        param.Leads{1}.epsilon = param.Leads{1}.epsilon - mu_N; %  dópold
        param.Leads{2}.epsilon = param.Leads{2}.epsilon - mu_N;        
        param.Leads{1}.M = round(40/sqrt(3)); 
        param.Leads{2}.M = round(40/sqrt(3)); 
        param.Leads{3}.epsilon = param.Leads{3}.epsilon - dope;
        param.Leads{4}.epsilon = param.Leads{4}.epsilon - dope;
        param.Leads{3}.M = lead_width;
        param.Leads{4}.M = lead_width;
    end
    
    param.scatter.epsilon = param.scatter.epsilon - mu_N;
    
    %take into account the lattice constant
    R = sqrt((width*1.5)^2+900); 

    %coordinates of inner circle
    cCircle_in = structures('circle');
    cCircle_in.center.x = (lead_width)/2*1.5;
    cCircle_in.center.y = (height)*sqrt(3)/2;
    cCircle_in.radius = Circ_in*R;

    cCircle_in2=Circ_in2*R;
    
    % outher circle
    cCircle_out = structures('circle');
    cCircle_out.center = cCircle_in.center;
    cCircle_out.radius = R;
    
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    % creating the output directory
    setOutputDir()
    
    % the pairing potential
    pair_potential =  max( [param.Leads{3}.pair_potential, param.Leads{4}.pair_potential]);

    mu_leads = [0 0 0 0];
    mu = 0; 
    T = 0;  

    hLeadModel = @LeadModel;

    param.Leads{3}.pair_potential = abs(pair_potential)*exp(1i*DeltaPhi);
    param.Leads{4}.pair_potential = abs(pair_potential);

    % creating the two layer class for Keldysh calculations
    cRibbon_K = Ribbon_hole_Keldysh('width', width, 'height', height, 'Opt', Opt, 'param', param, 'filenameOut', fullfile( outputdir, [outfilename, '.xml']), ...
             'leadmodel', hLeadModel, 'cCircle_in', cCircle_in, 'cCircle_in2', cCircle_in2, 'cCircle_out', cCircle_out, 'middle_width', middle_width, 'lead_width', lead_width, 'WorkingDir',  workingdir, 'bias_leads', mu_leads, 'EF', mu, 'T', T);                 

    ScatterPlot()

    % calculate the Josephson current
    diffcond_summed = CalculateDiffCond(DeltaPhi);

    %save( [outputdir,'/',outfilename, '.mat'], 'diffcond_summed', 'mu_vec', 'Opt', 'param', 'EF', 'Circ_in', 'lead_width', ...
    %                'height', 'PotentialStrength_lower', 'PotentialStrength_upper' );



%% CalculateJosephson_Specq
    function diffcond_summed = CalculateDiffCond(DeltaPhi)
 
        Opt.BdG = 1;

        % chemical potentials in the leads
        mu_leads = [0 0 0 0];

        % chemical potential in the central device
        mu = 0; % CreateHamiltonians is still not adopted to the Keldyhs BdG model

        % set new temperature value
        T = 0;  

        % setting the superconducting pair potentials in the leads
        param.Leads{1}.pair_potential = 0;
        param.Leads{2}.pair_potential = 0;
        param.Leads{3}.pair_potential = abs(pair_potential)*exp(1i*DeltaPhi);
        param.Leads{4}.pair_potential = abs(pair_potential);


        % function handle for the leadmondel
        hLeadModel = @LeadModel;
        
        % create handle for the scatter potential
        hScatterPotential = @ScatterPot;  
        
        % creating the two layer class for Keldysh calculations
        cRibbon_K = Ribbon_hole_Keldysh('width', width, 'height', height, 'Opt', Opt, 'param', param, 'filenameOut', fullfile( outputdir, [outfilename, '.xml']), ...
                 'leadmodel', hLeadModel, 'cCircle_in', cCircle_in, 'cCircle_in2', cCircle_in2, 'cCircle_out', cCircle_out, 'middle_width', middle_width, 'lead_width', lead_width, 'WorkingDir',  workingdir, 'bias_leads', mu_leads, 'EF', mu, 'T', T);                 

        % creating class to calculate the Josephson effect
        %cDiffcond = Diffcond_SSNN( Opt, 'junction', cRibbon_K, 'T', T, 'gfininvfromHamiltonian', true, 'scatterPotential', hScatterPotential );
        cDiffcond = Diffcond_SSNN( Opt, 'junction', cRibbon_K, 'T', T, 'gfininvfromHamiltonian', true );

        diffcond = cDiffcond.diffCond( 'mu_vec', mu_vec );

        diffcond_summed = diffcond;
        
    end

%% LeadModel
%> @brief insert Bz magnetic field into the leads
    function cLead = LeadModel( lead_idx, E, varargin )
        
        p_inner = inputParser;
        p_inner.addParameter('createCore', 0);
        p_inner.addParameter('Just_Create_Hamiltonians', 0);
        p_inner.addParameter('shiftLead', 0);
        p_inner.addParameter('coordinates_shift', 0);
        p_inner.addParameter('transversepotential',[]);
        p_inner.addParameter('gauge_field', [] );% gauge field for performing gauge transformation
        p_inner.addParameter('SelfEnergy',false); % set true to calculate the self energy of the semi-infinite lead
        p_inner.addParameter('SurfaceGreensFunction', true );% set true to calculate the surface Greens function of the semi-infinite lead
        p_inner.addParameter('leadmodel', []); %function handle for an individual physical model for the contacts
        p_inner.addParameter('CustomHamiltonian', []); 
        p_inner.addParameter('q', [] ) %transverse momentum
        p_inner.parse(varargin{:});
        createCore            = p_inner.Results.createCore;
        Just_Create_Hamiltonians = p_inner.Results.Just_Create_Hamiltonians;
        shiftLead             = p_inner.Results.shiftLead;
        coordinates_shift     = p_inner.Results.coordinates_shift;
        transversepotential             = p_inner.Results.transversepotential;
        SelfEnergy            = p_inner.Results.SelfEnergy;
        SurfaceGreensFunction = p_inner.Results.SurfaceGreensFunction;
        q                     = p_inner.Results.q;
        
        Opt2 = Opt;
        Opt2.Silent = true;        
  
        % create class constructing the Hamiltonians of the first layer (composed of two atomic sheets)
        cLead = Lead_Keldysh( Opt2, param, 'hanyadik_lead', lead_idx);
        
        if createCore
            return
        end
         
        % creating Hamiltonians
        cLead.CreateHamiltonians();
        
        % shifting the coordinates according to the position of the lead
        cLead.ShiftCoordinates( coordinates_shift );
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The normal lead 1 - LEFT SIDE
         if( lead_idx == 1 )

            coordinates = cLead.Read('coordinates');

            % Attach lead to the side of the sample 
            % 400=20^2 is the width of the weakly doped region
            tmp = coordinates.y - min(coordinates.y);
            coordinates.y = coordinates.x;
            coordinates.x = tmp-0.5;

            coordinates.x = coordinates.x - floor(( 2*width - lead_width )/2*1.5/3)*3 - 3;
            coordinates.y = coordinates.y + round( ( cCircle_out.center.y - 20)/sqrt(3) )*sqrt(3);

            cLead.Write('coordinates', coordinates);
        
            % "u" and "v" like components must be transformed differently in BdG theory
            if isprop( coordinates, 'BdG_u' ) && ~isempty(coordinates.BdG_u)
                fact = -(-1).^coordinates.BdG_u;
            else
                fact = ones( size(coordinates.x) );
            end
            
            % adding doping potential to the bottom sheet
            H0 = cLead.Read('H0');
            H0 = H0 + sparse( 1:size(H0,1), 1:size(H0,1), fact*PotentialStrength_lower, size(H0,1), size(H0,1));
            cLead.Write('H0', H0);                           

         end
        % RIGHT SIDE
        if( lead_idx == 2 )

            coordinates = cLead.Read('coordinates');

            % Attach lead to the side of the sample 
            tmp = coordinates.y-min( coordinates.y );
            coordinates.y = coordinates.x;
            coordinates.x = tmp-0.5;

            coordinates.x = coordinates.x - floor(( 2*width - lead_width )/2*1.5/3)*3 + middle_width*1.5;% - 3;
            coordinates.y = coordinates.y + round( ( cCircle_out.center.y - 20)/sqrt(3) )*sqrt(3); 

            cLead.Write('coordinates', coordinates);
        
            % "u" and "v" like components must be transformed differently in BdG theory
            if isprop( coordinates, 'BdG_u' ) && ~isempty(coordinates.BdG_u)
                fact = -(-1).^coordinates.BdG_u;
            else
                fact = ones( size(coordinates.x) );
            end
            
            % adding doping potential to the bottom sheet
            H0 = cLead.Read('H0');
            H0 = H0 + sparse( 1:size(H0,1), 1:size(H0,1), fact*PotentialStrength_upper, size(H0,1), size(H0,1));
            cLead.Write('H0', H0);            
        
        end

        
        if Just_Create_Hamiltonians
            return;
        end
        
        
        % Solve the eigen problem
        cLead.TrukkosSajatertekek(E);
        
        % group velocities
        cLead.Group_Velocity();
               
        
        % retarded surface Green operator
        if SurfaceGreensFunction
            cLead.SurfaceGreenFunction();
        end
        
        % retarded SelfEnergy
        if SelfEnergy
            cLead.SelfEnergy();
        end
        
    end 

    function ret = ScatterPot( CreateH , Energy)
    
        coordinates = CreateH.Read('coordinates');
        x = coordinates.x;
        y = coordinates.y;
        
        Hscatter = CreateH.Read('Hscatter');

        % sites on the right side of the scattering r.
        sites2shift_left  = x  < cCircle_in.center.x & abs( y - cCircle_in.center.y ) < 20; 
        sites2shift_right = x >= cCircle_in.center.x & abs( y - cCircle_in.center.y ) < 20; % 1 unit dist. is 0.142 nm
        
        % shift up the on-site energy on the right side and down on the
        % left side with EF
        fact = -(-1).^coordinates.BdG_u;

        Hscatter = Hscatter + sparse(1:length(x),1:length(x), sites2shift_right*( EF + mu_N ).*fact,length(x),length(x)) ...
                            + sparse(1:length(x),1:length(x), sites2shift_left* (-EF + mu_N ).*fact,length(x),length(x));                                

        %no shift
        %Hscatter = Hscatter + sparse(1:length(x),1:length(x), ( sites2shift_right & BdG_u )*( EF + mu_N ) - ( sites2shift_right & ~BdG_u )*( EF + mu_N ),length(x),length(x)) ...
        %                    + sparse(1:length(x),1:length(x), ( sites2shift_left & BdG_u )*(  EF + mu_N ) - ( sites2shift_left & ~BdG_u )*(  EF + mu_N ),length(x),length(x));
                        
        CreateH.Write('Hscatter', Hscatter);
        
        ret = zeros(1,length(coordinates.x));
    end

    function ScatterPlot()

        figure1 = figure( 'Units', 'Pixels', 'Visible', 'on');

        cRibbon_K.getCoordinates();
        CreateH = cRibbon_K.CreateH();
        coord = CreateH.Read('coordinates');
        x=coord.x;
        y=coord.y;

        plot(x,y,'x');
        xlabel('x'); ylabel('y');
        daspect([1 1 1]);

        print('-dpng', [outputdir,'/scatterplot.png'])
        close(figure1);

    end


%% sets the output directory
    function setOutputDir()
        resultsdir = ['Diffcond_H',num2str(height),'_W',num2str(width),'_Cin',num2str(Circ_in),'_Cin2',num2str(Circ_in2),'_EF',num2str(EF),'_res',num2str(resolution)];
        mkdir(resultsdir );
        outputdir = resultsdir;   
        
        workingdir = fullfile(pwd);
        
    end
end
