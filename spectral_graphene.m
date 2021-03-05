function spectral_graphene( height, width, flux , Circ_in, Circ_out , EF , location, resolution)

if ~exist('filenum', 'var')
    filenum = 1;
end
	disp(['calculating ABS spectra as a function of the transverse momentum with file extension ', num2str(filenum)]);

    filename = mfilename('fullpath');
    [directory, fncname] = fileparts( filename );

    % geometry
    if ~exist('height', 'var')
        height = 40;% 0.25 micron
    end
    
    if ~exist('width', 'var')
        width = 10;% 0.25 micron
    end
    
    if ~exist('flux', 'var')
        flux = 0;
    end
    
    if ~exist('Circ_in', 'var')
        Circ_in = 0.7;
    end
    
    if ~exist('Circ_out', 'var')
        Circ_out = 1.5;
    end
    
    if ~exist('EF', 'var')
        EF = 0.2;
    end
    
    if ~exist('location', 'var')
        location = 'outer';
    end
    
    if ~exist('resolution', 'var')
        resolution = '80';
    end
    
    height=str2num(height);
    width=str2num(width);
    flux=str2num(flux);
    Circ_in=str2num(Circ_in);
    Circ_out=str2num(Circ_out);
    EF=str2num(EF);
    resolution=str2num(resolution);
    
    % imaginary part of the Fermi energy
    eta = 1e-7;
    
    % The outfilename
    outfilename = [fncname, '_',num2str( filenum )];
    outputdir   = [];
    
    % The input and output XML files
    inputXML = 'Graphene_Input_4leads.xml';
    outputXML = [outfilename,'.xml'];
    
    % Loading the input parameters
    [Opt, param] = parseInput( inputXML );
    
    
    %dope the leads, usally 1
    dope = 1;%.5;
    mu_N = 0.5;

    % class representing the two terminal ribbon structure

    % the width in the middle in l.c.
    middle_width = 2*width;
    lead_width = 1.5*width;

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

    %take into account the lattice constant
    R = sqrt((width*1.5)^2+900); 

    %coordinates of inner circle
    cCircle_in = structures('circle');
    cCircle_in.center.x = lead_width/2*1.5;
    cCircle_in.center.y = (height)*sqrt(3)/2;
    cCircle_in.radius = Circ_in*R;

    % outher circle
    cCircle_out = structures('circle');
    cCircle_out.center = cCircle_in.center;
    cCircle_out.radius = R;

    param.scatter.epsilon = param.scatter.epsilon - mu_N;

    Delta = max( [param.Leads{3}.pair_potential, param.Leads{4}.pair_potential]);
    
    setOutputDir();       
       
    % class representing the two terminal ribbon structure
    cRibbon = [];
    
    % setting the energy and tranverse momentum grids for the calculations
    [Evec, phivec] = setVectors();

    
    % The calculated density of states
    density_of_states_right_electron = zeros( length(Evec), length(phivec));
    density_of_states_right_hole = zeros( length(Evec), length(phivec));
    density_of_states_left_electron = zeros( length(Evec), length(phivec));
    density_of_states_left_hole = zeros( length(Evec), length(phivec));
    
    polarization_right = zeros( length(Evec), length(phivec));
    polarization_left = zeros( length(Evec), length(phivec));

    % creating function handle for the Hamiltonians
    Opt.BdG = false;

    cRibbon=[];
    
    hLeadModel = @LeadModel;

    % creating the NTerminal class
    cRibbon = Ribbon_hole('width', width, 'height', height, 'Opt', Opt, 'param', param, 'filenameOut', fullfile( outputdir, [outfilename, '.xml']), ...
             'leadmodel', hLeadModel, 'cCircle_in', cCircle_in, 'cCircle_out', cCircle_out, 'middle_width', middle_width, 'lead_width', lead_width);     %{
                       
    %{
    cRibbon.Transport(1e-2, 'gfininvfromHamiltonian', true);                       
    
    CreateH = cRibbon.CreateH();
    coord = CreateH.Read('coordinates');
    x=coord.x;
    y=coord.y;
    
    plot(x,y,'x');
    daspect([1 1 1]);
    
    [hLead, hScatter, hgauge_field ] = createVectorPotential( flux );
    test_continuity( hLead, hScatter, hgauge_field, x , y ,cCircle_in, flux);
    %}
                                 
    % plot the coordinates of the scattering region
%    ScatterPlot();
    %return            
    
    % calculate the spectra of a slice with the give paramteres
    %Spectra(); return
    
    Opt.BdG = true;
    
    % open parallel pool (if not opened already)
    parallelmanager = Parallel( Opt );
    %parallelmanager.openPool();
    
    tic
    
    % calculate the density of states
    CalculateSpectralDensity();
    
    toc
    
    % close parallel pool
    %parallelmanager.closePool();  
    
    save( [outputdir,'/',outfilename, '.mat'], 'Evec', 'height', 'EF', 'density_of_states_right_electron', 'phivec', ...
                'density_of_states_right_hole','density_of_states_left_electron','density_of_states_left_hole','flux' );

    %EgeszAbraP(); 

%% SetVectors
%> @brief Sets the energy and transverse quantum number grid points
%> @return [1] 1D array of the energy points.
%> @return [2] 1D array of the transverse momentum points.
    function [Evec, phivec] = setVectors()
        
        %Delta_E = 0.001;
        %E_window = max(abs(Delta))/5;
        %Evec = EF - E_window:2*E_window/resolution:EF + E_window;
        %if max(abs(Delta)) < 0.1
        Evec = 0:max(abs(Delta))*1.05/resolution:1.05*max(abs(Delta));
        %end
        
        phivec_length = pi;
        phivec = 0:phivec_length/2:phivec_length;
        %phivec = 0:phivec_length/(resolution-1):phivec_length;
        %phivec = [pi/3,2*pi/3];
    end

%% CalculateSpectralDensity
%> @brief Calculate the density of states
    function CalculateSpectralDensity()

        % Do the calculations
        for idx = 1:length(phivec)
            phi = phivec(idx);
            
            % setting the phase difference
            param.Leads{2}.pair_potential = abs(param.Leads{2}.pair_potential)*exp(1i*phi);
            param.Leads{1}.pair_potential = abs(param.Leads{1}.pair_potential);
			
            disp(' ')
            disp(' ')
            disp(' ----------------------------------')
            disp(['Calculating Greens function for the phase difference phi=', num2str(phi)]) 
                
            % temporary array to store the calculated DOS
            Rho_right_electron = zeros( length(Evec), 1);
            Rho_right_hole     = zeros( length(Evec), 1);
            Rho_left_electron = zeros( length(Evec), 1);
            Rho_left_hole     = zeros( length(Evec), 1);
            
            % creating the NTerminal class
            cRibbon = Ribbon_hole('width', width, 'height', height, 'Opt', Opt, 'param', param, 'filenameOut', fullfile( outputdir, [outfilename, '.xml']), ...
                     'leadmodel', hLeadModel, 'cCircle_in', cCircle_in, 'cCircle_out', cCircle_out, 'middle_width', middle_width, 'lead_width', lead_width);                  
            % functio handle to pick the central sites in the scattering region
            switch(location)
                case 'outer'
                    hChoseSites = @ChoseSitesOuter;
                case 'center'
                    hChoseSites = @ChoseSitesCenter;
                case 'inner'
                    hChoseSites = @ChoseSitesInner;
            end

            hScatterPot = @ScatterPot;

            for jdx=1:length(Evec)
            %parfor jdx=1:length(Evec)
                 Energy = Evec(jdx);                 

                 % create an instance of class DOS to calculate the density of states along the whole scattering region
                 DOS_handles = DOS( Opt, 'junction', cRibbon, 'useSelfEnergy', false, 'scatterPotential', hScatterPot);

                 % calculate the local DOS along the whole scattering region
                 cLocalDOS = DOS_handles.LocalDOSCalc( Energy+1i*eta, 'ChoseSites', hChoseSites  );

                 % determine the sites of the left side
                 left_sites = cLocalDOS.coordinates.x < cCircle_in.center.x;              

                 % Right side
                 Rho_right_electron(jdx) = sum(cLocalDOS.DOSvec( (~left_sites) & cLocalDOS.coordinates.BdG_u ));                   
                 Rho_right_hole(jdx)     = sum(cLocalDOS.DOSvec( (~left_sites) & (~cLocalDOS.coordinates.BdG_u) ));                   

                 % Left side
                 Rho_left_electron(jdx) = sum(cLocalDOS.DOSvec( (left_sites) & cLocalDOS.coordinates.BdG_u ));                   
                 Rho_left_hole(jdx)     = sum(cLocalDOS.DOSvec( (left_sites) & (~cLocalDOS.coordinates.BdG_u) ));  

%                 CreateH = cRibbon.CreateH();
%                 H = CreateH.Read('Hscatter');

%                % setting the current energy
%                cRibbon.setEnergy( Energy+1i*eta );
%                
%                % Calculates the DOS for the given energy and transverse momentum quantum number
%                [Atot, G] = cRibbon.CalcSpectralFunction( Energy+1i*eta, ...
%                'gfininvfromHamiltonian', true, 'PotInScatter', hScatterPotential, ...
%                'decimateDyson', false);
%                
%                % determine the coordinates of the surface points of the scattering region
%                coordinates = cRibbon.getCoordinates();
%                
%                % determine the DOS on the points of the scattering reigon
%                DOSvec = -imag(diag(G))/pi;
%                DOSvec = DOSvec( 1:length(coordinates.z) );
%                
%                % determine the sites of the bottom layer
%                left_sites = abs(coordinates.z) < 0.5*max(coordinates.z); 
%                
%                % determine the Density of states for the right/left electrons/holes
%                Rho_right_electron(jdx) = sum(DOSvec( (~left_sites) & coordinates.BdG_u ));                   
%                Rho_right_hole(jdx)     = sum(DOSvec( (~left_sites) & (~coordinates.BdG_u) ));                   
%                Rho_left_electron(jdx) = sum(DOSvec( (left_sites) & coordinates.BdG_u ));                   
%                Rho_left_hole(jdx)     = sum(DOSvec( (left_sites) & (~coordinates.BdG_u) ));   

            end

            density_of_states_right_electron(:,idx) = Rho_right_electron;
            density_of_states_right_hole(:,idx)     = Rho_right_hole;
%            denominator = Rho_right_electron + Rho_right_hole;
%            denominator ( denominator == 0 ) = 1;
%            polarization_right(:,idx) = ( Rho_right_electron - Rho_right_hole )./ denominator ;

            density_of_states_left_electron(:,idx) = Rho_left_electron;
            density_of_states_left_hole(:,idx)     = Rho_left_hole;
%            denominator = Rho_left_electron + Rho_left_hole;
%            denominator ( denominator == 0 ) = 1;
%            polarization_left(:,idx) = ( Rho_left_electron - Rho_left_hole )./ denominator ;

            % calculated data are plotted at each step
            EgeszAbra()

            %EgeszAbra2(idx);

            % exporting the calculated data
            save( [outputdir,'/',outfilename, '.mat'], 'Evec', 'height', 'EF', 'density_of_states_right_electron', 'phivec', ...
                'density_of_states_right_hole','density_of_states_left_electron','density_of_states_left_hole','flux' );

            % displaying the status of the calculations
            disp([ num2str(idx/length(phivec)*100),' % calculated of the density of states.'])
        end
     
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
        transversepotential   = p_inner.Results.transversepotential;
        SelfEnergy            = p_inner.Results.SelfEnergy;
        SurfaceGreensFunction = p_inner.Results.SurfaceGreensFunction;
        q                     = p_inner.Results.q;
        CustomHamiltonian     = p_inner.Results.CustomHamiltonian;
        leadmodel             = p_inner.Results.leadmodel;
    
        Opt2 = Opt;
        Opt2.magnetic_field = false;
        Opt2.Silent = true;
        
        if (lead_idx > 2)
            
            cTransport_Interface = Transport_Interface(E, Opt2, param );
                
            cLead = cTransport_Interface.SurfaceGreenFunctionCalculator( lead_idx, 'createCore', createCore, ...
                        'Just_Create_Hamiltonians', Just_Create_Hamiltonians, ...
                        'shiftLead', shiftLead, ...
                        'coordinates_shift', coordinates_shift, ...
                        'transversepotential', transversepotential, ...
                        'SelfEnergy', SelfEnergy, ...
                        'SurfaceGreensFunction', SurfaceGreensFunction);
                    
        else
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
                H0 = H0 + sparse( 1:size(H0,1), 1:size(H0,1), fact*0, size(H0,1), size(H0,1));
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
                H0 = H0 + sparse( 1:size(H0,1), 1:size(H0,1), fact*0, size(H0,1), size(H0,1));
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
        
    end

    function ret = ScatterPot( CreateH, Energy )
    
        coordinates = CreateH.Read('coordinates');
        x = coordinates.x;
        y = coordinates.y;
        
        Hscatter = CreateH.Read('Hscatter');

        % sites on the right side of the scattering r.
        sites2shift_left  = x  < cCircle_in.center.x & abs( y - cCircle_in.center.y ) < 20;
        sites2shift_right = x >= cCircle_in.center.x & abs( y - cCircle_in.center.y ) < 20;

        % shift up the on-site energy on the right side and down on the
        % left side with EF
        fact = -(-1).^coordinates.BdG_u;
        
        Hscatter = Hscatter + sparse(1:length(x),1:length(x), sites2shift_right*( EF + mu_N ).*fact,length(x),length(x)) ...
                            + sparse(1:length(x),1:length(x), sites2shift_left* (-EF + mu_N ).*fact,length(x),length(x));                                

        %no shift
        %Hscatter = Hscatter + sparse(1:length(x),1:length(x), ( sites2shift_right & BdG_u )*( EF + mu_N ) - ( sites2shift_right & ~BdG_u )*( EF + mu_N ),length(x),length(x)) ...
        %                    + sparse(1:length(x),1:length(x), ( sites2shift_left & BdG_u )*(  EF + mu_N ) - ( sites2shift_left & ~BdG_u )*(  EF + mu_N ),length(x),length(x));
                        
        CreateH.Write('Hscatter', Hscatter);
%{
        figure1 = figure('rend','painters','pos',[10 10 900 400]);
  
        plot(coordinates.x,coordinates.y,'x','MarkerSize',5,'color','k');
        hold on
        plot(coordinates.x( sites2shift_right ),coordinates.y(sites2shift_right),'x','MarkerSize',5,'color','r');
        plot(coordinates.x( sites2shift_left ),coordinates.y(sites2shift_left),'x','MarkerSize',5,'color','b');
        daspect([1 1 1])

        fontsize=12;    
        xlabel('x [a]','FontSize', fontsize,'FontName','Times New Roman');
        ylabel('y [a]','FontSize', fontsize,'FontName','Times New Roman'); 
        
        legend('Scattering region','Shift by +EF', 'Shift by -EF');
        
        print('-dpng', fullfile(outputdir,['scatterplot']))
        %close(figure1);
%}        
        ret = zeros(1,length(coordinates.x));
    end
 

%% CreateHandlesForMagneticField
%> @brief Creates and set function handles of the magnetic vector potentials in the Ribbon class
%> @param B The magnetic field
    function CreateHandlesForMagneticField( flux )        
        [hLead, hScatter, gauge_field ] = createVectorPotential( flux );
        cRibbon.setHandlesForMagneticField( 'scatter', hScatter );
    end

%% createVectorPotential
%> @brief Creates the function handles of the magnetic vector potentials and the gauge field
%> @param B The magnetic field
%> @return [1] Function handle of the vector potential in the leads .
%> @return [2] Function handle of the vector potential in the scattering region.
%> @return [3] Function handle of the scalar potential that transforms vector potential in the scattering center into a vector potential in the leads.
    function [hLead, hScatter, hgauge_field ] = createVectorPotential( flux )
        
        % creting the funciton handles of the vector potentials
        hLead = @(x,y)(CircularVectorPotential(x,y, 0, cCircle_in));
        hScatter = @(x,y)(ConstantVectorPotential(x,y, flux , cCircle_in ));
        hgauge_field = @(x,y)(CircularGaugeField(x,y, flux, cCircle_in));

    end

%% plotfunction
    function EgeszAbra()
        
% ********************** plot the DOS ***********************

        % determine points to be plotted
        indexes = logical( density_of_states_right_electron(1,:));
        
        if length(phivec(indexes)) < 2
            return
        end
        
        % creating figure in units of pixels
        figure1 = figure( 'Units', 'Pixels', 'Visible', 'off', 'Colormap', ... % summer, 'pos',[10 10 900 700] );%, ...
    [0.0416666679084301 0 0;0.0833333358168602 0 0;0.125 0 0;0.16666667163372 0 0;0.20833332836628 0 0;0.25 0 0;0.291666656732559 0 0;0.333333343267441 0 0;0.375 0 0;0.416666656732559 0 0;0.458333343267441 0 0;0.5 0 0;0.541666686534882 0 0;0.583333313465118 0 0;0.625 0 0;0.666666686534882 0 0;0.708333313465118 0 0;0.75 0 0;0.791666686534882 0 0;0.833333313465118 0 0;0.875 0 0;0.916666686534882 0 0;0.958333313465118 0 0;1 0 0;1 0.0416666679084301 0;1 0.0833333358168602 0;1 0.125 0;1 0.16666667163372 0;1 0.20833332836628 0;1 0.25 0;1 0.291666656732559 0;1 0.333333343267441 0;1 0.375 0;1 0.416666656732559 0;1 0.458333343267441 0;1 0.5 0;1 0.541666686534882 0;1 0.583333313465118 0;1 0.625 0;1 0.666666686534882 0;1 0.708333313465118 0;1 0.75 0;1 0.791666686534882 0;1 0.833333313465118 0;1 0.875 0;1 0.916666686534882 0;1 0.958333313465118 0;1 1 0;1 1 0.0625;1 1 0.125;1 1 0.1875;1 1 0.25;1 1 0.3125;1 1 0.375;1 1 0.4375;1 1 0.5;1 1 0.5625;1 1 0.625;1 1 0.6875;1 1 0.75;1 1 0.8125;1 1 0.875;1 1 0.9375;1 1 1]);
        
        % font size on the figure will be 16 points
        fontsize = 12;        

        % define the colorbar limits
        colbarlimits = [min(min(real(log(density_of_states_right_electron)))) max(max(real(log(density_of_states_right_electron))))];
        
        % define the axis limits
        x_lim = [min(phivec) max(phivec)];
        y_lim = [min(Evec) max(Evec)]/min(abs(Delta));
        %y_lim = [min(Evec) max(Evec)]/0.1;%the usual value of Delta
        
        
        axes_DOS_right_electron = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...      
                'CLim', colbarlimits, ...
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        [X,Y] = meshgrid( phivec(indexes), Evec );
        levelnum = 50;
        density_of_states2plot = real(log(density_of_states_right_electron(:,indexes)));
        db=1;
        contour(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_right_electron);
       
        
        % Create xlabel
        xlabel('\Delta\theta RE','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_right_electron);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_right_electron);
        
        
        
        %---------------------------------------------------------------

        axes_DOS_right_hole = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...      
                'CLim', colbarlimits, ...
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        [X,Y] = meshgrid( phivec(indexes), Evec );
        levelnum = 50;
        density_of_states2plot = real(log(density_of_states_right_hole(:,indexes)));
        db=1;
        contour(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_right_hole);
        
        % Create xlabel
        xlabel('\Delta\theta RH','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_right_hole);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_right_hole);
        
        
        
        %---------------------------------------------------------------

        axes_DOS_left_electron = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...      
                'CLim', colbarlimits, ...
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        [X,Y] = meshgrid( phivec(indexes), Evec );
        levelnum = 50;
        density_of_states2plot = real(log(density_of_states_left_electron(:,indexes)));
        db=1;
        contour(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_left_electron);
        
        % Create xlabel
        xlabel('\Delta\theta LE','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left_electron);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left_electron);
        
        
        
        %---------------------------------------------------------------

        axes_DOS_left_hole = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...      
                'CLim', colbarlimits, ...
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        [X,Y] = meshgrid( phivec(indexes), Evec );
        levelnum = 50;
        density_of_states2plot = real(log(density_of_states_left_hole(:,indexes)));
        db=1;
        contour(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_left_hole);
        
        % Create xlabel
        xlabel('\Delta\theta LH','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left_hole);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left_hole);
       
        
        
        % setting figure position
        figure_pos = get( figure1, 'Position' );
        
        %set the position of the axes_DOS
        OuterPosition = get(axes_DOS_right_electron, 'OuterPosition');
        OuterPosition(1) = 0;
        OuterPosition(2) = figure_pos(4)/2;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_right_electron, 'OuterPosition', OuterPosition);  
        
        %set the position of the axes_pot
        OuterPosition = get(axes_DOS_right_hole, 'OuterPosition');
        OuterPosition(1) = figure_pos(3)/2;
        OuterPosition(2) = figure_pos(4)/2;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_right_hole, 'OuterPosition', OuterPosition);  
        
        
        
        %set the position of the axes_DOS
        OuterPosition = get(axes_DOS_left_electron, 'OuterPosition');
        OuterPosition(1) = 0;
        OuterPosition(2) = 0;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_left_electron, 'OuterPosition', OuterPosition); 
        
        %set the position of the axes_pot
        OuterPosition = get(axes_DOS_left_hole, 'OuterPosition');
        OuterPosition(1) = figure_pos(3)/2;
        OuterPosition(2) = 0;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_left_hole, 'OuterPosition', OuterPosition);      

        
        print('-depsc2', [outputdir,'/',outfilename,'.eps'])
        close(figure1);
        
        
    end
    function EgeszAbra2(idx)
        
% ********************** plot the DOS ***********************

        % determine points to be plotted
        indexes = logical( density_of_states_right_electron(1,:));
        
        % creating figure in units of pixels
        figure1 = figure( 'Units', 'Pixels', 'Visible', 'on','pos',[0 0 900 400]);%, 'Colormap', ... % summer, 'pos',[10 10 900 700] );%, ...
    %[0.0416666679084301 0 0;0.0833333358168602 0 0;0.125 0 0;0.16666667163372 0 0;0.20833332836628 0 0;0.25 0 0;0.291666656732559 0 0;0.333333343267441 0 0;0.375 0 0;0.416666656732559 0 0;0.458333343267441 0 0;0.5 0 0;0.541666686534882 0 0;0.583333313465118 0 0;0.625 0 0;0.666666686534882 0 0;0.708333313465118 0 0;0.75 0 0;0.791666686534882 0 0;0.833333313465118 0 0;0.875 0 0;0.916666686534882 0 0;0.958333313465118 0 0;1 0 0;1 0.0416666679084301 0;1 0.0833333358168602 0;1 0.125 0;1 0.16666667163372 0;1 0.20833332836628 0;1 0.25 0;1 0.291666656732559 0;1 0.333333343267441 0;1 0.375 0;1 0.416666656732559 0;1 0.458333343267441 0;1 0.5 0;1 0.541666686534882 0;1 0.583333313465118 0;1 0.625 0;1 0.666666686534882 0;1 0.708333313465118 0;1 0.75 0;1 0.791666686534882 0;1 0.833333313465118 0;1 0.875 0;1 0.916666686534882 0;1 0.958333313465118 0;1 1 0;1 1 0.0625;1 1 0.125;1 1 0.1875;1 1 0.25;1 1 0.3125;1 1 0.375;1 1 0.4375;1 1 0.5;1 1 0.5625;1 1 0.625;1 1 0.6875;1 1 0.75;1 1 0.8125;1 1 0.875;1 1 0.9375;1 1 1]);
        
        % font size on the figure will be 16 points
        fontsize = 12;

        DOS_sum_right  = real(density_of_states_right_electron(:,idx)) + real(density_of_states_right_hole(:,idx));
        DOS_diff_right = real(density_of_states_right_electron(:,idx)) - real(density_of_states_right_hole(:,idx));
        DOS_sum_left  = real(density_of_states_left_electron(:,idx)) + real(density_of_states_left_hole(:,idx));
        DOS_diff_left = real(density_of_states_left_electron(:,idx)) - real(density_of_states_left_hole(:,idx));
        % define the axis limits        
        
        x_lim = [min(Evec) max(Evec)];%/min(abs(Delta));
        y_lim = [min([min(DOS_diff_right),min(DOS_diff_left)]) 1.05*max([max(DOS_sum_right),max(DOS_sum_left)])];
        
        axes_DOS_right = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...                      
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        plot(Evec,DOS_sum_right,'LineStyle','-','Color','red','Parent', axes_DOS_right)
        plot(Evec,DOS_diff_right,'LineStyle','-','Color','blue','Parent', axes_DOS_right)
        
        % Create xlabel
        xlabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_right);
        
        % Create ylabel
        ylabel('Left side','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_right);
        
        legend({'$\rho_{e}+\rho_{h}$','$\rho_{e}-\rho_{h}$'},'Interpreter','Latex','FontSize',fontsize,'Location','best');
        
        %---------------------------------------------------------------

        axes_DOS_left = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...      
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        plot(Evec,DOS_sum_left,'LineStyle','-','Color','red','Parent', axes_DOS_left)
        plot(Evec,DOS_diff_left,'LineStyle','-','Color','blue','Parent', axes_DOS_left)
        
        % Create xlabel
        xlabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left);
        
        % Create ylabel
        ylabel('Right side','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left);
        
        legend({'$\rho_{e}+\rho_{h}$','$\rho_{e}-\rho_{h}$'},'Interpreter','Latex','FontSize',fontsize,'Location','best');

        %---------------------------------------------------------------
%{
        axes_DOS_left_electron = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...      
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        density_of_states2plot = real(density_of_states_left_electron(:,indexes));
        plot(Evec,density_of_states2plot,'LineStyle','-','Color','k','Parent', axes_DOS_left_electron)
        
        
        % Create xlabel
        xlabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left_electron);
        
        % Create ylabel
        ylabel('DOS','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left_electron);
        
        
        
        %---------------------------------------------------------------

        axes_DOS_left_hole = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...      
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        density_of_states2plot = real(density_of_states_left_hole(:,indexes));
        plot(Evec,density_of_states2plot,'LineStyle','-','Color','k','Parent', axes_DOS_left_hole)
        
        % Create xlabel
        xlabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left_hole);
        
        % Create ylabel
        ylabel('DOS','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left_hole);
%}
        
        
        % setting figure position
        figure_pos = get( figure1, 'Position' );
        
        %set the position of the axes_DOS
        OuterPosition = get(axes_DOS_right, 'OuterPosition');
        OuterPosition(1) = 0;
        OuterPosition(2) = 0;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4);
        set(axes_DOS_right, 'OuterPosition', OuterPosition);
        
        %set the position of the axes_pot
        OuterPosition = get(axes_DOS_left, 'OuterPosition');
        OuterPosition(1) = figure_pos(3)/2;
        OuterPosition(2) = 0;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4);
        set(axes_DOS_left, 'OuterPosition', OuterPosition);
        
        
%{
        %set the position of the axes_DOS
        OuterPosition = get(axes_DOS_left_electron, 'OuterPosition');
        OuterPosition(1) = 0;
        OuterPosition(2) = 0;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_left_electron, 'OuterPosition', OuterPosition); 
        
        %set the position of the axes_pot
        OuterPosition = get(axes_DOS_left_hole, 'OuterPosition');
        OuterPosition(1) = figure_pos(3)/2;
        OuterPosition(2) = 0;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_left_hole, 'OuterPosition', OuterPosition);      
%}
        
        print('-dpng', [outputdir,'/',outfilename,'phi_',num2str(phivec(idx)),'.png'])
        close(figure1);
        
        
    end
    function EgeszAbraP()

% ********************** plot the DOS ***********************

        % determine points to be plotted
        indexes = logical( polarization_right(1,:));
        
        if length(phivec(indexes)) < 2
            return
        end
        
        % creating figure in units of pixels
        figure1 = figure( 'Units', 'Pixels', 'Visible', 'off', 'Colormap', summer, 'pos',[10 10 900 700] );
    %[0.0416666679084301 0 0;0.0833333358168602 0 0;0.125 0 0;0.16666667163372 0 0;0.20833332836628 0 0;0.25 0 0;0.291666656732559 0 0;0.333333343267441 0 0;0.375 0 0;0.416666656732559 0 0;0.458333343267441 0 0;0.5 0 0;0.541666686534882 0 0;0.583333313465118 0 0;0.625 0 0;0.666666686534882 0 0;0.708333313465118 0 0;0.75 0 0;0.791666686534882 0 0;0.833333313465118 0 0;0.875 0 0;0.916666686534882 0 0;0.958333313465118 0 0;1 0 0;1 0.0416666679084301 0;1 0.0833333358168602 0;1 0.125 0;1 0.16666667163372 0;1 0.20833332836628 0;1 0.25 0;1 0.291666656732559 0;1 0.333333343267441 0;1 0.375 0;1 0.416666656732559 0;1 0.458333343267441 0;1 0.5 0;1 0.541666686534882 0;1 0.583333313465118 0;1 0.625 0;1 0.666666686534882 0;1 0.708333313465118 0;1 0.75 0;1 0.791666686534882 0;1 0.833333313465118 0;1 0.875 0;1 0.916666686534882 0;1 0.958333313465118 0;1 1 0;1 1 0.0625;1 1 0.125;1 1 0.1875;1 1 0.25;1 1 0.3125;1 1 0.375;1 1 0.4375;1 1 0.5;1 1 0.5625;1 1 0.625;1 1 0.6875;1 1 0.75;1 1 0.8125;1 1 0.875;1 1 0.9375;1 1 1]);
        
        % font size on the figure will be 16 points
        fontsize = 12;        
        
        % define the colorbar limits
        %mins_tmp = [ min(min(real(polarization_right))), min(min(real(polarization_left)))];
        %maxs_tmp = [ max(max(real(polarization_right))), max(max(real(polarization_left)))];
        
        colbarlimits = [-1 ,1];
        
        % define the axis limits
        x_lim = [min(phivec) max(phivec)];
        y_lim = [min(Evec) max(Evec)]/min(abs(Delta));
        %y_lim = [min(Evec) max(Evec)]/0.1;
        
        
        axes_DOS_right = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...      
                'CLim', colbarlimits, ...
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        [X,Y] = meshgrid( phivec(indexes), Evec );
        levelnum = 50;
        density_of_states2plot = real(polarization_right(:,indexes));
        db=1;
        contourf(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_right);
       
        
        % Create xlabel
        xlabel('\Delta\theta','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_right);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_right);
        
        %---------------------------------------------------------------

        axes_DOS_left = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...      
                'CLim', colbarlimits, ...
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        [X,Y] = meshgrid( phivec(indexes), Evec );
        levelnum = 50;
        density_of_states2plot = real(polarization_left(:,indexes));
        db=1;
        contourf(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_left);
        
        % Create xlabel
        xlabel('\Delta\theta','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_left);
        
        cb=colorbar;
        cb.Position = [ 0.46 0.61 0.04 0.315 ];
        
        % setting figure position
        figure_pos = get( figure1, 'Position' );
        
        %set the position of the axes_DOS
        OuterPosition = get(axes_DOS_right, 'OuterPosition');
        OuterPosition(1) = 0;
        OuterPosition(2) = figure_pos(4)/2;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_right, 'OuterPosition', OuterPosition);  
        
        %set the position of the axes_pot
        OuterPosition = get(axes_DOS_left, 'OuterPosition');
        OuterPosition(1) = figure_pos(3)/2;
        OuterPosition(2) = figure_pos(4)/2;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_left, 'OuterPosition', OuterPosition);  
        
        print('-dpng', [outputdir,'/',outfilename,'_polarization.png'])
        close(figure1);
        
    end
% plot the scattering region 
    function ScatterPlot()

        figure1 = figure( 'Units', 'Pixels', 'Visible', 'on');

        cRibbon.Transport(1e-2, 'gfininvfromHamiltonian', true);
    
        CreateH = cRibbon.CreateH();
        coord = CreateH.Read('coordinates');
        x=coord.x;
        y=coord.y;

        plot(x,y,'x');
        xlabel('x'); ylabel('y');
        daspect([1 1 1]);

        print('-dpng', [outputdir,'/scatterplot.png'])
        close(figure1);

    end
% plot the conductance of the sample
    function ConductancePlot( flux )

        figure1 = figure( 'Units', 'Pixels', 'Visible', 'on' );

        plot(Evec,Conductance);
        ylabel('Conductance');
        xlabel('E');
        title(['Conductance of H=',num2str(height),', W=',num2str(width),' samples']);

        print('-dpng', [outputdir,'/conductance_plot_flux=',num2str(flux),'.png'])
        close(figure1);

    end
% calculate a basic spectra with the given parameters
    function Spectra( )

        figure1 = figure( 'Units', 'Pixels', 'Visible', 'on' ,'pos',[ 0 0 500 300]);

        HLead=CreateLeadHamiltonians(Opt, param, 'Hanyadik_Lead', 1, 'q',0);
        HLead.CreateHamiltonians();

        ka_vec=2*pi/3-0.1:0.2/300:2*pi/3+0.1;
        k_points=300;
        db=20;
        k_lim=int16(db*(k_points+1)/2);

        [SpectrumData2]=HLead.CalcSpektrum('ka_vec',ka_vec,'db',db,'offset',-dope,'calcWaveFnc',false);

        plot(SpectrumData2(1:k_lim*2),SpectrumData2(k_lim*2+1:k_lim*4),'.','MarkerSize', 7);%, 'color', [1 0 0]);
        %ylim([-2 2]);

        Opt.BdG = true; 

        hold on

        pause(1)

        HLead=CreateLeadHamiltonians(Opt, param, 'Hanyadik_Lead', 1, 'q',0);
        HLead.CreateHamiltonians();
        %db=30;

        [SpectrumData2]=HLead.CalcSpektrum('ka_vec',ka_vec,'db',db,'offset',-dope,'calcWaveFnc',false);

        plot(SpectrumData2(1:k_lim*2),SpectrumData2(k_lim*2+1:k_lim*4),'.','MarkerSize', 7);%, 'color', [1 0 0]);

        %xlim([2 2.188]);
        xlim([1.944 2.244]);
        %xticks([2,2.094,2.188]);
        xticks([1.944,2.094,2.244]);
        %xticklabels(["7\pi/12","2\pi/3","9\pi/12"]);
        %ylim([-1*dope-1 -1*dope+1]);
        ylim([-0.2 0.2]);
        yticks([-0.2,-0.1,0,0.1,0.2]);
        yticklabels(["-2\Delta","-\Delta","0","\Delta","2\Delta"]);
        xlabel('$k$','Interpreter','latex');
        ylabel('$E\,[\Delta]$','Interpreter','latex');
        legend('Normal graphene','Superconducting graphene');
        print('-dpng', [outputdir,'/spectra_ribbon','.png'])
        close(figure1);
    end

%% sets the output directory
    function setOutputDir()
        resultsdir = ['ABS_spectral2_H',num2str(height),'_W',num2str(width),'_flux',num2str(flux),'_Cin',num2str(Circ_in),'_Cout',num2str(Circ_out),'_EF',num2str(EF),'_LOC',location,'_res',num2str(resolution)];
        mkdir(resultsdir );
        outputdir = resultsdir;        
    end
end
