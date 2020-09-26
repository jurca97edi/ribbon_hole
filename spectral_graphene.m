function ret = spectral_graphene( height, width, flux , Circ_in, Circ_out , EF )

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
    
    height=str2num(height);
    width=str2num(width);
    flux=str2num(flux);
    Circ_in=str2num(Circ_in);
    Circ_out=str2num(Circ_out);
    EF=str2num(EF);
    
    % The Fermi level in the contacts
    %EF = 0.2e-0;%-PotentialStrength in eV
    % imaginary part of the Fermi energy
    eta = 1e-8;

	% the chosen transverse momentum
	%q = 5e-3;
    
    
    % Planck contant
    h = 6.626e-34;
    % The charge of the electron
    qe = 1.602e-19;
    % atomic distance
    rCC = 1.42*1e-10; %In Angstrom
    % flux quanta
    flux0 = h/(2*qe); 
    
    % The outfilename
    outfilename = [fncname, '_',num2str( filenum )];
    outputdir   = [];
    
    % The input and output XML files
    inputXML = 'Graphene_Input.xml';
    outputXML = [outfilename,'.xml'];
    
    % Loading the input parameters
    [Opt, param] = parseInput( inputXML );
    
    %dope the leads
    dope = 0.5;
    
    % setting the Fermi energy in the superconducting lead
    param.Leads{1}.epsilon = param.Leads{1}.epsilon - dope;
    param.Leads{2}.epsilon = param.Leads{2}.epsilon - dope;
    
    %shift all epsilons to Fermi energy EF
    param.scatter.epsilon = param.scatter.epsilon - EF;
    %param.scatter.epsilon = param.scatter.epsilon - 0.1;
    
    param.Leads{1}.M = width;
    param.Leads{2}.M = width;
    
    Delta = max( [param.Leads{1}.pair_potential, param.Leads{2}.pair_potential]);
       
    
    setOutputDir();       
       
    % class representing the two terminal ribbon structure
    cRibbon = [];
    
    % the width in the middle
    middle_width = 2*width;

    % the inner hole
    cCircle_in = structures('circle');
    cCircle_in.center.x = (width+1)/2;
    cCircle_in.center.y = (height+1)*sqrt(3)/2;
    cCircle_in.radius = Circ_in*width;
    
    % outher circle
    cCircle_out = structures('circle');
    cCircle_out.center = cCircle_in.center;
    cCircle_out.radius = Circ_out*width;
    
    % setting the energy and tranverse momentum grids for the calculations
    [Evec, phivec] = setVectors();

    
    % The calculated density of states
    density_of_states_upper_electron = zeros( length(Evec), length(phivec));
    density_of_states_upper_hole = zeros( length(Evec), length(phivec));
    density_of_states_lower_electron = zeros( length(Evec), length(phivec));
    density_of_states_lower_hole = zeros( length(Evec), length(phivec));

    % creating function handle for the Hamiltonians
    Opt.BdG = false;

    hLeadModel = @LeadModel;

    cRibbon=[];

    % creating the NTerminal class
    cRibbon = Ribbon_hole('width', width, 'height', height, 'Opt', Opt, 'param', param, 'filenameOut', fullfile( outputdir, [outfilename, '.xml']), ...
                       'leadmodel', hLeadModel, 'cCircle_in', cCircle_in, 'cCircle_out', cCircle_out, 'middle_width', middle_width); 
    %{
    Flux=[0; pi/2; pi; pi*3/2; 2*pi];

    Conductance = zeros(length(Evec),1);
               
    for jdx=1:length(Flux)

        % create handle for the piercing magnetic flux
        flux=Flux(jdx);
        CreateHandlesForMagneticField( flux );
        parfor kdx=1:length(Evec)
            Conductance(kdx) = cRibbon.Transport(Evec(kdx), 'gfininvfromHamiltonian', true);
        end
        disp([ num2str(jdx/length(Flux)*100),' % calculated of the conductance.'])
        
        ConductancePlot(flux);

    end
    %}
                       
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
    %ScatterPlot();
    %return            
    
    % calculate the spectra of a slice with the give paramteres
    %Spectra();return
    
    Opt.BdG = true;
    
    % open parallel pool (if not opened already)
    parallelmanager = Parallel( Opt );
    %parallelmanager.openPool();     
    
    % calculate the density of states
    CalculateSpectralDensity();
    
    % close parallel pool
    %parallelmanager.closePool();  
    

    save( [outputdir,'/',outfilename, '.mat'], 'Evec', 'height', 'EF', 'density_of_states_upper_electron', 'phivec', ...
                'density_of_states_upper_hole','density_of_states_lower_electron','density_of_states_lower_hole','flux' );
    
    EgeszAbra();    


%% SetVectors
%> @brief Sets the energy and transverse quantum number grid points
%> @return [1] 1D array of the energy points.
%> @return [2] 1D array of the transverse momentum points.
    function [Evec, phivec] = setVectors()
        %Delta = 1e-3;
        %Evec = (0*min(abs(Delta)):max(abs(Delta))/100:1.05*max(abs(Delta)));
        Evec = (0*min(abs(Delta)):max(abs(Delta))/50:1.05*max(abs(Delta)));
        %Evec = (0.5*min(abs(Delta)):max(abs(Delta))/100:0.65*max(abs(Delta)));
        
        %phivec = 0:pi/100:pi;
        phivec = 0:pi/50:pi;

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
            Rho_upper_electron = zeros( length(Evec), 1);
            Rho_upper_hole     = zeros( length(Evec), 1);
            Rho_lower_electron = zeros( length(Evec), 1);
            Rho_lower_hole     = zeros( length(Evec), 1);
            
            % creating the NTerminal class
            cRibbon = Ribbon_hole('width', width, 'height', height, 'Opt', Opt, 'param', param, 'filenameOut', fullfile( outputdir, [outfilename, '.xml']), ...
                     'leadmodel', hLeadModel, 'cCircle_in', cCircle_in, 'cCircle_out', cCircle_out, 'middle_width', middle_width); 
            
            % modify Hamiltonian of the scattering r.
            %CreateH = cRibbon.CreateH();
            %cRibbon.ApplyPotentialInScatter( CreateH, @ScatterPot);   
               
            % creating funcfion handles for the magnetic vector potentials
            CreateHandlesForMagneticField( flux )

            % functio handle to pick the central sites in the scattering region
            hChoseSites = @ChoseSites;           
            hScatterPot = @ScatterPot;
            
            for jdx=1:length(Evec)
            %parfor jdx=1:length(Evec)
                 Energy = Evec(jdx);                 
                
                 % create an instance of class DOS to calculate the density of states along the whole scattering region
                 % modify Hamiltonian of the scattering r.
                 DOS_handles = DOS( Opt, 'junction', cRibbon, 'useSelfEnergy', false, 'scatterPotential', hScatterPot);
                 %DOS_handles = DOS( Opt, 'junction', cRibbon, 'useSelfEnergy', false);
                 
                 %CreateH = cRibbon.CreateH();
                 %Hscatter = CreateH.Read('Hscatter');
                 
                 % calculate the local DOS along the whole scattering region
                 cLocalDOS = DOS_handles.LocalDOSCalc( Energy+1i*eta, 'ChoseSites', hChoseSites  );
                         
                 % determine the sites of the bottom layer
                 lower_sites = cLocalDOS.coordinates.x < mean(cLocalDOS.coordinates.x);              
                             
                 % determine the Density of states for the upper/lower electrons/holes
                 Rho_upper_electron(jdx) = sum(cLocalDOS.DOSvec( (~lower_sites) & cLocalDOS.coordinates.BdG_u ));                   
                 Rho_upper_hole(jdx)     = sum(cLocalDOS.DOSvec( (~lower_sites) & (~cLocalDOS.coordinates.BdG_u) ));                   
                 Rho_lower_electron(jdx) = sum(cLocalDOS.DOSvec( (lower_sites) & cLocalDOS.coordinates.BdG_u ));                   
                 Rho_lower_hole(jdx)     = sum(cLocalDOS.DOSvec( (lower_sites) & (~cLocalDOS.coordinates.BdG_u) ));   
                              
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
%                lower_sites = abs(coordinates.z) < 0.5*max(coordinates.z); 
%                
%                % determine the Density of states for the upper/lower electrons/holes
%                Rho_upper_electron(jdx) = sum(DOSvec( (~lower_sites) & coordinates.BdG_u ));                   
%                Rho_upper_hole(jdx)     = sum(DOSvec( (~lower_sites) & (~coordinates.BdG_u) ));                   
%                Rho_lower_electron(jdx) = sum(DOSvec( (lower_sites) & coordinates.BdG_u ));                   
%                Rho_lower_hole(jdx)     = sum(DOSvec( (lower_sites) & (~coordinates.BdG_u) ));   
                
            end
            
            density_of_states_upper_electron(:,idx) = Rho_upper_electron;
            density_of_states_upper_hole(:,idx)     = Rho_upper_hole;
            density_of_states_lower_electron(:,idx) = Rho_lower_electron;
            density_of_states_lower_hole(:,idx)     = Rho_lower_hole;
            
            % calculated data are plotted at each step
            EgeszAbra()
            
            % exporting the calculated data
            save( [outputdir,'/',outfilename, '.mat'], 'Evec', 'height', 'EF', 'density_of_states_upper_electron', 'phivec', ...
                'density_of_states_upper_hole','density_of_states_lower_electron','density_of_states_lower_hole','flux' );
            
            % displaying the status of the calculations
            disp([ num2str(idx/length(phivec)*100),' % calculated of the density of states.'])
        end
        
        
                
    end


%% LeadModel
%> @brief insert Bz magnetic field into the leads
    function cLead = LeadModel( lead_idx, E, varargin )
        
    p = inputParser;
    
    p.addParameter('createCore', 0);
    p.addParameter('Just_Create_Hamiltonians', 0);
    p.addParameter('shiftLead', 0);
    p.addParameter('coordinates_shift', 0);
    p.addParameter('transversepotential',[]);
    p.addParameter('Lead',[]);
    p.addParameter('gauge_field', [] );% gauge field for performing gauge transformation
    p.addParameter('SelfEnergy',false); % set true to calculate the self energy of the semi-infinite lead
    p.addParameter('SurfaceGreensFunction', true );% set true to calculate the surface Greens function of the semi-infinite lead
    p.addParameter('leadmodel', []); %function handle for an individual physical model for the contacts
    p.addParameter('CustomHamiltonian', []); 
    p.addParameter('q', [] ) %transverse momentum
    p.parse(varargin{:});
    createCore            = p.Results.createCore;
    Just_Create_Hamiltonians = p.Results.Just_Create_Hamiltonians;
    shiftLead             = p.Results.shiftLead;
    coordinates_shift     = p.Results.coordinates_shift;
    transversepotential   = p.Results.transversepotential;
    Lead_ret              = p.Results.Lead;
    gauge_field           = p.Results.gauge_field; % gauge field for performing gauge transformation
    SelfEnergy            = p.Results.SelfEnergy;
    SurfaceGreensFunction = p.Results.SurfaceGreensFunction;
    q                     = p.Results.q;
    CustomHamiltonian     = p.Results.CustomHamiltonian;
    leadmodel             = p.Results.leadmodel;
    
    
    Opt2 = Opt;
    Opt2.magnetic_field = false;
    cTransport_Interface = Transport_Interface(E, Opt2, param );
    
    %[Opt2, param] = ValidateStructures( Opt2, param );
    
    cLead = cTransport_Interface.SurfaceGreenFunctionCalculator( lead_idx, 'createCore', createCore, ...
                            'Just_Create_Hamiltonians', Just_Create_Hamiltonians, ...
                            'shiftLead', shiftLead, ...
                            'coordinates_shift', coordinates_shift, ...
                            'transversepotential', transversepotential, ...
                            'SelfEnergy', SelfEnergy, ...
                            'SurfaceGreensFunction', SurfaceGreensFunction);

    % perform gauge transformation
    %[~, ~, hgauge_field ] = createVectorPotential( -flux(jdx) );
    [~, ~, hgauge_field ] = createVectorPotential( -flux );

    % retriving the corrdinates
    coordinates = cLead.Read('coordinates');
    kulso_szabfokok = cLead.Read('kulso_szabfokok');
    
    % creating objet to perform gauge transformation
    cPeierls = Peierls(Opt);

    gsurf_inv = cLead.Read('gsurfinv');

    coordinates = coordinates.KeepSites(kulso_szabfokok);

    gsurf_inv = cPeierls.gaugeTransformation( gsurf_inv, coordinates, hgauge_field );

    % saving the gauge transformed Greens function
    
    cLead.Write('gsurfinv', gsurf_inv);
    
    end 
    function ret = ScatterPot( CreateH , Energy)
    
        coordinates = CreateH.Read('coordinates');
        x = coordinates.x;
        y = coordinates.y;
        BdG_u = coordinates.BdG_u;

        Hscatter = CreateH.Read('Hscatter');
        
        % sites on the right side of the scattering r.
        sites2shift = x >= cCircle_in.center.x; 
        
        % shift the on-site energy - take into account the BdG Ham.
        Hscatter = Hscatter + sparse(diag ( ( sites2shift & BdG_u )*2*EF - ( sites2shift & ~BdG_u )*2*EF));
        diagl= diag(Hscatter);
        %CreateH.Write('Hscatter', Hscatter);

        %figure1 = figure('rend','painters','pos',[10 10 900 400]);
  
        %plot(coordinates.x,coordinates.y,'x','MarkerSize',5,'color','k');
        %daspect([1 1 1])

        %fontsize=12;    
        %xlabel('x [a]','FontSize', fontsize,'FontName','Times New Roman');
        %ylabel('y [a]','FontSize', fontsize,'FontName','Times New Roman'); 
        
        %print('-dpng', fullfile(outputdir,['scatterplot']))
        %close(figure1);
        
        ret = zeros(1,length(x));      
        %ret = sparse(diag ( ( sites2shift & BdG_u )*2*EF - ( sites2shift & ~BdG_u )*2*EF));
    end
 

%% CreateHandlesForMagneticField
%> @brief Creates and set function handles of the magnetic vector potentials in the Ribbon class
%> @param B The magnetic field
    function CreateHandlesForMagneticField( flux )        
        [hLead, hScatter, gauge_field ] = createVectorPotential( flux );
        %cRibbon.setHandlesForMagneticField('scatter', hScatter, 'lead', hLead, 'gauge_field', gauge_field );
        cRibbon.setHandlesForMagneticField('scatter', hScatter, 'lead', hLead );
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
        
        %hScatter = @(x,y)(CircularVectorPotential(x,y, flux, cCircle_in));
        hScatter = @(x,y)(ConstantVectorPotential(x,y, flux , cCircle_in ));
        hgauge_field = @(x,y)(CircularGaugeField(x,y, flux, cCircle_in));
        %hgauge_field = @(x,y)(ConstantGaugeField(x,y));

    end

%% plotfunction
    function EgeszAbra()
        
% ********************** plot the DOS ***********************

        % determine points to be plotted
        indexes = logical( density_of_states_upper_electron(1,:));
        
        if length(phivec(indexes)) < 2
            return
        end
        
        % creating figure in units of pixels
        figure1 = figure( 'Units', 'Pixels', 'Visible', 'off', 'Colormap', ...
    [0.0416666679084301 0 0;0.0833333358168602 0 0;0.125 0 0;0.16666667163372 0 0;0.20833332836628 0 0;0.25 0 0;0.291666656732559 0 0;0.333333343267441 0 0;0.375 0 0;0.416666656732559 0 0;0.458333343267441 0 0;0.5 0 0;0.541666686534882 0 0;0.583333313465118 0 0;0.625 0 0;0.666666686534882 0 0;0.708333313465118 0 0;0.75 0 0;0.791666686534882 0 0;0.833333313465118 0 0;0.875 0 0;0.916666686534882 0 0;0.958333313465118 0 0;1 0 0;1 0.0416666679084301 0;1 0.0833333358168602 0;1 0.125 0;1 0.16666667163372 0;1 0.20833332836628 0;1 0.25 0;1 0.291666656732559 0;1 0.333333343267441 0;1 0.375 0;1 0.416666656732559 0;1 0.458333343267441 0;1 0.5 0;1 0.541666686534882 0;1 0.583333313465118 0;1 0.625 0;1 0.666666686534882 0;1 0.708333313465118 0;1 0.75 0;1 0.791666686534882 0;1 0.833333313465118 0;1 0.875 0;1 0.916666686534882 0;1 0.958333313465118 0;1 1 0;1 1 0.0625;1 1 0.125;1 1 0.1875;1 1 0.25;1 1 0.3125;1 1 0.375;1 1 0.4375;1 1 0.5;1 1 0.5625;1 1 0.625;1 1 0.6875;1 1 0.75;1 1 0.8125;1 1 0.875;1 1 0.9375;1 1 1]);
        
        % font size on the figure will be 16 points
        fontsize = 12;        
        
        % define the colorbar limits
        colbarlimits = [min(min(real(log(density_of_states_upper_electron)))) max(max(real(log(density_of_states_upper_electron))))];
        
        % define the axis limits
        x_lim = [min(phivec) max(phivec)];
        y_lim = [min(Evec) max(Evec)]/min(abs(Delta));
        
        
        axes_DOS_upper_electron = axes('Parent',figure1, ...
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
        density_of_states2plot = real(log(density_of_states_upper_electron(:,indexes)));
        db=1;
        contour(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_upper_electron);
       
        
        % Create xlabel
        xlabel('\Phi','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_upper_electron);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_upper_electron);
        
        
        
        
        %---------------------------------------------------------------

        axes_DOS_upper_hole = axes('Parent',figure1, ...
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
        density_of_states2plot = real(log(density_of_states_upper_hole(:,indexes)));
        db=1;
        contour(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_upper_hole);
        
        % Create xlabel
        xlabel('\Phi','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_upper_hole);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_upper_hole);
        
        
        
        %---------------------------------------------------------------

        axes_DOS_lower_electron = axes('Parent',figure1, ...
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
        density_of_states2plot = real(log(density_of_states_lower_electron(:,indexes)));
        db=1;
        contour(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_lower_electron);
        
        % Create xlabel
        xlabel('\Phi','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_lower_electron);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_lower_electron);
        
        
        
        %---------------------------------------------------------------

        axes_DOS_lower_hole = axes('Parent',figure1, ...
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
        density_of_states2plot = real(log(density_of_states_lower_hole(:,indexes)));
        db=1;
        contour(X(1:db:end,1:db:end), Y(1:db:end,1:db:end)/min(abs(Delta)), density_of_states2plot(1:db:end,1:db:end), levelnum ,'LineStyle','none','Fill','on',...   'LevelList',[0.2391 0.2866 0.3342 0.3817 0.4293 0.4769 0.5244 0.572 0.6195 0.6671 0.7147 0.7622 0.8098 0.8573 0.9049 0.9524],...
           'Parent', axes_DOS_lower_hole);
        
        % Create xlabel
        xlabel('\Phi','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_lower_hole);
        
        % Create ylabel
        ylabel('E [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_lower_hole);
       
        
        
        % setting figure position
        figure_pos = get( figure1, 'Position' );
        
        %set the position of the axes_DOS
        OuterPosition = get(axes_DOS_upper_electron, 'OuterPosition');
        OuterPosition(1) = 0;
        OuterPosition(2) = figure_pos(4)/2;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_upper_electron, 'OuterPosition', OuterPosition);  
        
        %set the position of the axes_pot
        OuterPosition = get(axes_DOS_upper_hole, 'OuterPosition');
        OuterPosition(1) = figure_pos(3)/2;
        OuterPosition(2) = figure_pos(4)/2;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_upper_hole, 'OuterPosition', OuterPosition);  
        
        
        
        %set the position of the axes_DOS
        OuterPosition = get(axes_DOS_lower_electron, 'OuterPosition');
        OuterPosition(1) = 0;
        OuterPosition(2) = 0;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_lower_electron, 'OuterPosition', OuterPosition); 
        
        %set the position of the axes_pot
        OuterPosition = get(axes_DOS_lower_hole, 'OuterPosition');
        OuterPosition(1) = figure_pos(3)/2;
        OuterPosition(2) = 0;
        OuterPosition(3) = figure_pos(3)/2;
        OuterPosition(4) = figure_pos(4)/2;
        set(axes_DOS_lower_hole, 'OuterPosition', OuterPosition);      

        
        print('-depsc2', [outputdir,'/',outfilename,'.eps'])
        close(figure1);
        
        
    end

% plot the scattering region 
    function ScatterPlot()

        figure1 = figure( 'Units', 'Pixels', 'Visible', 'on' );

        cRibbon.Transport(1e-2, 'gfininvfromHamiltonian', true);
    
        CreateH = cRibbon.CreateH();
        coord = CreateH.Read('coordinates');
        x=coord.x;
        y=coord.y;

        plot(x,y,'x');
        xlabel('x'); ylabel('y');
        daspect([1 1 1]);

        print('-dpng', [outputdir,'/scatterplot.png'])
        %print('-dpng', ['scatterplot.png'])
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

        HLead=CreateLeadHamiltonians(Opt, param, 'Hanyadik_Lead', 1, 'q',0);
        HLead.CreateHamiltonians();

        k_points=100;
        db=20;
        k_lim=int16(db*(k_points+1)/2);

        [SpectrumData2]=HLead.CalcSpektrum('ka_num',k_points,'db',db,'offset',0,'calcWaveFnc',false);

        plot(SpectrumData2(1:k_lim*2),SpectrumData2(k_lim*2+1:k_lim*4),'.','MarkerSize', 10);%, 'color', [1 0 0]);  
        ylim([-2 2]);

        Opt.BdG = true; 

        hold on

        pause(1)

        HLead=CreateLeadHamiltonians(Opt, param, 'Hanyadik_Lead', 1, 'q',0);
        HLead.CreateHamiltonians();

        [SpectrumData2]=HLead.CalcSpektrum('ka_num',k_points,'db',db,'offset',0,'calcWaveFnc',false);

        plot(SpectrumData2(1:k_lim*2),SpectrumData2(k_lim*2+1:k_lim*4),'.','MarkerSize', 5);%, 'color', [1 0 0]);  
 
    end

%% sets the output directory
    function setOutputDir()
        resultsdir = ['ABS_spectral_H',num2str(height),'_W',num2str(width),'_flux',num2str(flux),'_Cin',num2str(Circ_in),'_Cout',num2str(Circ_out),'_EF',num2str(EF)];
        mkdir(resultsdir );
        outputdir = resultsdir;        
    end
end
