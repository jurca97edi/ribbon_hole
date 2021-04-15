function diffcond_graphene( height, width, flux , Circ_in, Circ_in2 , EF , resolution)

if ~exist('filenum', 'var')
    filenum = 1;
end

filename = mfilename('fullpath');
[directory, fncname] = fileparts( filename );

% The outfilename
outfilename = [fncname, '_',num2str( filenum )];
%outputdir   = ['results_diffcond'];

% geometry
if ~exist('height', 'var')
    height = 40;% 0.25 micron
end

if ~exist('width', 'var')
    width = 10;% 0.25 micron
end

if ~exist('Circ_in', 'var')
    Circ_in = 0.7;
end

if ~exist('Circ_in', 'var')
    Circ_in2 = 0.7;
end

if ~exist('EF', 'var')
    EF = 0.2;
end

if ~exist('resolution', 'var')
    resolution = '80';
end


height=str2num(height);
width=str2num(width);
Circ_in=str2num(Circ_in);
Circ_in2=str2num(Circ_in2);
EF=str2num(EF);
resolution=str2num(resolution);

% Right side lead potential
top_gate = 0;

%set Delta to 10 meV
pair_potential = 1e-2;

% Left side lead
%bottom_gate_vec = -1e-3:2*1e-3/(resolution-1):1e-3;
%bottom_gate_vec = -1.5e-3:1.5e-3/150:1.5e-3;
bottom_gate = 0;

% The bias vector
mu_min = -pair_potential*1.05;
mu_max = pair_potential*1.05;
mu_db = resolution;

mu_vec = mu_min:(mu_max-mu_min)/(mu_db-1):mu_max;

setOutputDir()

% The phase of the superconducting conatct
%DeltaPhi =[0, pi/2, pi];
DeltaPhi =[pi, pi/2, 0];

diffcond = zeros( 3, length(DeltaPhi), length(mu_vec) ); %diffcond contains total diffcond, electron part and hole part in (~,:,~) indexes

tic

for idx = 1:length(DeltaPhi)
     
     diffcond_summed = DiffCond_ribbon_hole( DeltaPhi(idx), height, width, Circ_in , Circ_in2, EF, resolution, bottom_gate, top_gate, mu_vec, outputdir);
     
     diffcond(:,idx, :) = diffcond_summed;
     
     EgeszAbra();

     save( [outputdir,'/',outfilename, '.mat'], 'diffcond', 'mu_vec', 'EF', 'Circ_in', 'Circ_in2', 'height', 'bottom_gate', 'mu_vec' );
     
     % displaying the status of the calculations
    disp([ num2str(idx/length(DeltaPhi)*100),' % calculated of the density of states.'])

end

toc

%% plotfunction
    function EgeszAbra()
         
        diffcond_electron(:) = diffcond(2,idx,:);
        diffcond_hole(:) = diffcond(3,idx,:);
        
% ********************** plot the DOS ***********************

        % determine points to be plotted
        %indexes = logical( diffcond_to_plot);
        
        %if length(bottom_gate_vec(indexes)) < 2
        %    return
        %end
        
        % creating figure in units of pixels
        figure1 = figure( 'Units', 'Pixels', 'Visible', 'on');%, 'Colormap',...
        %[0.0416666679084301 0 0;0.0833333358168602 0 0;0.125 0 0;0.16666667163372 0 0;0.20833332836628 0 0;0.25 0 0;0.291666656732559 0 0;0.333333343267441 0 0;0.375 0 0;0.416666656732559 0 0;0.458333343267441 0 0;0.5 0 0;0.541666686534882 0 0;0.583333313465118 0 0;0.625 0 0;0.666666686534882 0 0;0.708333313465118 0 0;0.75 0 0;0.791666686534882 0 0;0.833333313465118 0 0;0.875 0 0;0.916666686534882 0 0;0.958333313465118 0 0;1 0 0;1 0.0416666679084301 0;1 0.0833333358168602 0;1 0.125 0;1 0.16666667163372 0;1 0.20833332836628 0;1 0.25 0;1 0.291666656732559 0;1 0.333333343267441 0;1 0.375 0;1 0.416666656732559 0;1 0.458333343267441 0;1 0.5 0;1 0.541666686534882 0;1 0.583333313465118 0;1 0.625 0;1 0.666666686534882 0;1 0.708333313465118 0;1 0.75 0;1 0.791666686534882 0;1 0.833333313465118 0;1 0.875 0;1 0.916666686534882 0;1 0.958333313465118 0;1 1 0;1 1 0.0625;1 1 0.125;1 1 0.1875;1 1 0.25;1 1 0.3125;1 1 0.375;1 1 0.4375;1 1 0.5;1 1 0.5625;1 1 0.625;1 1 0.6875;1 1 0.75;1 1 0.8125;1 1 0.875;1 1 0.9375;1 1 1]);
        
        % font size on the figure will be 16 points
        fontsize = 12;        
        
        % define the colorbar limits
        %colbarlimits = [min(min(diffcond_to_plot)) max(max(diffcond_to_plot))];
        
        % define the axis limits        
        Max=max([max(abs(diffcond_electron)),max(abs(diffcond_hole))]);
        x_lim = [min(mu_vec) max(mu_vec)]/pair_potential;
        y_lim = [-1.05*Max 1.05*Max];        
        
        axes_DOS_upper_electron = axes('Parent',figure1, ...
                'Visible', 'on',...
                'FontSize', fontsize,... 
                'xlim', x_lim,...           
                'ylim', y_lim,...
                'Box', 'on',...
                'Units', 'Pixels', ...
                'FontName','Times New Roman');
        hold on; 
        
        % plot the data
        %[X,Y] = meshgrid( bottom_gate_vec(indexes), mu_vec );
        X = mu_vec;
%       diffcond2plot = diffcond_to_plot';
        % plot the data
        plot(X/pair_potential,diffcond_electron,'LineStyle','-','Color','red','Parent', axes_DOS_upper_electron)
        hold on;
        plot(X/pair_potential,diffcond_hole,'LineStyle','-','Color','black','Parent', axes_DOS_upper_electron)        
        plot([x_lim(1),x_lim(2)],[0,0],'k--');
        
        legend({'EE','EH'},'Interpreter','Latex','FontSize',fontsize,'Location','best');        
        
        % Create xlabel
        xlabel('eV [\Delta]','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_upper_electron);
        
        % Create ylabel
        ylabel('dI/dV','FontSize', fontsize,'FontName','Times New Roman', 'Parent', axes_DOS_upper_electron);    
        
        print('-dpng', [outputdir,'/',outfilename,'dPhi',num2str(DeltaPhi(idx)),'.png'])
        close(figure1);
       
    end
    function setOutputDir()
        resultsdir = ['Diffcond_H',num2str(height),'_W',num2str(width),'_Cin',num2str(Circ_in),'_Cin2',num2str(Circ_in2),'_EF',num2str(EF),'_res',num2str(resolution)];
        mkdir(resultsdir );
        outputdir = resultsdir;        
    end
end
