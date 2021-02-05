function average_polarization()
    
    Height=ones(1,3)*200;
    Width=ones(1,length(Height))*300;
    Flux=zeros(1,length(Height));
    EF=0.06;
    
    Cin=[0,0,0, 0.1,0.1,0.1, 0.2,0.2,0.2, 0.3,0.3,0.3, 0.4,0.4,0.4, 0.5,0.5,0.5, 0.55,0.55,0.55];
    Cin=[0.5,0.5,0.5];
    Cout=ones(1,length(Height))*0.7;
    Location=["inner","outer","center","center","outer","inner","center","outer","inner","center","outer","inner","center","outer","inner","center","outer","inner","center","outer"];
    Resolution=ones(1,length(Height))*1000;
    
    for i=1:length(Height)
        
        switch Location(i)
            case "inner"
                location = 'inner';
            case "center"
                location = 'center';
            case "outer"
                location = 'outer';
            otherwise
        end
        figure1 = figure( 'Units', 'Pixels', 'Visible', 'on' , 'pos',[10 10 800 350]);
        hold on;

        resultsdir = strjoin(['ABS_spectral_H',num2str(Height(i)),'_W',num2str(Width(i)),'_flux',num2str(Flux(i)),'_Cin',num2str(Cin(i)),'_Cout',num2str(Cout(i)),'_EF',num2str(EF),'_LOC',Location(i),'_res',num2str(Resolution(i))],"");
        %resultsdir
        load(strjoin([resultsdir,"/spectral_graphene_1.mat"],""))

        polarization_upper1 = ( density_of_states_upper_electron + density_of_states_upper_hole );
        polarization_upper2 = ( density_of_states_upper_electron - density_of_states_upper_hole );

        polarization_lower1 = ( density_of_states_lower_electron + density_of_states_lower_hole );
        polarization_lower2 = ( density_of_states_lower_electron - density_of_states_lower_hole );
        
        pos1 = [0.05 0.15 0.42 0.8];
        subplot('Position',pos1)
        
        fontsize=15;

        plot(Evec,polarization_upper1,'color','red','LineStyle','-','LineWidth',2);
        hold on;
        plot(Evec,polarization_upper2,'color','blue','LineStyle','-','LineWidth',2);
        plot([min(Evec),max(Evec)],[0,0],'k--');
        %plot([EF,EF],[Min, Max],'k-');
        title("Left branch, "+Location(i));
        legend({'$\rho_{e}+\rho_{h}$','$\rho_{e}-\rho_{h}$'},'Interpreter','Latex','FontSize',fontsize,'Location','best');
        xlim([min(Evec) max(Evec)]);%/min(abs(Delta));
        xticks([0.059,0.06,0.061]);
        xticklabels([0.59,0.6,0.61]);
        xlabel('E [\Delta]','FontSize',fontsize);
        Max = 1.05*max([abs(min([min(polarization_upper2),min(polarization_lower2)]) ),max([max(polarization_upper1),max(polarization_lower1)])]);
        ylim([-Max Max]);       
        
        pos2 = [0.55 0.15 0.42 0.8];
        subplot('Position',pos2)

        plot(Evec,polarization_lower1,'color','red','LineStyle','-','LineWidth',2);
        hold on;
        plot(Evec,polarization_lower2,'color','blue','LineStyle','-','LineWidth',2);   
        plot([min(Evec),max(Evec)],[0,0],'k--');        
        %plot([EF,EF],[Min, Max],'k-');
        title("Right branch, "+Location(i));
        legend({'$\rho_{e}+\rho_{h}$','$\rho_{e}-\rho_{h}$'},'Interpreter','Latex','FontSize',fontsize,'Location','best');
        xlim([min(Evec) max(Evec)]);%/min(abs(Delta));
        xticks([0.059,0.06,0.061]);
        xticklabels([0.59,0.6,0.61]);
        xlabel('E [\Delta]','FontSize',fontsize);
        ylim([-Max Max]);   
        
        name = ['polarization_H',num2str(Height(i)),'_W',num2str(Width(i)),'_flux',num2str(Flux(i)),'_Cin',num2str(Cin(i)),'_Cout',num2str(Cout(i)),'_EF',num2str(EF),'_LOC',location,'_res',num2str(Resolution(i)),'.png'];
        print('-dpng', name);
        close(figure1);
    end
%{
    subplot(1,2,1);    
 
    plot([-1,max(Location_index)+1],[0,0],'k--');
    plot([-1,max(Location_index)+1],[1,1],'k--');
    plot([-1,max(Location_index)+1],[-1,-1],'k--');
    
    xlim([min(Location_index)-0.5,max(Location_index)+0.5]);
    ylim([-1.2,1.2]); 
    xticks(1:length(Height)/3);
    xticklabels(string(strsplit(num2str(Cin(1:3:end))," ")));
    xlabel("Size of inner circle");
    ylabel("Electron polarization");
    title("Left branch");
    legend({'inner','center','outer'},'Location','west');

    denominator = density_of_states_upper_electron + density_of_states_upper_hole;
    denominator ( denominator == 0 ) = 1;
    polarization_upper = ( density_of_states_upper_electron - density_of_states_upper_hole )./ denominator ;

    denominator = density_of_states_lower_electron + density_of_states_lower_hole;
    denominator ( denominator == 0 ) = 1;
    polarization_lower = ( density_of_states_lower_electron - density_of_states_lower_hole )./ denominator ; 
    
    
    
    subplot(1,2,2);
    
    plot(Location_index(1:3:length(Height)),average_lower(1:3:length(Height)),'color','red','Marker','x','MarkerSize',20,'LineStyle','None','LineWidth',2)
    hold on;
    
    plot(Location_index(2:3:length(Height)),average_lower(2:3:length(Height)),'color','green','Marker','*','MarkerSize',20,'LineStyle','None','LineWidth',2)
    
    plot(Location_index(3:3:length(Height)),average_lower(3:3:length(Height)),'color','yellow','Marker','+','MarkerSize',20,'LineStyle','None','LineWidth',2)
    
    plot([-1,max(Location_index)+1],[0,0],'k--');
    plot([-1,max(Location_index)+1],[1,1],'k--');
    plot([-1,max(Location_index)+1],[-1,-1],'k--');
    
    
    xlim([min(Location_index)-0.5,max(Location_index)+0.5]);
    ylim([-1.2,1.2]);
    xticks(1:length(Height)/3);
    xticklabels(string(strsplit(num2str(Cin(1:3:end))," ")));
    xlabel("Size of inner circle");
    ylabel("Electron polarization");
    title("Right branch");
    legend({'inner','center','outer'},'Location','west');
    %}
end