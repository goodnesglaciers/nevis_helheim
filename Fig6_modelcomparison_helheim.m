%% Figure 7 COMPARE TWO MODEL RUNS WITH AVERAGED VARIABLE VALUES
% LAS 8 June 2017
% LAS 18 MAY 2020: update for Helheim runs
clear all;
% load first model of variables saved from nevis_summary3
load varh22222_ubspatial_R227_lakerampM_4tiles_Ks100_s1e4H.mat
Q_outQ = ps.Q*[tt.Q_outQ];
Q_outq = ps.Q*[tt.Q_outq];

%% Set up figure
figure(6); clf; set(gcf,'PaperPositionMode','auto','Units','centimeters','Position',1.5.*[1 1 10 13]);
set(0,'DefaultAxesFontSize',9);

axe1(1) = axes('Position',[0.14 0.78 0.84 0.2],'Box','on','NextPlot','add','TickDir','out','XTick',[0 50 80 150 200 250 300 350],'XtickLabel',[]); 
axe1(2) = axes('Position',[0.14 0.55 0.84 0.2],'Box','on','NextPlot','add','TickDir','out','XTick',[0 50 80 150 200 250 300 350],'XtickLabel',[]);
axe1(3) = axes('Position',[0.14 0.32 0.84 0.2],'Box','on','NextPlot','add','TickDir','out','XTick',[0 50 80 150 200 250 300 350],'XtickLabel',[]); 
axe1(4) = axes('Position',[0.14 0.09 0.84 0.2],'Box','on','NextPlot','add','TickDir','out','XTick',[0 50 80 150 200 250 300 350]);
m=8;

axes(axe1(1)); xlim([0 8]); ylim([0 700])
%text(50, 6000, 'a.', 'FontSize',m+1,'FontWeight','bold');
ylabel('E and Q_{all} [m^{3} s^{-1}]');
        plot(t-228,E,'Color',[0 0.6 0],'LineWidth',1.2);
        hold on; plot(t-228,Q_out,'Color',[0.8 0 0],'LineWidth',1.2);        
        legend('E','Q_{out}');
        hold all;
        
axes(axe1(2)); xlim([0 8]); ylim([0 400])
%text(50, 700, 'a.', 'FontSize',m+1,'FontWeight','bold');
ylabel('Q Partitioned [m^{3} s^{-1}]');
        plot(t-228,Q_outQ,'Color',[0.8 0 0],'LineWidth',1.2);
        hold on; plot(t-228,Q_outq,'--','Color',[0.8 0 0],'LineWidth',1.2);        
        legend('M227 Q_{outQ} (channels)','M227 Q_{outq} (sheet)');
        hold all;
        
axes(axe1(3)); xlim([0 8]); ylim([0 1])
ylabel('Q Fraction [  ]');
        plot(t-228,Q_outQ./Q_out,'-','Color',[0.8 0 0],'LineWidth',1.2);
        hold on; plot(t-228,Q_outq./Q_out,'--','Color',[0.8 0 0],'LineWidth',1.2);        
%         legend('M227 Q_{outQ} (channels)','M227 Q_{outq} (sheet)',...
%             'M67 Q_{outQ} (channels)','M67 Q_{outq} (sheet)');
        hold all;
    
axes(axe1(4)); xlim([0 8]); %ylim([0 0.35])
%text(50, 0.3, 'b.', 'FontSize',m+1,'FontWeight','bold');
ylabel('h_{all}/A [m]');
%     plot(t,(he + hs),'Color',[0.8 0 0],'LineWidth',1.2);
%     hold on;
%     plot(t,hs ,'--','Color',[0.8 0 0],'LineWidth',1.2);
%     plot(t,he ,':','Color',[0.8 0 0],'LineWidth',1.2);
%     hold all;
    
    plot(t-228,(he + hs)./A,'Color',[0.8 0 0],'LineWidth',1);
    hold on;
    plot(t-228,hs./A ,'--','Color',[0.8 0 0],'LineWidth',1.6);
    plot(t-228,he./A ,':','Color',[0.8 0 0],'LineWidth',1.2);
    hold all;
    legend('h_{all}','h_{cav}','h_{el}','h_{all}','h_{cav}','h_{el}');
    xlabel('time [ days ]');
    
% axes(axe1(4)); xlim([0 8]); %ylim([0 8]);
% text(50, 0.01, 'c.', 'FontSize',m+1,'FontWeight','bold');
% ylabel('S_{all}/A [m]');
%         %plot(t,S,'Color',[0.8 0 0],'LineWidth',1.2); hold on;
%         %plot(t,S./A,'Color',[0.8 0 0],'LineWidth',1.2); hold on
%     
%         [AX, h1, h2] = plotyy(t-228,S./A,[1:1:365],zeros(1,365));
%         hold(AX(1)); hold(AX(2)); hold on
%         ylabel(AX(1),'S_{all}/A [m]','FontSize',m);
%         %ylabel(AX(2),'Efficient Drainage Area [%]','FontSize',m);
%         set(AX(1),'xlim',[0 8],'xtick',[0:2:300],...
%             'ylim',[0 0.01],'ytick',...
%             [0 0.005 0.01],'yticklabel',[0 0.005 0.01],'FontSize',m);
%         set(AX(2),'xlim',[0 8],'xtick',[0:2:300],...
%             'ylim',[0 80],'FontSize',m);
%         set(AX,{'ycolor'},{[0 0 0];[0 0 0]}) 
%         set(h1,'Color',[0.8 0 0],'LineWidth',1.2)
%         set(h2,'Color',[0.8 0 0],'LineWidth',1.2,'LineStyle',':')
%         
%     hold all;
    
%% load and plot second model
load varh22222_ubspatial_R67_lakerampM_4tiles_Ks100_s1e4H.mat
Q_outQ = ps.Q*[tt.Q_outQ];
Q_outq = ps.Q*[tt.Q_outq];

m=8;

axes(axe1(1)); 
hold on; plot(t-68,E,'-','Color',[0.8 0.6 0],'LineWidth',1.1);
plot(t-68,Q_out,'Color',[0 0.2 0.8],'LineWidth',1.2);
        legend('M227 E','M227 Q_{out}',...
            'M67 E','M67 Q_{out}');
        set(gca,'xlim',[0 8],'xtick',[0:1:300]);grid on; 

axes(axe1(2)); xlim([0 8]); ylim([0 400])
%text(50, 6000, 'a.', 'FontSize',m+1,'FontWeight','bold');
ylabel('Q Partitioned [m^{3} s^{-1}]');
        plot(t-68,Q_outQ,'-','Color',[0 0.2 0.8],'LineWidth',1.2);
        hold on; plot(t-68,Q_outq,'--','Color',[0 0.2 0.8],'LineWidth',1.2);        
        legend('M227 Q_{outQ} (channels)','M227 Q_{outq} (sheet)',...
            'M67 Q_{outQ} (channels)','M67 Q_{outq} (sheet)');
        set(gca,'xlim',[0 8],'xtick',[0:1:300]);grid on; hold all;
        
axes(axe1(3)); xlim([0 8]); ylim([0 1])
ylabel('Q Fraction [  ]');
        plot(t-68,Q_outQ./Q_out,'-','Color',[0 0.2 0.8],'LineWidth',1.2);
        hold on; plot(t-68,Q_outq./Q_out,'--','Color',[0 0.2 0.8],'LineWidth',1.2);        
%         legend('M227 Q_{outQ} (channels)','M227 Q_{outq} (sheet)',...
%             'M67 Q_{outQ} (channels)','M67 Q_{outq} (sheet)');
        set(gca,'xlim',[0 8],'xtick',[0:1:300]);grid on; hold all;
         
         
axes(axe1(4));     plot(t-68,(he./A) + (hs./A),'Color',[0 0.2 0.8],'LineWidth',1);
    plot(t-68,(hs./A) ,'--','Color',[0 0.2 0.8],'LineWidth',1.8);
    plot(t-68,(he./A) ,':','Color',[0 0.2 0.8],'LineWidth',1.2);
    legend('h_{all}','h_{cav}','h_{el}','h_{all}','h_{cav}','h_{el}');
    set(gca,'xlim',[0 8],'xtick',[0:1:300])
    grid on
    
% axes(axe1(4)); 
% channel_area_percentage4 = zeros(365,1);
% EDA= zeros(365,1);
%     [AX2, h3, h4] = plotyy(t-68,S./A,[1:1:365],[channel_area_percentage4]);
%     hold on; hold(AX2(1)); hold(AX2(2));
%     plot(AX2(2), [1:1:365], EDA, ':', 'Color',[0.8 0 0],'LineWidth',1.2)
%     set(AX2(1),'xlim',[0 8],'xtick',[0:2:300],...
%             'ylim',[0 0.01],'ytick',...
%             [0 0.005 0.01],'yticklabel',[0 0.005 0.01],'FontSize',m);
%         set(AX2(2),'xlim',[0 8],'xtick',[0:2:300],...
%             'ylim',[0 80],'FontSize',m);
%     set(AX2,{'ycolor'},{[0 0 0];[0 0 0]}) 
%     set(h3,'Color',[0 0.2 0.8],'LineWidth',1.2)
%     set(h4,'Color',[0 0.2 0.8],'LineWidth',1.2,'LineStyle',':')

% not averaged over area A
% axes(axe1(1)); plot(t,Q_out,'Color',[0 0.2 0.8],'LineWidth',1.2);
% axes(axe1(2));     plot(t,(he + hs) ,'Color',[0 0.2 0.8],'LineWidth',1.2);
%     plot(t,hs,'--','Color',[0 0.2 0.8],'LineWidth',1.2);
%     plot(t,he ,':','Color',[0 0.2 0.8],'LineWidth',1.2);
% axes(axe1(3)); plot(t,S,'Color',[0 0.2 0.8],'LineWidth',1.2);
% axes(axe1(4)); plot(t,N,'Color',[0 0.2 0.8],'LineWidth',1.2);

print(gcf,'-dpng','-r500',['suppfig_modelcomparison_Ks1e0_s1e4H_19May2020.png']);
