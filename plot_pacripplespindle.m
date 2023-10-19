
eeg.chanlocs = EEG.chanlocs(1);
fs=1250;
data_t=Ch2;
raw_data=Ch70;
IsSOs=locs_locs;
ripple_frequency=[150 200];
[ SIIndex, SIIndexAll, MeanTFR,T,F, MeanripplePower, Meanspindle, Numspindle,spindle ] = ripplespi_pac(fs, raw_data,data_t, IsSOs, ripple_frequency,EEG)
% figure;
% subplot(2,2,1);
% h=imagesc( T, F,MeanTFR);hold on;
% plot(T,Meanspindle/4+100*ones(size(T)),T,MeanripplePower/20+100*ones(size(T)),'LineWidth',3)
% set(gca,'YDir', 'normal');colormap jet;h=colorbar;ax=gca;ax.FontWeight ='bold';ax.FontSize =30;
% xlabel('Time(s)','FontSize',30,'FontWeight','Bold');ylabel('Freq(Hz)','FontSize',30,'FontWeight','Bold');
% 
% subplot(2,2,3);
% circ_plot(angle(SIIndexAll),'hist',[],20,true,true,'linewidth',2,'color','r');hold on;
% compass(SIIndex*0.4,'r');
% 
% 




figure;
% subplot(2,2,1);
h=imagesc( T, F,squeeze(MeanTFR));hold on;
plot(T,Meanspindle/4+160*ones(size(T)),'Color',[0.2 0.6 0.2],'LineWidth',3);hold on;
set(gca,'YDir', 'normal');colormap jet;h=colorbar;ax=gca;ax.FontWeight ='bold';ax.FontSize =10;
xlabel('Time(s)','FontSize',10,'FontWeight','Bold');ylabel('Frequency(Hz)','FontSize',10,'FontWeight','Bold');
ylim([100 300]);
% legend('SO (Filtered)','Fast Spindle Power(Filtered)','FontSize',10,'FontWeight','Bold');
caxis([-40 40]);
ylim([100 250]);
 h=colorbar;h.Label.String='Baselineincrease (%)';
 h=line([0,0],[ax.YLim(1),ax.YLim(2)],'color','k','LineStyle','--','Linewidth',2)
  pos = get(gca,'Position');
 set(gca, 'position',[pos(1) pos(2)-0.01 pos(3) pos(4) ])
% subplot(2,2,3);
% figure;
% circ_plot(angle(SIIndexAll),'hist',[],20,true,true,'linewidth',2,'color','w');hold on;
% compass(SIIndex*0.4,'r');

[mu, ul, ll] = circ_mean(angle(SIIndexAll),[]);
circ_plot(angle(SIIndexAll),'hist',[],20,true,true,'linewidth',2,'color','w');
hold on;
circ_plot(mu,'hist',[],1,true,true,'linewidth',2,'color','r');
title('phase ripple spindle','Fontsize',16)


% % % 
% % % figure;
% % % Color_axis=[.7 .7 .7];r = 0.5;
% % %  c = [0 0];pos = [c-r 2*r 2*r];h=rectangle('Position',pos,'Curvature',[1 1]);hold on;h.EdgeColor=Color_axis;
% % % h=plot([-r r], [0,0]);h.Color=Color_axis;h=plot([0 0], [-r,r]);h.Color=Color_axis;
% % % r2 = r/2;pos = [c-r2 2*r2 2*r2];h=rectangle('Position',pos,'Curvature',[1 1]);hold on;h.EdgeColor=Color_axis;
% % %  axis equal;xlim([-r r]);ylim([-r r]);axis off
% % %   text(1.02, 0.5, '0','FontSize',12,'FontWeight','Bold','Unit','Normalized');
% % %   text(0.44, 1.06, '90','FontSize',12,'FontWeight','Bold','Unit','Normalized')
% % %    text(-0.29, 0.5, '±180','FontSize',12,'FontWeight','Bold','Unit','Normalized');
% % %     text(0.44, -0.05, '-90','FontSize',12,'FontWeight','Bold','Unit','Normalized');
% % %     text(0.7, 0.7, num2str(r/4),'FontSize',12,'FontWeight','Bold','Unit','Normalized');
% % %     text(.87, .87, num2str(r/2),'FontSize',12,'FontWeight','Bold','Unit','Normalized');
% % %     hold on;
% % % x=nanmean(SIIndexAll);
% % %  h=plot([0 real(nanmean(x))], [0, imag(nanmean(x))]);
% % %  h.LineWidth=3;
% % %  h.Color=[0 0 0];
% % %  hold on;
% % %  
% % %  
% % % % dumvar=0;
% % % % save(temp,'dumvar','-v7.3');
% % % % save(temp,'A', '-append');
% % % % save('temp.mat','iceChart','-v7.3');
% % % % 