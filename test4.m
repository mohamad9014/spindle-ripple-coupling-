%%

load ch2.mat
load ch70.mat
load g.EEGlabLOCATIONfile.mat
%%
% load states.mat
states.rem = [3905.5	4005.4
4651.5	4782.5
5382.8	5485.8
6039.8	6112.4
6700.3	6866.2
11430	11572
12637	12651
12930	13031
13804	13926
14512	14520
14762	14777
14903	14935];
states.nrem =[823.3	1015.1
1327.6	1559.7
3415.9	3899.3
4018.1	4650.6
4833.1	5379.7
5513.2	5527.3
5527.3	6038.1
6133.3	6688.2
9880.2	11424
11911	12632
12672	12818
12839	12927
13042	13791
13971	14761
14783	14900];
states.wake=[6881.3	9070.2];
fs=1250;
frequenc_filter=[7 16];
state=zeros(length(ch2),1);
s_nrem=round((states.nrem)*fs);
s_rem=round((states.rem)*fs);
s_wake=round((states.wake)*fs);
%%
for i=1:size(states.nrem,1)
    state(s_nrem(i,1): s_nrem(i,2))=3;
end
for i=1:size(states.rem,1)
    state(s_rem(i,1): s_rem(i,2))=2;
end
for i=1:size(states.wake,1)
    state(s_wake(i,1): s_wake(i,2))=1;
end
% ch2(state == 0)=[];
% ch70(state == 0)=[];
% state(state == 0)=[];
%%
f_ch2_test = eegfilt(double(ch2).', fs, 0.5, 300, 0, fix(fs/0.5)*3, 0, 'fir1');
ch2 = zscore(f_ch2_test);

f_ch70_test = eegfilt(double(ch70).', fs, 0.5, 300, 0, fix(fs/0.5)*3, 0, 'fir1');
ch70 = zscore(f_ch70_test);
%%
time=(1:length(ch2))/fs;
ch2_nrem_data=ch2(state==3);
ch70_nrem_data=ch70(state==3);

nrem_time=time(state==3);
nrem_time_con = (1:length(ch2_nrem_data))/fs;

% nrem_onset=time(find(diff(state==3)==1)+1);
% rem_ind=find(diff(state==2)==1)+1;
% rem_time=time(rem_ind);
% wake_ind=find(diff(state==1)==1)+1;
% wake_time=time(wake_ind);
%%
f_data = eegfilt(ch2_nrem_data, fs, frequenc_filter(1), frequenc_filter(2), 0, fix(fs/frequenc_filter(1))*3, 0, 'fir1');
filtered=[nrem_time_con.' f_data.'];
spindles2=FindSpindles(filtered,'threshold',3,'durations',[500 3000]);

%%

spindle_event = zeros(length(spindles2), 3);
spindle_event(:, 1) = round(spindles2(:, 1) * fs);
spindle_event(:, 3) = round(spindles2(:, 3) * fs);

for i = 1:length(spindle_event)
    [~, max_ind] = max(f_data(spindle_event(i, 1):spindle_event(i, 3)));
    spindle_event(i, 2) = max_ind + spindle_event(i, 1) - 1; 
end

%%
figure
plot(nrem_time, ch2_nrem_data, 'k');
hold on
plot(nrem_time, f_data + 5, 'b');
for i = 1:length(spindle_event)
    plot(nrem_time([spindle_event(i,1) spindle_event(i,3)]), [6 6], 'r');
    plot(nrem_time(spindle_event(i,2)), 6,'*r');
    xlabel('TIME (s)')
    yticks([0 5 6])
yticklabels(["AD" "filtered [7 16] HZ" "spindle"])
    g=gca;
    g.FontSize=14;
    
end

%%

nrem_before_rem = [3899.3 4650.6
5379.7 5527.3 
6038.1 6688.2
11424 12632
12818 12927
13791 14761
14900 14902];
spindle_ind = [];
for i = 1:length(nrem_before_rem)
    stop_N2R = nrem_before_rem(i);
    start_N2R = stop_N2R - 25;

    start_ind_N2R = find(nrem_time >= start_N2R, 1, "first");
    stop_ind_N2R = find(nrem_time <= stop_N2R, 1, "last");

    spindle_ind = [spindle_ind; find(spindle_event(:, 1) >= start_ind_N2R & spindle_event(:, 3) <= stop_ind_N2R)];
end
%%

% nrem_before_wake = [3899.3 4650.6
% 5379.7 5527.3 
% 6038.1 6688.2
% 11424 12632
% 12818 12927
% 13791 14761
% 14900 14902];
% spindle_ind_2 = [];
% for j = 1:length(nrem_before_wake)
%     stop_N2W = nrem_before_wake(j);
%     start_N2W = stop_N2W - 25;
% 
%     start_ind_N2W = find(nrem_time >= start_N2W, 1, "first");
%     stop_ind_N2W = find(nrem_time <= stop_N2W, 1, "last");
% 
%     spindle_ind_2 = [spindle_ind_2; find(spindle_event(:, 1) >= start_ind_N2W & spindle_event(:, 3) <= stop_ind_N2W)];
% end
% 


%% mean

window = 1*fs;

data_mean = zeros(1, 2*window+1);
time_mean = ((0:length(data_mean)-1)-window)/fs;

bad_ind = find(spindle_event(:,2)-window < 0);
if ~isempty(bad_ind)
    spindle_ind(spindle_ind == bad_ind) = [];
end

bad_ind = find(spindle_event(:,2)+window > length(ch2_nrem_data));
if ~isempty(bad_ind)
    spindle_ind(spindle_ind == bad_ind) = [];
end

for i = 1:length(spindle_ind)
    data_mean = data_mean + f_data(spindle_event(spindle_ind(i),2)-window:spindle_event(spindle_ind(i),2)+window) / length(spindle_ind);
end

figure
plot(time_mean, data_mean)
xlim([time_mean(1) time_mean(end)])
xlabel("Time (s)")
ylabel("Mean")

%% PAC

pac = PAC(0.25, [-2.5 -1.5], false, true);
pac.Run(ch70_nrem_data, [150 200], ch2_nrem_data, [7 16], spindle_event(spindle_ind(:),2), fs);

%% plot pac

figure
imagesc(pac.Time, pac.Freq, pac.MeanTFR*100)
set(gca,'YDir', 'normal');
colormap("jet");
cb = colorbar;
ylabel(cb,'Baseline increase(%)');
hold on
plot(pac.Time, pac.MeanData/max(abs(pac.MeanData))*(pac.Freq(end)-pac.Freq(1))/2+mean(pac.Freq), "LineWidth", 2, "Color", "k")

figure
phases = angle(pac.SI);
mu = circ_mean(phases,[]);
[pval, z] = circ_rtest(phases);

polarhistogram(phases, 20, "FaceColor", "r", "EdgeColor", "k", "Normalization", "probability");
hold on;
polarplot([mu mu], rlim, "LineWidth", 2, "Color", "b")

figure
histogram(abs(pac.SI), 20, "Normalization", "probability");
med = median(abs(pac.SI));


%%
% figure('Name','Raw data and filtered data')
% subplot(2,1,1);
% plot(nrem_time,ch2_nrem_data,'linewidth',2)
% xlim([0 10]);
% xlabel('time(s)')
% ylabel('data')
% g=gca;
% g.FontSize=14;
% subplot(2,1,2);
% plot(nrem_time,f_data,'linewidth',2)
% xlim([0 10]);
% xlabel('time(s)')
% ylabel('data')
% g=gca;
% g.FontSize=14;
% figure
% plot(time,state,'linewidth',2)
% hold on
% plot(spindles2(:,2),ones(length(spindles2(:,2)),1)*4,'*k')
% xlabel('time(s)')
% yticks([1 2 3 4])
% yticklabels(["WAKE" "REM" "NREM" "Peak spindle"])
% ylim([-1 7])
% g=gca;
% g.FontSize=14;
% rate_N2R=zeros(1,length(rem_time));
% rate_N2R(1)=(sum(spindles2(:,2)<rem_time(1))/rem_time(1))*60;
% for i=2:size(rem_time.',1)
%     N2R=sum(spindles2(:,2)<rem_time(i) & spindles2(:,2)>rem_time(i-1));
%     start=nrem_onset(nrem_onset>rem_time(i-1) & nrem_onset<rem_time(i));
%     stop=rem_time(i);
%     duration=stop-start;
%     rate_N2R(i)=(N2R/duration)*60;
% end
% ave_N2R=(sum(rate_N2R))/25;
% N2W=sum(spindles2(:,2)<wake_time(1) & spindles2(:,2)<nrem_onset(10));
% duration2=wake_time(1);
% rate_N2W=(N2W/duration2)*60;
% %% 25
% locs_locs = zeros(length(ch2), 1);
% 
% for j=1:size(rem_time.',1)
% 
%     sp_ind = spindles2(:,2)>rem_time(j)-25 & spindles2(:,2)<rem_time(j);
% 
%     locs_locs(round(spindles2(sp_ind,2)*fs)) = 1;
%     
% end
% 
% figure
% plot(locs_locs)
% 
% 
% eeg.chanlocs = EEG.chanlocs(1);
% data_t=ch2_nrem_data;
% raw_data=ch70_nrem_data;
% IsSOs=locs_locs.';
% ripple_frequency=[150 200];
% [ SIIndex, SIIndexAll, MeanTFR,T,F, MeanripplePower, Meanspindle, Numspindle,spindle ] = ripplespi_pac(fs, raw_data,data_t, IsSOs, ripple_frequency,EEG);
% 
% %%
% 
% 
% figure;
% % subplot(2,2,1);
% h=imagesc( T, F,squeeze(MeanTFR));hold on;
% plot(T,Meanspindle*10+160*ones(size(T)),'Color',[0.2 0.6 0.2],'LineWidth',3);hold on;
% set(gca,'YDir', 'normal');colormap jet;h=colorbar;ax=gca;ax.FontWeight ='bold';ax.FontSize =10;
% xlabel('Time(s)','FontSize',10,'FontWeight','Bold');ylabel('Frequency(Hz)','FontSize',10,'FontWeight','Bold');
% ylim([100 300]);
% % legend('SO (Filtered)','Fast Spindle Power(Filtered)','FontSize',10,'FontWeight','Bold');
% caxis([-50 50]);
% ylim([100 250]);
%  h=colorbar;h.Label.String='Baselineincrease (%)';
%  h=line([0,0],[ax.YLim(1),ax.YLim(2)],'color','k','LineStyle','--','Linewidth',2)
%   pos = get(gca,'Position');
%  set(gca, 'position',[pos(1) pos(2)-0.01 pos(3) pos(4) ])
% % subplot(2,2,3);
% % figure;
% % circ_plot(angle(SIIndexAll),'hist',[],20,true,true,'linewidth',2,'color','w');hold on;
% % compass(SIIndex*0.4,'r');
% 
% 
