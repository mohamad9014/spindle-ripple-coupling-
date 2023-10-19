function [ SIIndex, SIIndexAll, MeanTFR,T,F, MeanripplePower, Meanspindle, Numspindle,spindle ] = ripplespi_pac(fs, raw_data,data_t, IsSOs, ripple_frequency,EEG)

% Initialize PAC variable
win=0.25;% interval for computing the SI was -1 to 1 s around the SO center
win_total=3;
PadWin=win_total*2+2;

raw_data(isnan(raw_data))=0;
BlipIndex=find(IsSOs(1,:)==1);
SI=nan(5000,1);% Maximum of Slow Oscillation event is assumed 3000 event

cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'ROI';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 5:0.25:300; % analysis 5 to 30 Hz in steps of 0.25 Hz
cfg.t_ftimwin    = 5./cfg.foi;  % 5 cycles per time window
cfg.trials=1;
cfg.channel=1;
cfg.fsample=fs;
cfg.toi=[0:1/fs:2*win_total];


pp=0;
for j=1:length(BlipIndex)
    if  BlipIndex(j)-win_total*fs>0 & BlipIndex(j)+win_total*fs<length(raw_data)
         pp=pp+1;
        d=raw_data(1,BlipIndex(j)-floor(win_total*fs):BlipIndex(j)+floor(win_total*fs));
        spin=data_t(1,BlipIndex(j)-floor(win_total*fs):BlipIndex(j)+floor(win_total*fs));
        
        
        eeg = pop_importdata('data',d,'srate',fs);
        eeg.nbchan = 1;eeg.trials = 1;
        eeg.chanlocs = EEG.chanlocs(1);
        fteeg = eeglab2fieldtrip(eeg,'preprocessing');
        
        tfr = ft_freqanalysis(cfg,fteeg);
        
        powers=squeeze(tfr.powspctrm);
        time_total=tfr.time-win_total;
        
        t1 = find(time_total<=-win, 1, 'last' );
        t2 = find(time_total>=win, 1, 'first');% The interval for computing the SI was -0.25 s to	0.25 s around the SO center.
        t3 = find(time_total<=-2.5, 1, 'last' );
        t4 = find(time_total>=-1.5, 1, 'first');% TFRs were then normalized as difference to pre-event baseline (-1.5 to -1 s of the epoch)
        
        f11=ripple_frequency(1);
        f22=ripple_frequency(2);
        
        f1 = find(tfr.freq<=f11, 1, 'last' );
        f2 = find(tfr.freq>=f22, 1, 'first');
        
        
        if pp==1
            temp=nan(length(BlipIndex),length(cfg.foi),t2-t1+1);
            ripplePower=nan(length(BlipIndex),t2-t1+1);
            spindle=nan(length(BlipIndex),t2-t1+1);
            
        end
        T=time_total(t1:t2);
        F=tfr.freq;
        temp(j,:,:)=powers(:,t1:t2)-repmat(nanmean(powers(:,t3:t4),2),1,size(powers(:,t1:t2),2));
        ripplePower(j,:)=squeeze((nanmean(temp(j,f1:f2,:),2)));
        spindle(j,:)=squeeze(spin(t1:t2));
        
        %             figure;h=imagesc( time_total, tfr.freq,powers);hold on;
        %             set(gca,'YDir', 'normal');colormap jet;h=colorbar;ax=gca;ax.FontWeight ='bold';ax.FontSize =30;
        %             xlabel('Time(s)','FontSize',30,'FontWeight','Bold');ylabel('Freq(Hz)','FontSize',30,'FontWeight','Bold');
        %
        %             figure;plot(time_total(t1:t2),d(t1:t2),'LineWidth',3)
        %             ax=gca;ax.FontWeight ='bold';ax.FontSize =30;
        %             xlabel('Time(s)','FontSize',30,'FontWeight','Bold');
        %
        %             figure;plot(time_total(t1:t2),squeeze(SpindlePower(j,:)),'LineWidth',3)
        %             ax=gca;ax.FontWeight ='bold';ax.FontSize =30;
        %             xlabel('Time(s)','FontSize',30,'FontWeight','Bold');
        %
        ripplePower(j,isnan(ripplePower(j,:)))=0;
        temp_ripplePower=padarray(ripplePower(j,:),[1,PadWin*fs],0);
        ripple_Frequency=eegfilt(squeeze(temp_ripplePower(2,:)),fs,13,15,0,floor(fs/13)*3,0,'fir1');
        ripple_Frequency=ripple_Frequency(PadWin*fs:PadWin*fs+win*2*fs);
        
        temp_spin=padarray(spin(t1:t2),[1,PadWin*fs],0);
        spindle_Frequency = eegfilt(temp_spin(2,:),fs,13,15,0,floor(fs/13)*3,0,'fir1');
       spindle_Frequency =spindle_Frequency(PadWin*fs:PadWin*fs+win*2*fs);
        
        %             figure;plot(time_total(t1:t2),SO_Frequency,time_total(t1:t2),Spindle_Frequency,'LineWidth',3)
        %             ax=gca;ax.FontWeight ='bold';ax.FontSize =30;
        %             xlabel('Time(s)','FontSize',30,'FontWeight','Bold');
        %             legend('SO (Filtered 0.5-1.25)','Fast Spindle Power(Filtered 0.5-1.25)')
        %
        ph=angle(hilbert(spindle_Frequency))-angle(hilbert(ripple_Frequency));
        SI(j) =nanmean(exp(1i * ph));rad2deg(angle(SI(j)));
    end
end
SIIndex=nanmean(SI);
SIIndexAll=(SI);
Numspindle=length(BlipIndex);
MeanTFR=nanmean(temp(:,:,:),1);
MeanripplePower=nanmean(ripplePower);
Meanspindle=nanmean(spindle);
end
