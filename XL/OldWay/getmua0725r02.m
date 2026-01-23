
tic
jsonFile = '\\vs03\VS03-VandC-6\Neuropixels_NHP\Data_collection\m032\20240725\recording2\structure.oebin';

% list_open_ephys_binary(jsonFile, 'continuous')
% 
% list_open_ephys_binary(jsonFile, 'events')

settingsFile = '\\vs03\VS03-VandC-6\Neuropixels_NHP\Data_collection\m032\20240725\settings_2.xml';

savename = '\\vs03\VS03-VandC-6\Neuropixels_NHP\Data_collection\m032\20240725\recording2\tempmua';

saveout = 1;

[s] = xml2struct(settingsFile);

chypos = s.SETTINGS.SIGNALCHAIN{1,1}.PROCESSOR{1,1}.EDITOR.NP_PROBE.ELECTRODE_YPOS.Attributes; %!processor order in sighnal chain might change
ch = fieldnames(chypos);
ypos = zeros(length(ch),2);
for i=1:length(ch)
    ypos(i,2)=str2num(getfield(chypos,ch{i})); %yposition
    ypos(i,1)=str2num(extractAfter(ch{i},"CH")); %channel number
end
ypos_sort=sortrows(ypos,2); %channel ordered by yposition
CHorder=ypos_sort(:,1); %channel number ordered by yposition

%apfilename = 'E:\OpenEphys\m032-2024-08-26_15-21-10\Record Node 110\experiment1\recording2\continuous\Neuropix-PXI-100.ProbeA-AP\continuous.dat'


% lfp = load_open_ephys_binary(jsonFile, 'continuous', 2);
%gain 250, fs2.5
clear ap
ap = load_open_ephys_binary(jsonFile, 'continuous', 1);

%[ap, apTS, apInfo] = load_open_ephys_data_faster(apfilename);
%gain500, fs30
clear ana
ana = load_open_ephys_binary(jsonFile, 'continuous', 3);
%fs2.5
diode = ana.Data(5,:);
toc

diodeTS = ana.Timestamps;

% diodeTrialDur = find(diff(diode)>8000)'; %!check threshold for different tasks
% diodeTrialOn = 1+find(diff(diodeTrialDur)>2000); %find the first timepoint of the trial, vary by task
% % diodeTrialOn = 1+find(diff(diodeTrialDur)>250); % natural images? 
% diodeTrialOnIdx = [diodeTrialDur(1);diodeTrialDur(diodeTrialOn)];
% diodeTrialOnTS = diodeTS(diodeTrialOnIdx);

% figure;plot(diode)
% hold on
% scatter(diodeTrialOnIdx,diode(diodeTrialOnIdx)) %check if onset time is sensible


clear bits
bits = load_open_ephys_binary(jsonFile, 'events', 3);

BitsTS = bits.Timestamps;

bdata = bits.Data;
bdCo = find(bdata==6); %find correct trials
BTrialDur = [bdCo-1;bdCo-2]; %varies by task
BTrialOn = find(bdata(BTrialDur)==4); %find trial start
BTrialOnIdx = sort(BTrialDur(BTrialOn));
BitsTrialOnTS = BitsTS(BTrialOnIdx);


figure;plot(bdata)
hold on
scatter(BTrialOnIdx,bdata(BTrialOnIdx)) %check if onset time is sensible


apTimestamps = ap.Timestamps;

%find ap index closest to trial bit onset
Trls = nan(length(BitsTrialOnTS),1);
for i = 1:length(BitsTrialOnTS)
    [val,idx] = min(abs(apTimestamps-BitsTrialOnTS(i)));
    Trls(i) = idx;
end
Trials = Trls;


Fs = 30000;
trial_length = 2;
pre_trial = 0.2;
post_trial = trial_length-pre_trial;
tb =-pre_trial*1000:1:(post_trial)*1000;

%diode timestamps start before np, cut data from the start
diode_cutfront = round((apTimestamps(1) - ana.Timestamps(1))*ana.Header.sample_rate);
diode_addsample = zeros(1,round((apTimestamps(end) - ana.Timestamps(end))...
    *ana.Header.sample_rate));
diode_align = [diode(diode_cutfront+1:end),diode_addsample];
diodeU = upsample(diode_align,12);

%get diode onset for each trial
for trl=1:length(Trials)
    clear diode_trl
    temp_int = Trials(trl):Trials(trl)+Fs*trial_length-1;
    diode_trl=diodeU(temp_int);
    f=find(abs(diode_trl)>5000,1,'first'); %!threshold might vary by task
    if isempty(f)
        dio_onset(trl) = 0;
    else
        dio_onset(trl)= f;
    end
end

% diodeT = BitsTrialOnTS + diode_onset'/30000;
clear RAW
RAW = ap.Data;%/500; %adjust for gain

%common average ref
mean_RAW = int16(mean(RAW,1));
RAW_car = RAW - mean_RAW;

downsample_factor = 30;

bp_low_cutoff = 300;
bp_high_cutoff = 5000;

% band-pass filter
bpFilt = designfilt('bandpassiir', 'FilterOrder', 4, ...
                    'HalfPowerFrequency1', bp_low_cutoff, ...
                    'HalfPowerFrequency2', bp_high_cutoff, ...
                    'SampleRate', Fs);

lp_cutoff = 50;

% low-pass filter
lpFilt = designfilt('lowpassiir', 'FilterOrder', 4, ...
                    'HalfPowerFrequency', lp_cutoff, ...
                    'SampleRate', Fs);

clear temp_mua raw                
temp_mua = nan([length(CHorder) length(Trials) length(tb)]);
for chn = 1:4:length(CHorder)
    %raw=RAW(CHorder(chn)+1,:);
    raw=RAW_car(CHorder(chn)+1,:);
    raw=double(raw);
    
    raw_bpfilt = filtfilt(bpFilt, raw);
    
    rec_bpfilt = abs(raw_bpfilt);
    
    muafilt = filtfilt(lpFilt, rec_bpfilt);

    % remove 50Hz line noise (+ harmonics)
    buttLoop = muafilt;
    for i = 1:5
        d = designfilt('bandstopiir','FilterOrder',2, ...
            'HalfPowerFrequency1',50*i-5,'HalfPowerFrequency2',50*i+5, ...
            'DesignMethod','butter','SampleRate',Fs);
        buttLoop = filtfilt(d,buttLoop);
    end

    % remove 60Hz screen refresh rate (+ harmonics)
    for i = 1:5
        d = designfilt('bandstopiir','FilterOrder',2, ...
            'HalfPowerFrequency1',60*i-5,'HalfPowerFrequency2',60*i+5, ...
            'DesignMethod','butter','SampleRate',Fs);
        buttLoop = filtfilt(d,buttLoop);
    end
    muafilt = buttLoop';

    % select trials and downsample
    clear  temp_chn
    temp_chn = nan([length(Trials) length(tb)]);
    z = 1;
    for trl=1:length(Trials)
        if Trials(trl,end)>0
            clear temp_int temp_data
            temp_int = Trials(trl)+dio_onset(trl)-Fs*pre_trial:Trials(trl)+dio_onset(trl)+Fs*post_trial;
%             temp_int = Trials(trl)-Fs*pre_trial:Trials(trl)+Fs*post_trial;%-1;
            temp_chn(z,:) = downsample(muafilt(temp_int),downsample_factor);
            z=z+1;
        end
    end
    temp_mua(chn,:,:) = temp_chn;
end


% figure;
% for i=1:100
%     subplot(10,10,i)
%     plot((temp_chn(i,:)))
% end


% figure;
% for i=1:100
%     subplot(10,10,i)
%     plot((temp_chn(i+100,:)))
% end
% 
% figure;plot(temp_chn(1,:))
% 
% figure;
% for i=1:100
%     subplot(10,10,i)
%     plot(squeeze(temp_mua(311,i+100,:)))
% end



% 
% % figure;plot(mean(temp_chn))
% 
% figure;
% for i=1:100
%     subplot(10,10,i)
%     plot(squeeze(mean(temp_mua(i,:,:))))
% end
% 
% figure;
% plot(squeeze(mean(temp_mua(125,:,:))))
% 
% % figure;plot(squeeze(nanmean(temp_mua,[1,2])));
% 
% % figure;
% % for i=40:60
% %     hold on
% %     plot(squeeze(mean(temp_block_mua(i,:,:))))
% % end
% 

if saveout
save(savename,'temp_mua','-v7.3')
end

toc

% baseT = tb > -100 & tb <= 0;
% 
% normMUA = nan(size(temp_mua));
% for chn = 1:384
%     temp_chn_all = [];
%     base = nanmean(nanmean(temp_mua(chn,:,baseT),3),2);
%     temp_trial_mua = squeeze(temp_mua(chn,:,:))-base;
%     noise = nanstd(nanmean(temp_trial_mua(:,baseT),2));
%     signal = smooth(squeeze(nanmean(temp_trial_mua(:,tb>0.02))),25,'lowess');
%     temp_trial_mua = temp_trial_mua./max(signal);
%     temp_chn_all = [temp_chn_all; temp_trial_mua];
%     maxsignal = max(signal);
%     SNR(chn) = maxsignal/noise;
%     normMUA(chn,:,:) = temp_chn_all;
%     chn;
% end

% figure;
% for i=1:100
%     subplot(10,10,i)
%     plot(squeeze(mean(normMUA(i+100,:,:))))
% end