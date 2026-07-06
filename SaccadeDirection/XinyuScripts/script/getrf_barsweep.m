uselocal = 1; % if not, use data on server, log files always on server
getmua = 1;
tic;

Monkey = 'm032';
Day = '20251016'; % YYYYMMDD
ExpN = 1;
RecN = 1;
Dataloc = '\\vs03\VS03-VandC-6\Neuropixels_NHP\Data_collection\';
%Dataloc = '\\vs03.herseninstituut.knaw.nl\VS03-VandC-6\Neuropixels_NHP\Data_collection\';
localdayDir = 'C:\OpenEphys\m032-2025-10-16_14-24-46\Record Node 110'; %change day and timestamp

figuredir = [Dataloc,Monkey,'\',Day,'\recording',num2str(RecN),'\']; %Location for the figures

Logdir = [Dataloc,Monkey,'\',Day,'\run-00',num2str(RecN),'*\'];
logfile = dir([Logdir,'sub-*.mat']);
load([logfile.folder,'\',logfile.name],'Log','Stm','Par');

if uselocal
    jsonFile = [localdayDir,'\experiment',num2str(ExpN),'\recording',num2str(RecN),'\structure.oebin'];
    settingsFile = [localdayDir,'\settings.xml'];
    ttlDir = [localdayDir,'\experiment',num2str(ExpN),'\recording',num2str(RecN),'\events\NI-DAQmx-114.PXIe-6341\TTL\'];
    apDir = [localdayDir,'\experiment',num2str(ExpN),'\recording',num2str(RecN),'\continuous\Neuropix-PXI-100.ProbeA-AP'];
    muasave = [localdayDir,'\experiment',num2str(ExpN),'\recording',num2str(RecN),'\mua',num2str(RecN)];
    savename = [localdayDir,'\experiment',num2str(ExpN),'\recording',num2str(RecN),'\bar',num2str(RecN)];
else
    jsonFile = [Dataloc,Monkey,'\',Day,'\experiment',num2str(ExpN),'\recording',num2str(RecN),'\structure.oebin'];
    settingsFile = [Dataloc,Monkey,'\',Day,'\experiment',num2str(ExpN),'\settings.xml'];
    ttlDir = [Dataloc,Monkey,'\',Day,'\experiment',num2str(ExpN),'\recording',num2str(RecN),'\events\NI-DAQmx-114.PXIe-6341\TTL\'];
    apDir = [Dataloc,Monkey,'\',Day,'\experiment',num2str(ExpN),'\recording',num2str(RecN),'\continuous\Neuropix-PXI-100.ProbeA-AP'];
    muasave = [Dataloc,Monkey,'\',Day,'\experiment',num2str(ExpN),'\mua',num2str(RecN)];
    savename = [Dataloc,Monkey,'\',Day,'\experiment',num2str(ExpN),'\bar',num2str(RecN)]; %Name of data file
end


info = jsondecode(fileread(jsonFile));
info = info.continuous;
info = info(1);
disp(['Reading ',info.folder_name])
numChannels= info.num_channels; 

saveout = 1;

%% get ttl
lines = readNPY(fullfile(ttlDir, 'states.npy'));
sampleNumbers = readNPY(fullfile(ttlDir, 'sample_numbers.npy'));
timestamps = readNPY(fullfile(ttlDir, 'timestamps.npy'));
numEvents = length(lines);
ttls = DataFrame(lines, sampleNumbers, timestamps, lines > 0, ...
    'VariableNames', {'line','sample_number','timestamp', 'state'});


%% get trial information


bdReward = find(ttls.line==8); %reward bit trial, contains extra streak rewards
%bdError = find(ttls.line==10); %error bit trial

% rewardtrial = find(strcmp({Log.events.type},'reward_start'));
% for i=1:length(rewardtrial) %find streak reward
%     streakreward(i) = strcmp({Log.events(rewardtrial(i)-1).type},'reward_stop');
% end
% 
% bdReward(streakreward) = []; %remove streak rewards from reward bit summary

% bdAll = [bdReward;bdError];
% bdAll = sort(bdAll);
bdAll = bdReward;

bdAll = bdAll(2:202); %!!!!

correcttrial = find(strcmp({Log.events.type},'correct'));
rewardtrial = find(strcmp({Log.events.type},'reward_start'));

assert(length(bdAll)==length(correcttrial),'There are (%d) recorded trials while (%d) trials were logged',length(bdAll),length(correcttrial));

% BTrialDur = [bdAll;bdAll-2;bdAll-3;bdAll-4]; %crt
BTrialDur = [bdAll;bdAll-3;bdAll-4;bdAll-5;bdAll-6]; %barsweep
BTrialOn = find(ttls.line(BTrialDur)==4); %find trials onset
BTrialOnIdx = sort(BTrialDur(BTrialOn));
BitsTrialOnTS = ttls.timestamp(BTrialOnIdx); %trial onset timestamps

nBitsTrialOn = numel(BitsTrialOnTS);


% %%ttl timestamp not saved, infer from sample number
% apTimestamps = data.timestamps;
% ttls.sample_number = double(ttls.sample_number);
% ttls.atime = ttls.sample_number/2500;
% timedifference = ttls.atime(1) - apTimestamps(1); %correct for ttl and ap onset difference
% ttls.atime = ttls.atime - timedifference;
% BitsTrialOnTS = ttls.atime(BTrialOnIdx);

%% get AP
buf = memmapfile(fullfile(apDir, 'continuous.dat'), 'Format', 'int16');
data.samples = reshape(buf.Data, [numChannels, length(buf.Data)/numChannels]);
data.sampleNumbers = readNPY(fullfile(apDir, 'sample_numbers.npy'));
data.timestamps = readNPY(fullfile(apDir, 'timestamps.npy'));
Fs = info.sample_rate;

%%
apTimestamps = data.timestamps;

%find ap index closest to stim onset
Trls = nan(length(BitsTrialOnTS),1);
for i = 1:length(BitsTrialOnTS)
    [val,idx] = min(abs(apTimestamps-BitsTrialOnTS(i)));
    Trls(i) = idx;
end
Trials = Trls;


%% get channel information    
[s] = xml2struct(settingsFile);

chypos = s.SETTINGS.SIGNALCHAIN{1,1}.PROCESSOR{1,1}.EDITOR.CUSTOM_PARAMETERS.NP_PROBE.ELECTRODE_YPOS.Attributes; %!processor order in sighnal chain might change
ch = fieldnames(chypos);
ypos = zeros(length(ch),2);
for i=1:length(ch)
    ypos(i,2)=str2num(getfield(chypos,ch{i})); %yposition
    ypos(i,1)=str2num(extractAfter(ch{i},"CH")); %channel number
end
ypos_sort=sortrows(ypos,2); %channel ordered by yposition
CHorder=ypos_sort(:,1); %channel number ordered by yposition

%%
[b_hp,a_hp] = butter(3,300/(Fs/2),'high'); 
[b_lp,a_lp] = butter(3,5000/(Fs/2),'low'); %remove high frequency noise
[b_mualp,a_mualp] = butter(3,200/(Fs/2),'low');
[b_ds,a_ds] = butter(3,0.01/1,'high'); %Spatial filter (currently tuned by eye)

%ADC selection
%Samples that are obtained together in the interleave will have common
%noise, so e.g. [1,2,25,26,49,50,...361,362];
%but we dont reref neighbouring channles to avoid subtracting spikes
%There seems to be disagreement about which channels are sampled together
%IBL says that there are 12 ADCs and each samples 32 channels e.g.
%[1,13,..etc]
%Others say also 32 channels, but grouped e.g. [1,2,25,26,etc]
z = 0;
a = 1:24:384;
b = 2:24:384;
cycle = NaN(1,384);
for s = 1:12
    cycle(a+(s-1)*2) = s;
    cycle(b+(s-1)*2) = s;
end

%Convert cycle to Phase shift per channel, relative to channel 1.
%1/13 because the LFP channel also gets sampled
ph = (rem(cycle-1,12)./13).*(1/Fs);

%%

trial_length = 2; %barsweep duration 1.8s
pre_trial = 0.2;
post_trial = trial_length - pre_trial;
pre_trialstart = round(pre_trial.*Fs); %number of samples to extract before the TTL
post_trialstart = round(post_trial.*Fs); %Number of samples to extract after the TTL

if getmua
    tmbs = -pre_trialstart:1:post_trialstart;
    tbi = tmbs./Fs;
    samps_per_trial = numel(tmbs);
    downs = 30;
    downsamps= length(1:downs:samps_per_trial);
    tbds = tbi(1:downs:end);
    MUA = zeros(downsamps,384,nBitsTrialOn);
    first_sample = data.sampleNumbers(1);
    L = samps_per_trial/Fs; %Length in seconds
    smps = 0:1:(samps_per_trial-1); %Sample numbers
    f = smps/L; %frequencies

    %Vectorize
    fmat = repmat(f',1,384);
    phmat = repmat(ph,samps_per_trial,1);

    for k = 1:nBitsTrialOn
        %tic
        samp_st = Trials(k)-pre_trialstart+1;
        samp_ed = Trials(k)+post_trialstart+1;
        cutdata = double(data.samples(1:384,samp_st:samp_ed))'*0.195;

        %% HP filter
        cutdata = filtfilt(b_hp,a_hp,cutdata);

        %% Phase align
        fbuf=fft(cutdata);
        fbuf = exp(-1j.*2.*pi.*fmat.*phmat).*fbuf;
        aligndata=ifft(fbuf,'symmetric');

        %% Low-pass at 5000hz
        aligndata = filtfilt(b_lp,a_lp,aligndata);

        %% Desripe chunck
        destriped = filtfilt(b_ds,a_ds,aligndata')';

        %% Convert to MUA
        buf = filtfilt(b_mualp,a_mualp,abs(destriped));
        for ch = 1:384
            MUA(:,ch,k) = decimate(buf(:,ch),30);
        end
    end

    save(muasave,'MUA','-v7.3'); % save mua file, trial duration x channel x trial
    disp(['MUA processed in ',num2str(toc),' seconds'])
else
    load(muasave)
end
    
%%
analyzedata = 1;

MUA = permute(MUA,[2,3,1]); % reshape mua into channel x trial x trial duration
MUA = MUA(CHorder+1,:,:);

% To plot dots on top of data:
%Useful for testing potential positions for stimuli
%Set to empty e.g. RFx = []; if you don;t want this
RFx = []; %
RFy = []; % 

fig = 1; % plot rf fit
notch = 0;

%split channels into sensible banks, or not
bank = zeros(1,384);
bank(1:48) = 1;
bank(49:96) = 2;
bank(97:144) = 3;
bank(145:192) = 4;
bank(193:240) = 5;
bank(241:288) = 6;
bank(289:336) = 7;
bank(337:384) = 8;
bankcol = [1 0 0;0 1 0;0 0 1;0 0 0;1 1 0;0 1 1;1 0 1;1 0.5 0];
% bank(1:150) = 1;
% bank(151:384) = 2;
% bankcol = [1 0 0;0 1 0];

speed = Stm.BarSpeed;
bardur = Stm.BarDur;
bardist = Stm.BarDist;
xo = Stm.RFx;
yo = Stm.RFy;
pixperdeg = Par.PixPerDeg;

trlcidx = find(strcmp({Log.events.type},'correct'));
trldir = [Log.events(trlcidx).dir];

%0 = left to right, 90 = down-to-up, 180 = right to left, 270 = up to down
trl1 = find(trldir==0);
trl2 = find(trldir==90);
trl3 = find(trldir==180);
trl4 = find(trldir==270);
trlid = trldir;
trlid(trl1) = 1;
trlid(trl2) = 2;
trlid(trl3) = 3;
trlid(trl4) = 4; %assign trial id

%adjust based on actual trial information
% Fs = 30000; %obtained 
% trial_length = 2;
% pre_trial = 0.2;
% post_trial = trial_length-pre_trial;
tb =-pre_trial*1000:1:(post_trial)*1000;


nchans = size(MUA,1);
ModS = zeros(nchans,4);
ModM = zeros(nchans,4);
Ons = zeros(nchans,4);
Offs = zeros(nchans,4);

clear SNR
if analyzedata
    if notch
        mon = 60; %Frequency of minitor
        no = 2;
        wn = [mon-2 mon+2]./(Fs./2);
        [b1,a1] = butter(no,wn,'stop');
        wn = [mon*2-2 mon*2+2]./(Fs./2);
        [b2,a2] = butter(no,wn,'stop');
        wn = [mon*3-2 mon*3+2]./(Fs./2);
        [b3,a3] = butter(no,wn,'stop');
    end

    if fig
        figure,hold on
        colind = jet(96);
    end

    for ch = 1:size(MUA,1)
        clear mua
        mua = [];
        buf = squeeze(MUA(ch,:,:));

        if notch
            %Filter out the monitor frequency
            buf = filtfilt(b1,a1,buf);
            buf = filtfilt(b2,a2,buf);
            buf = filtfilt(b3,a3,buf);
        end
        mua = [mua,buf];

        for n = 1:4
            %Get trials with this motion direction
            f = find(trlid == n);

            %Average them
            MUAm(n,:) = mean(mua(f,:),'omitnan');
            MUAs(n,:) = std(mua(f,:),'omitnan')./sqrt(length(f));

            %Get noise levels before smoothing
            BaseT = find(tb >-100 & tb < 0);
            Base = mean(MUAm(n,BaseT),'omitnan');
            BaseS = std(MUAm(n,BaseT),'omitnan');

            %Smooth it to get a maximum...
            gt = find(tb>0 & tb<2000);
            sm = smooth(MUAm(n,gt),10);
            [mx,mi] = max(sm);
            Scale = mx-Base;

            %Is the max significantly different to the base?
            SigDif(ch,n) = mx > (Base+(1.*BaseS));


             %Automatic RF fitting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             %Fit a guassian to it...
             %Get fitting ranges
             RespT = find(tb > 0);
             %time
             td = tb(RespT);
             index = 1:5:length(td);
             pxm = td(index);
             sdind = [10:10:200];

             z = 0;
             for i = 1:length(pxm)
                 for j = 1:length(sdind)
                     norm = normpdf(td,pxm(i),sdind(j));
                     norm = norm./(max(norm));
                     model = (norm.*Scale)+Base;
                     diff = sum((model-MUAm(n,RespT)).^2);
                     z = z+1;
                     ind(z,:) = [diff,pxm(i),sdind(j)];
                 end
             end
             [x,y] = min(ind(:,1));
             ModM(ch,n) = ind(y,2);
             ModS(ch,n) = ind(y,3);
             Ons(ch,n) = ModM(ch,n)-(1.95.*ModS(ch,n));
             Offs(ch,n) = ModM(ch,n)+(1.95.*ModS(ch,n));
             %Bet model
             norm = normpdf(td, ModM(ch,n), ModS(ch,n));
             norm = norm./(max(norm));
             bestmodel = (norm.*Scale)+Base;


            if fig
                subplot(2,2,n),plot(td,MUAm(n,RespT),'b',td,bestmodel,'r')
            end


            if fig
                h = line([Ons(ch,n),Ons(ch,n)],get(gca,'YLim'));
                set(h,'Color',[1 0 1])
                h = line([Offs(ch,n),Offs(ch,n)],get(gca,'YLim'));
                set(h,'Color',[1 0 1])
            end

            SNR(ch,n)=Scale/BaseS;
        end
        ch % pause here for individual curve fit

    end

    if saveout
    save(savename,'Ons','Offs','SigDif','speed','SNR','Log','pixperdeg')
    end
else
    load(savename)
end
    
    
%SKIP TO HERE
SNR = mean(SNR,2);
SNRcutoff = 4;
figure;
set(gcf,'Position',[0 0 1024 768])
for ch = 1:size(MUA,1)
    %Only plot channels where all directon were signifcant and SNR is high enough
    if sum(SigDif(ch,:)) == 4 && SNR(ch)>SNRcutoff %SNR(ch)>SNRcutoff
        %Now distance = speed*time
        %This gives distanbce travelled by bar in pixels before the onset and
        %offset
        onsdist = speed.*(Ons(ch,:)./1000);
        offsdist = speed.*(Offs(ch,:)./1000);
        
        %Stimuli 1-4 go
        %1 = horizontal left-to-right (180 deg),
        %2 = bottom to top 270
        %3 = right to left 0
        %4 = top to bottom 90
        angles = [180 270 0 90];
        
        %Get starting position of bars
        sx = xo+(bardist./2).*cosd(angles);
        sy = yo+(bardist./2).*sind(angles);
        
        %Angular distance moved
        %(direction is opposite to angle of starting
        %position)
        on_angx = onsdist.*cosd(180-angles);
        on_angy = onsdist.*sind(angles);
        off_angx = offsdist.*cosd(180-angles);
        off_angy = offsdist.*sind(angles);
        
        %So the on and off points are starting position + angular distance...
        onx = sx+on_angx;
        ony = sy-on_angy;
        offx = sx+off_angx;
        offy = sy-off_angy;
        
        %get RF vboundaries
        bottom = (ony(2)+offy(4))./2;
        right = (onx(1)+offx(3))./2;
        top =   (ony(4)+offy(2))./2;
        left =   (onx(3)+offx(1))./2;
        
        RF.centerx(ch) = (right+left)./2;
        RF.centery(ch) = (top+bottom)./2;
        
        RF.sz(ch) = sqrt(abs(top-bottom).*abs(right-left));
        RF.szdeg(ch) = sqrt(abs(top-bottom).*abs(right-left))./pixperdeg;
        
        XVEC1 = [left  right  right  left  left];
        YVEC1 = [bottom bottom  top top  bottom];
        
        RF.XVEC1(ch,:) = XVEC1;
        RF.YVEC1(ch,:) = YVEC1;
        
        h = line(XVEC1,YVEC1);
        set(h,'Color',bankcol(bank(ch),:))
%         axis square
        hold on
        scatter(0,0,'r','f')
        scatter(RF.centerx(ch),RF.centery(ch),36,bankcol(bank(ch),:),'f')
        axis([-512 512 -384 384])
        xticks([-20*pixperdeg -15*pixperdeg -10*pixperdeg -5*pixperdeg 0 ...
            5*pixperdeg 10*pixperdeg 15*pixperdeg 20*pixperdeg])
        xticklabels({'-20','-15','-10','-5','0','5','10','15','20'})
        yticks([-15*pixperdeg -10*pixperdeg -5*pixperdeg 0 ...
            5*pixperdeg 10*pixperdeg 15*pixperdeg])
        yticklabels({'-15','-10','-5','0','5','10','15'})
        hold on,scatter(sx,sy)
%         disp(['channel: ' ,num2str(ch)])
%         disp(['centerx = ',num2str(RF.centerx(ch))])
%         disp(['centery = ',num2str(RF.centery(ch))])
        %position in degrees
        RF.ang(ch)= atand(RF.centery(ch)./RF.centerx(ch));
        
        %pix2deg conversion
        RF.ecc(ch) = sqrt(RF.centerx(ch).^2+RF.centery(ch).^2)./pixperdeg;
        
        % disp(['Angle = ',num2str(RF_ang(ch))])
        disp(['Ecc = ',num2str(RF.ecc(ch))])
        disp(' ')
        
        text(RF.centerx(ch),RF.centery(ch),num2str(ch))
        %Save out centx
    else
        %Bad channels get set to NaN
        RF.centerx(ch)=NaN;
        RF.centery(ch)=NaN;
        RF.sz(ch)=NaN;
        RF.szdeg(ch)=NaN;
        RF.XVEC1(ch,:)=NaN(1,5);
        RF.YVEC1(ch,:)=NaN(1,5);
        RF.ang(ch) = NaN;
        RF.ecc(ch) = NaN;
    end
end

if ~isempty(RFx)
    %SCatter on markers
    hold on,scatter(RFx,RFy,'MarkerFaceColor',[0.8 0.8 0.8])
    for i = 1:length(RFx)
        text(RFx(i),RFy(i),(['x=' num2str(RFx(i)) ', y=' num2str(RFy(i))]))
    end
end

figure;
for a = 1:8
    subplot(2,4,a),hold on
    for ch = find(bank==a)
        h = line(RF.XVEC1(ch,:)./pixperdeg,RF.YVEC1(ch,:)./pixperdeg);
        set(h,'Color',bankcol(a,:))
        scatter(RF.centerx(ch)./pixperdeg,RF.centery(ch)./pixperdeg,36,bankcol(a,:),'f')
        text(RF.centerx(ch)./pixperdeg,RF.centery(ch)./pixperdeg,num2str(ch))
    end
    scatter(0,0,'r','f')
    axis square
    axis([-20 20 -20 20])
end

disp(['%%% ... done in: ',num2str(round(toc)),'s']) 
