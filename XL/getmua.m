uselocal = 1; % if not, use data on server, log files always on server
getmua = 1;
tic;

Monkey = 'm032';
Day = '20251022'; % YYYYMMDD
ExpN = 1;
RecN = 12;
RunN = 11; % two digits
Dataloc = '\\vs03.herseninstituut.knaw.nl\VS03-VandC-6\Neuropixels_NHP\Data_collection\';
localdayDir = 'E:\OpenEphys\m032-2025-10-22_12-43-22\Record Node 110'; %change day and timestamp

figuredir = [Dataloc,Monkey,'\',Day,'\recording',num2str(RecN),'\']; %Location for the figures

Logdir = [Dataloc,Monkey,'\',Day,'\run-0',num2str(RunN),'*\'];
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
bdStart = find(ttls.line==4);

bdAll = bdStart;

logtrialn = Log.trial(length(Log.trial)).Trlnum;

assert(length(bdAll)==logtrialn,'There are (%d) recorded trials while (%d) trials were logged',length(bdAll),length(Log.trial));

BitsTrialOnTS = ttls.timestamp(bdStart); %trial onset timestamps
nBitsTrialOn = numel(BitsTrialOnTS);

%% get AP
buf = memmapfile(fullfile(apDir, 'continuous.dat'), 'Format', 'int16');
data.samples = reshape(buf.Data, [numChannels, length(buf.Data)/numChannels]);
data.sampleNumbers = readNPY(fullfile(apDir, 'sample_numbers.npy'));
data.timestamps = readNPY(fullfile(apDir, 'timestamps.npy'));
Fs = info.sample_rate;

%%
apTimestamps = data.timestamps;
%ap index closest to stim onset
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

%ADC
z = 0;
a = 1:24:384;
b = 2:24:384;
cycle = NaN(1,384);
for s = 1:12
    cycle(a+(s-1)*2) = s;
    cycle(b+(s-1)*2) = s;
end

%Convert cycle to phase shift per channel
ph = (rem(cycle-1,12)./13).*(1/Fs);

%%

trial_length = 1.1; %tdct, stim 100ms
pre_trial = 0.2;
post_trial = trial_length - pre_trial;
pre_trialstart = round(pre_trial.*Fs);
post_trialstart = round(post_trial.*Fs);

tmbs = -pre_trialstart:1:post_trialstart;
tb = tmbs./Fs;
samps_per_trial = numel(tmbs);
downs = 30;
downsamps= length(1:downs:samps_per_trial);
tbds = tb(1:downs:end);
MUA = zeros(downsamps,384,nBitsTrialOn);
first_sample = data.sampleNumbers(1);
L = samps_per_trial/Fs; 
smps = 0:1:(samps_per_trial-1);
f = smps/L;

%vectorize
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

%     for ch = 1:384
%         %bad channel
%         [f,p] = ez_powermeasure(buf,Fs);
%         high_power(k,ch) = mean(p(f>0.8*(Fs/2)));
%     end

    %% Low-pass at 5000hz
    aligndata = filtfilt(b_lp,a_lp,aligndata);

    %% Now desripe this chunck
    destriped = filtfilt(b_ds,a_ds,aligndata')';

    %% Convert to MUA
    buf = filtfilt(b_mualp,a_mualp,abs(destriped)); %Now Take abs value and low-pass
    for ch = 1:384
        MUA(:,ch,k) = decimate(buf(:,ch),30);
    end
    %t = toc;
    %disp(['Trial ',num2str(k),' processed in ',num2str(t),' seconds'])
end

disp(['MUA processed in ',num2str(toc),' seconds'])

%% Take mean response per trial and subtract channel baseline
base = mean(squeeze(mean(MUA(tbds>-0.15&tbds<0,:,:))),2);
MUA = MUA-repmat(base',size(MUA,1),1,size(MUA,3));

figure,plot(tbds,squeeze(mean(mean(MUA(:,:,:),3),2)))

MUA = permute(MUA,[2,3,1]); % reshape mua into channel x trial x trial duration
MUA = MUA(CHorder+1,:,:);

figure,imagesc(squeeze(mean(MUA,2))),caxis([-1 3]),xlabel('Time'),ylabel('Channels')

if saveout
save(muasave,'MUA','-v7.3')
end
