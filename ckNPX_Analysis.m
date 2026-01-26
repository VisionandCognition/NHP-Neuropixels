%% add tools
scriptfld = pwd();
openephys_fld = fullfile(scriptfld,'OpenEphys','analysis');
addpath(genpath(openephys_fld));

%% set paths
%data_fld = '/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_collection';
%data_fld = '/media/chris/CK4TB/Neuropixels_NHP/Data_collection';
if ismac
    data_fld = '/Users/chris/Dropbox/CURRENT_PROJECTS/NEUROPIXELS/Sample_Data';
else
    data_fld = '/home/chris/Documents/CURRENT_PROJECTS/NEUROPIXELS/Sample_Data';
end

log_fld = data_fld;
% which rec session
sn = 1;

%% add recordings info
ckRecList;

%% read json 
jsonFile = fullfile(data_fld,Session(sn).Monkey,Session(sn).Day,...
    ['recording' num2str(Session(sn).RecN)],'structure.oebin');
info = jsondecode(fileread(jsonFile));
info = info.continuous;
info_cont = info(1); info_bnc = info(3);
numChannels = info_cont.num_channels; 

%% get ttl
ttlDir = fullfile(data_fld,Session(sn).Monkey,Session(sn).Day,...
    ['recording',num2str(Session(sn).RecN)],'events',...
    'NI-DAQmx-114.PXIe-6341','TTL');
lines = readNPY(fullfile(ttlDir, 'states.npy'));
sampleNumbers = readNPY(fullfile(ttlDir, 'sample_numbers.npy'));
timestamps = readNPY(fullfile(ttlDir, 'timestamps.npy'));
numEvents = length(lines);
ttls = DataFrame(lines, sampleNumbers, timestamps, lines > 0, ...
    'VariableNames', {'line','sample_number','timestamp', 'state'});

%% get Tracker log
logfile = dir(...
    fullfile(data_fld,Session(sn).Monkey,Session(sn).Day,...
    ['run-',sprintf('%03d', Session(sn).RunN)  ,'*'],'sub-*.mat'));
load(fullfile(logfile.folder,logfile.name),'Log','Stm','Par');

%% get trial information
bdStart = find(ttls.line==4);
rewStart = find(ttls.line==8);
bdAll = bdStart;
logtrialn = Log.trial(length(Log.trial)).Trlnum;
assert(length(bdAll)==logtrialn,...
    'There are (%d) recorded trials while (%d) trials were logged',...
    length(bdAll),length(Log.trial));
BitsTrialOnTS = ttls.timestamp(bdStart); %trial onset timestamps
nBitsTrialOn = numel(BitsTrialOnTS);
BitsRewOnTS = ttls.timestamp(rewStart); %trial onset timestamps

%% get AP
apDir = fullfile(data_fld,Session(sn).Monkey,Session(sn).Day,...
    ['recording',num2str(Session(sn).RecN)],'continuous','Neuropix-PXI-100.ProbeA-AP');
buf = memmapfile(fullfile(apDir, 'continuous.dat'), 'Format', 'int16');
data.samples = reshape(buf.Data, [numChannels, length(buf.Data)/numChannels]);
data.sampleNumbers = readNPY(fullfile(apDir, 'sample_numbers.npy'));
data.timestamps = readNPY(fullfile(apDir, 'timestamps.npy'));
Fs = info_cont.sample_rate; FsAP=Fs;

%% get photodiode
bncDir = fullfile(data_fld,Session(sn).Monkey,Session(sn).Day,...
    ['recording',num2str(Session(sn).RecN)],'continuous','NI-DAQmx-114.PXIe-6341');
buf = memmapfile(fullfile(bncDir, 'continuous.dat'), 'Format', 'int16');
nchan = info_bnc.num_channels; 
bnc.samples = reshape(buf.Data, [nchan, length(buf.Data)/nchan]);
bnc.sampleNumbers = readNPY(fullfile(bncDir, 'sample_numbers.npy'));
bnc.timestamps = readNPY(fullfile(bncDir, 'timestamps.npy'));
Fs = info_bnc.sample_rate; FsBNC=Fs;

%% Get trial moments
apTimestamps = data.timestamps;
%ap index closest to stim onset
Trls = nan(length(BitsTrialOnTS),1);
for i = 1:length(BitsTrialOnTS)
    [val,idx] = min(abs(apTimestamps-BitsTrialOnTS(i)));
    Trls(i) = idx;
end
Trials = Trls;

%% Get trial starts from photodiode
pd = bnc.samples(5,:);
pd_ts = bnc.timestamps;

dd = abs(pd);
dd6k = dd>6000;
dd6kidx = find(dd6k>0);
tsd = diff(pd_ts(dd6k));
%figure;histogram(tsd,100)

dd6kidx = dd6kidx(2:end);
pdidx = [dd6kidx(1) dd6kidx(tsd>0.500)];
pd_trl = pd_ts(pdidx);

%%  get the ap sample idx for pd trial starts
% get the ap sample idx for correct trials
TrlsPD = nan(length(pd_trl),1);
for i = 1:length(pd_trl)
    [val,idx] = min(abs(apTimestamps-pd_trl(i)));
    TrlsPD(i) = idx;
end 

%% select correct trials pd start
cTrl = nan(length(BitsRewOnTS),1); 
for i = 1:length(BitsRewOnTS)
    sidx = find(pd_trl<BitsRewOnTS(i),1,'last');
    cTrl(i) = pd_trl(sidx);
end

% get the ap sample idx for correct trials
cTrlsPD = nan(length(cTrl),1);
for i = 1:length(cTrl)
    [val,idx] = min(abs(apTimestamps-cTrl(i)));
    cTrlsPD(i) = idx;
end

TrCorrBool = ismember(TrlsPD,cTrlsPD);

%% Plot for insight
figure; hold on
plot(pd_ts-(pd_ts(1)),dd,'Color',[0.6 0.6 0.6])
hold on
xline(pd_ts(pdidx)-pd_ts(1),'r')
xline(cTrl-pd_ts(1),'g')
xline(BitsRewOnTS-pd_ts(1),'m--')
yline(6000,'b')
set(gca,'xlim',[0 20])

%% get channel information 
settingsFile = fullfile(data_fld,Session(sn).Monkey,Session(sn).Day,'settings.xml');
[s] = xml2struct(settingsFile);
chypos = s.SETTINGS.SIGNALCHAIN{1,1}.PROCESSOR{1,1}.EDITOR.CUSTOM_PARAMETERS.NP_PROBE.ELECTRODE_YPOS.Attributes; %!processor order in signal chain might change
ch = fieldnames(chypos);
ypos = zeros(length(ch),2);
for i=1:length(ch)
    ypos(i,2)=str2num(getfield(chypos,ch{i})); %yposition
    ypos(i,1)=str2num(extractAfter(ch{i},"CH")); %channel number
end
ypos_sort=sortrows(ypos,2); %channel ordered by yposition
CHorder=ypos_sort(:,1); %channel number ordered by yposition

%% Get envelope MUA
Fs = FsAP;
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

%% Extract and filter trials
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
%MUA = zeros(downsamps,384,nBitsTrialOn);
MUA = zeros(downsamps,384,length(cTrlsPD));

first_sample = data.sampleNumbers(1);
L = samps_per_trial/Fs; 
smps = 0:1:(samps_per_trial-1);
f = smps/L;

%vectorize
fmat = repmat(f',1,384);
phmat = repmat(ph,samps_per_trial,1);

fprintf('Processing Trial:  ')
for k = 1:nBitsTrialOn
%for k = 1:length(cTrlsPD)
    
    fprintf([' ' num2str(k)])
    samp_st = Trials(k)-pre_trialstart+1;
    samp_ed = Trials(k)+post_trialstart+1;
    %samp_st = cTrlsPD(k)-pre_trialstart+1;
    %samp_ed = cTrlsPD(k)+post_trialstart+1;
    cutdata = double(data.samples(1:384,samp_st:samp_ed))'*0.195;

    % HP filter
    cutdata = filtfilt(b_hp,a_hp,cutdata);

    % Phase align
    fbuf=fft(cutdata);
    fbuf = exp(-1j.*2.*pi.*fmat.*phmat).*fbuf;
    aligndata=ifft(fbuf,'symmetric');

    % Low-pass at 5000hz
    aligndata = filtfilt(b_lp,a_lp,aligndata);

    % Now destripe this chunck
    destriped = filtfilt(b_ds,a_ds,aligndata')';

    % Convert to MUA
    buf = filtfilt(b_mualp,a_mualp,abs(destriped)); %Now Take abs value and low-pass
    for ch = 1:384
        MUA(:,ch,k) = decimate(buf(:,ch),30);
    end
end
fprintf('\nMUA processed\n');

%% Take mean response per trial and subtract channel baseline
base = mean(squeeze(mean(MUA(tbds>-0.15&tbds<0,:,:))),2);
MUA = MUA-repmat(base',size(MUA,1),1,size(MUA,3));
%figure,plot(tbds,squeeze(mean(mean(MUA(:,:,:),3),2)))

MUA = permute(MUA,[2,3,1]); % reshape mua into channel x trial x trial duration
MUA = MUA(CHorder+1,:,:);

figure; 
subplot(1,2,1);
imagesc(squeeze(mean(MUA,2))),caxis([-1 3]),xlabel('Time (ms)'),ylabel('Channels')
set(gca,'XTick',[0:100:1100],'XTickLabel',[-pre_trial*1000:100:900]);
hold on;
xline(pre_trial*1000,'w--','LineWidth', 2);
xline((pre_trial*1000)+100,'w-','LineWidth', 2);
xline((pre_trial*1000)+200,'w-','LineWidth', 2);

subplot(1,2,2)
imagesc(squeeze(mean(MUA(:,TrCorrBool),2))),caxis([-1 3]),xlabel('Time (ms)'),ylabel('Channels')
set(gca,'XTick',[0:100:1100],'XTickLabel',[-pre_trial*1000:100:900]);
hold on;
xline(pre_trial*1000,'w--','LineWidth', 2);
xline((pre_trial*1000)+100,'w-','LineWidth', 2);
xline((pre_trial*1000)+200,'w-','LineWidth', 2);

%% save out
muasave = fullfile(data_fld,Session(sn).Monkey,Session(sn).Day,...
    ['recording',num2str(Session(sn).RecN)],['mua',num2str(Session(sn).RecN)]);
save(muasave,'MUA','TrCorrBool','-v7.3');














