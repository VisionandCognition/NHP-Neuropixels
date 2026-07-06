%% add tools
fld.script = pwd();
fld.openephys = fullfile(fld.script,'OpenEphys','analysis');
addpath(genpath(fld.openephys));

%% set paths
%fld.data = '/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_collection';
%fld.data = '/media/chris/CK4TB/Neuropixels_NHP/Data_collection';
if ismac
    fld.data = '/Users/chris/Dropbox/CURRENT_PROJECTS/NEUROPIXELS/Sample_Data';
else
    fld.data = '/home/chris/Documents/CURRENT_PROJECTS/NEUROPIXELS/Sample_Data';
end

fld.log = fld.data;
% which rec session
sess.n = 1;

%% add recordings info
ckRecList;

%% read json 
sess.jsonFile = fullfile(fld.data,Session(sess.n).Monkey,Session(sess.n).Day,...
    ['recording' num2str(Session(sess.n).RecN)],'structure.oebin');
info = jsondecode(fileread(sess.jsonFile));
info = info.continuous;
info_ap = info(1); info_bnc = info(3);
ap.numChannels = info_ap.num_channels; 

%% get ttl
fld.ttlDir = fullfile(fld.data,Session(sess.n).Monkey,Session(sess.n).Day,...
    ['recording',num2str(Session(sess.n).RecN)],'events',...
    'NI-DAQmx-114.PXIe-6341','TTL');
ttl.lines = readNPY(fullfile(fld.ttlDir, 'states.npy'));
ttl.sampleNumbers = readNPY(fullfile(fld.ttlDir, 'sample_numbers.npy'));
ttl.timestamps = readNPY(fullfile(fld.ttlDir, 'timestamps.npy'));
ttl.numEvents = length(ttl.lines);
ttl.df = DataFrame(ttl.lines, ttl.sampleNumbers, ttl.timestamps, ttl.lines > 0, ...
    'VariableNames', {'line','sample_number','timestamp', 'state'});

%% get Tracker log
sess.logfile = dir(...
    fullfile(fld.data,Session(sess.n).Monkey,Session(sess.n).Day,...
    ['run-',sprintf('%03d', Session(sess.n).RunN)  ,'*'],'sub-*.mat'));
load(fullfile(sess.logfile.folder,sess.logfile.name),'Log','Stm','Par');

%% get trial information
ttl.bdStart = find(ttl.df.line==4);
ttl.rewStart = find(ttl.df.line==8);
ttl.bdAll = ttl.bdStart;
Log.trialn = Log.trial(length(Log.trial)).Trlnum;
assert(length(ttl.bdAll)==Log.trialn,...
    'There are (%d) recorded trials while (%d) trials were logged',...
    length(ttl.bdAll),length(Log.trial));
ttl.BitsTrialOnTS = ttl.df.timestamp(ttl.bdStart); %trial onset timestamps
ttl.nBitsTrialOn = numel(ttl.BitsTrialOnTS);
ttl.BitsRewOnTS = ttl.df.timestamp(ttl.rewStart); %trial onset timestamps

%% get AP
sess.apDir = fullfile(fld.data,Session(sess.n).Monkey,Session(sess.n).Day,...
    ['recording',num2str(Session(sess.n).RecN)],'continuous','Neuropix-PXI-100.ProbeA-AP');
buf = memmapfile(fullfile(sess.apDir, 'continuous.dat'), 'Format', 'int16');
ap.samples = reshape(buf.Data, [ap.numChannels, length(buf.Data)/ap.numChannels]);
ap.sampleNumbers = readNPY(fullfile(sess.apDir, 'sample_numbers.npy'));
ap.timestamps = readNPY(fullfile(sess.apDir, 'timestamps.npy'));
ap.Fs = info_ap.sample_rate; 

%% get photodiode
sess.bncDir = fullfile(fld.data,Session(sess.n).Monkey,Session(sess.n).Day,...
    ['recording',num2str(Session(sess.n).RecN)],'continuous','NI-DAQmx-114.PXIe-6341');
buf = memmapfile(fullfile(sess.bncDir, 'continuous.dat'), 'Format', 'int16');
bnc.nchan = info_bnc.num_channels; 
bnc.samples = reshape(buf.Data, [bnc.nchan, length(buf.Data)/bnc.nchan]);
bnc.sampleNumbers = readNPY(fullfile(sess.bncDir, 'sample_numbers.npy'));
bnc.timestamps = readNPY(fullfile(sess.bncDir, 'timestamps.npy'));
bnc.Fs = info_bnc.sample_rate;

%% Get trial moments
%ap index closest to stim onset
ttl.Trls = nan(length(ttl.BitsTrialOnTS),1);
for i = 1:length(ttl.BitsTrialOnTS)
    [val,idx] = min(abs(ap.timestamps-ttl.BitsTrialOnTS(i)));
    ttl.Trls(i) = idx;
end
%Trials = Trls;

%% Get trial starts from photodiode
pd.samples = bnc.samples(5,:);
pd.timestamps = bnc.timestamps;

dd = abs(pd.samples);
dd6k = dd>6000;
dd6kidx = find(dd6k>0);
tsd = diff(pd.timestamps(dd6k));
%figure;histogram(tsd,100)

dd6kidx = dd6kidx(2:end);
pd.pdidx = [dd6kidx(1) dd6kidx(tsd>0.500)];
pd.trl.timestamps = pd.timestamps(pdidx);

%%  get the ap sample idx for pd trial starts
pd.trl.apidx = nan(length(pd.trl.timestamps),1);
for i = 1:length(pd.trl.timestamps)
    [val,idx] = min(abs(ap.timestamps-pd.trl.timestamps(i)));
    pd.trl.apidx(i) = idx;
end 

%% select correct trials pd start
pd.ctrl.pdts = nan(length(ttl.BitsRewOnTS),1); 
for i = 1:length(ttl.BitsRewOnTS)
    sidx = find(pd.trl.timestamps<ttl.BitsRewOnTS(i),1,'last');
    pd.ctrl.pdts(i) = pd.trl.timestamps(sidx);
end

% get the ap sample idx for correct trials
pd.ctrl.apidx = nan(length(pd.ctrl.pdts),1);
for i = 1:length(pd.ctrl.pdts)
    [val,idx] = min(abs(ap.timestamps-pd.ctrl.pdts(i)));
    pd.ctrl.apidx(i) = idx;
end

pd.trl.TrCorrBool = ismember(pd.trl.apidx,pd.ctrl.apidx);

%% Plot for insight
figure; hold on
plot(pd.timestamps-(pd.timestamps(1)),abs(pd.samples),'Color',[0.6 0.6 0.6])
hold on
xline(pd.timestamps(pd.pdidx)-pd.timestamps(1),'r')
xline(pd.ctrl.pdts-pd.timestamps(1),'g')
xline(ttl.BitsRewOnTS-pd.timestamps(1),'m--')
yline(6000,'b')
set(gca,'xlim',[0 20])

%% get channel information 
sess.setFile = fullfile(fld.data,Session(sess.n).Monkey,...
    Session(sess.n).Day,'settings.xml');
[sess.s] = xml2struct(sess.setFile);
sess.chypos = sess.s.SETTINGS.SIGNALCHAIN{1,1}.PROCESSOR{1,1}.EDITOR.CUSTOM_PARAMETERS.NP_PROBE.ELECTRODE_YPOS.Attributes; %!processor order in signal chain might change
sess.ch = fieldnames(sess.chypos);
sess.ypos = zeros(length(sess.ch),2);
for i=1:length(sess.ch)
    sess.ypos(i,2) = str2num(getfield(sess.chypos,sess.ch{i})); %yposition
    sess.ypos(i,1) = str2num(extractAfter(sess.ch{i},"CH")); %channel number
end
sess.ypos_sort = sortrows(sess.ypos,2); %channel ordered by yposition
sess.CHorder = sess.ypos_sort(:,1); %channel number ordered by yposition

%% Get envelope MUA
[f.b_hp,f.a_hp] = butter(3,300/(ap.Fs/2),'high'); 
[f.b_lp,f.a_lp] = butter(3,5000/(ap.Fs/2),'low'); %remove high frequency noise
[f.b_mualp,f.a_mualp] = butter(3,200/(ap.Fs/2),'low');
[f.b_ds,f.a_ds] = butter(3,0.01/1,'high'); %Spatial filter (currently tuned by eye)

% NB! reconsider these filters [CK]

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
ap.ph = (rem(cycle-1,12)./13).*(1/ap.Fs);

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














