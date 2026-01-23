%% add tools
scriptfld = pwd();
openephys_fld = fullfile(scriptfld,'OpenEphys','open-ephys-matlab-tools');
addpath(genpath(openephys_fld));

%% set paths
%data_fld = '/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_collection';
%data_fld = '/media/chris/CK4TB/Neuropixels_NHP/Data_collection';
if ismac
    data_fld = '/Users/chris/Dropbox/CURRENT_PROJECTS/NEUROPIXELS/Sample_Data';
else
    data_fld = '/home/chris/Documents/CURRENT_PROJECTS/NEUROPIXELS/Sample_Data';
end

subject = 'm032';
day = '20240826';
expt = 1; reci = 1;
ProbeLabel = 'Record Node 101';

%% re-order data for the tools to work
data_root_org = fullfile(data_fld,subject,day);
data_root_probe = fullfile(data_fld,subject,day,ProbeLabel);
% put experiments under probe lable
if ~isdir(data_root_probe)
    mkdir(data_root_probe);
    movefile(fullfile(data_root_org,'experiment*'),data_root_probe);
end

% put it back
%movefile(fullfile(data_root_probe,'experiment*'),data_root_org);
%rmdir(data_root_probe)

% add the linked folder to the path so that the toolbox can search it
addpath(genpath(fullfile(data_root_probe,ProbeLabel)));
session = Session(fullfile(data_fld,subject,day));
% session now contain the recording info

%% node refers to recording node 
node = session.recordNodes{1};

% each entry of node.recordings{} is a a run 
% note that runs get indexed in ascending order; use the filepath under
% node.recordings{i}.directory to identify the real run number
rec = node.recordings{reci}; % refers to a run/recording
disp(rec); % see info

% identify what the datastreams are
disp 'Continuous'
for i = 1:length(rec.continuous)
    disp(['Stream ',num2str(i),': ' rec.info.continuous(i).stream_name]);
end

% Pick the BNC stream (PXIe-6341)
% in this case it's stream 3
cs = rec.continuous;
cskeys = keys(cs);
for i = 1:length(cskeys)
    disp(['Continous data ',num2str(i),': ' cskeys{i}]);
end
BNCstream = rec.continuous(cskeys{1});
APstream = rec.continuous(cskeys{2});

%% convert timebase
achan = 5;
ts = 0:1/BNCstream.metadata.sampleRate:...
    length(BNCstream.timestamps)*(1/BNCstream.metadata.sampleRate);
ts(end)=[];
ts2 = BNCstream.timestamps - BNCstream.timestamps(1);
%trange=(1:1e5);plot(ts2(trange),BNCstream.samples(achan,trange))

d=BNCstream.samples(achan,:);
dd = abs(diff(d));
ts3 = ts2(1:end-1);
%plot(ts3,dd)

%% find the first > 600 for each trial
dd6k = dd>6000;
dd6kidx = find(dd6k>0);
tsd = diff(ts3(dd6k));
%histogram(tsd,100)

dd6kidx = dd6kidx(2:end);
idx = [dd6kidx(1) dd6kidx(tsd>0.500)];

% plot(ts3,dd)
% hold on
% xline(ts3(idx),'r')
% yline(6000,'g')

%% check what info is available as events
% looked at XL's scripts
disp 'Events'
for i = 1:length(rec.ttlEvents)
    disp(['Stream ',num2str(i),': ' rec.info.events{i}.stream_name]);
end
D = rec.ttlEvents(cskeys{1});
bits = {D.line double(D.timestamp)-BNCstream.timestamps(1) ...
    D.sample_number D.state};

% correct trials == 6
% trial bit = 4

TrOn = bits{2}(bits{1}==4 & bits{4}==1); % trial
b8 = bits{2}(bits{1}==8  & bits{4}==1); % reward
b1 = bits{2}(bits{1}==1  & bits{4}==1);

%TrC = bits{2}(bits{1}==6);

figure; 
subplot(2,1,1); hold on;
plot(ts3,dd,'Color',[0.2 0.2 0.2])
xline(ts3(idx),'r','LineWidth',3)
yline(6000,'g','LineWidth',1)
xline(TrOn,'c','LineWidth',3)
xline(b8,'m','LineWidth',1)
%xline(b1,'y','LineWidth',1) %timer

set(gca,'xlim',[0 30]);

%% get the Tracker log
% load(fullfile(data_fld,subject,day,'run-006_T-143435',...
%     'sub-m032_ses-20251022_task-ct_stim-fgdots_run-006'));
load(fullfile(data_fld,subject,day,ProbeLabel,'experiment1',...
    'recording3','run-003_T-154348',...
    'sub-m032_ses-20240826_task-ct_stim-fgdots_run-003'));
%%
TrLog = []; toffset = 0;
MotionLog = []; DelayLog = [];
ITI0Log = []; ITI1Log = [];
TargLog = []; CorrLog = [];
RewStart = [];
for i=1:length(Log.events)
    if strcmp(Log.events(i).type,'trial_start')
        TrLog = [TrLog; Log.events(i).t + toffset Log.events(i).trialn];
        motion_logged = false;
    elseif strcmp(Log.events(i).type,'motion_start') 
        if ~motion_logged
            MotionLog = [MotionLog; Log.events(i).t + toffset];
        end
        motion_logged = true;
    elseif strcmp(Log.events(i).type,'delay_start') 
        DelayLog = [DelayLog; Log.events(i).t + toffset];
    elseif strcmp(Log.events(i).type,'reward_start') 
        RewStart = [RewStart; Log.events(i).t + toffset];    
    elseif strcmp(Log.events(i).type,'iti_start')
        ITI0Log = [ITI0Log; Log.events(i).t + toffset];
    elseif strcmp(Log.events(i).type,'iti_stop')
        ITI1Log = [ITI1Log; Log.events(i).t + toffset];
        toffset = toffset + Log.events(i).t;
    elseif strcmp(Log.events(i).type,'targ_start')
        TargLog = [TargLog; Log.events(i).t + toffset];
    elseif strcmp(Log.events(i).type,'correct')
        CorrLog = [CorrLog; Log.events(i).t + toffset];
    end
end

subplot(2,1,2); hold on;
% get timing offset OpenEphys / Tracker
dt = TrOn(1) - TrLog(1,1);

xline(TrLog(:,1)+dt,'c');
xline(MotionLog+dt,'g');
xline(DelayLog+dt,'c--');
%xline(TargLog+dt,'c');
xline(RewStart+dt,'m');
set(gca,'xlim',[0 30]);

% NB! onset times are not referenced to Experiment start
% check runstim and correct if desired

%% get the average response over all trials
basechan = 100; numchan = 10;
ch = basechan:basechan+numchan-1;
tw = [-0.200 0.800];

APt = single(APstream.timestamps);
APt = APt(1:30:end)-APt(1);

RAW = single(APstream.samples(ch,:));
RRAW = RAW-mean(RAW);
FAP = zeros(size(RRAW));

% get MUA envelope
Fs = APstream.metadata.sampleRate;
Fbp = [500, 5000];
N = 2;
Fn = Fs/2;
Fl = 200;
[B,A] = butter(N,[min(Fbp)/Fn max(Fbp)/Fn]);
[D,C] = butter(N,Fl/Fn,'low');

for i=1:size(RAW,1)
    disp(['Channel ' num2str(ch(i))])
    S=RAW(i,:);
    dum1 = filtfilt(B,A,S);
    dum2 = abs(dum1);
    muafilt = filtfilt(D,C,dum2);
    
    % remove 50 Hz (line) and 60 Hz (monitor)
    buttLoop = muafilt;
    for j=1:5
        d = designfilt('bandstopiir','FilterOrder',2,...
            'HalfPowerFrequency1',50*j-5,...
            'HalfPowerFrequency2',50*j+5,...
            'DesignMethod','butter','SampleRate',Fs);
        buttLoop = filtfilt(d,buttLoop);
    end
    for j=1:5
        d = designfilt('bandstopiir','FilterOrder',2,...
            'HalfPowerFrequency1',60*j-5,...
            'HalfPowerFrequency2',60*j+5,...
            'DesignMethod','butter','SampleRate',Fs);
        buttLoop = filtfilt(d,buttLoop);
    end
    FAP(i,:) = buttLoop';
end

APd = downsample(FAP',30)';
%clear RAW muafilt S dum1 dum2 FAP
ns = find(APt<=1,1,'last');
TRIALS=[];

% select correct trials
cTrOn=[]; targtype=1;
for r = 1:length(RewStart)
    tRew = RewStart(r);
    sT = find(TrLog(:,1)<tRew,1,'last');
    cTrOn = [cTrOn;TrLog(sT,1) sT];   
end
tn = targidx(targidx(:,2)==targtype,1);
scTrOn = cTrOn(ismember(cTrOn(:,2),tn),1);

for i=1:length(scTrOn)
    tr0 = scTrOn(i)+tw(1);
    si = find(APt>=tr0,1,'first');
    TRIALS = cat(3,TRIALS,APd(:,si:si+ns));
end

%size(TRIALS)
avgT = median(TRIALS,3);

%%
figure; hold on;
for c = 1:length(ch)
    t = APt(1:1+ns)+tw(1);
    y = avgT(c,:);
    tsel = t>-0.1 & t < 0;
    yy = (y-mean(y(tsel)))./std(y(tsel));
    plot(t,smooth(yy,33));
end
xline(0,'-');
xline(0.1,'--');
xline(0.2,'-');


%% put data back in the original file structure
rmpath(genpath(fullfile(data_root_probe,ProbeLabel)));
movefile(fullfile(data_root_probe,'experiment*'),data_root_org);
rmdir(data_root_probe)