%% add tools
scriptfld = '/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels';
%scriptfld = '/Users/chris/Documents/EPHYS_ANALYSIS/NHP-Neuropixels';
openephys_fld = fullfile(scriptfld,'OpenEphys','open-ephys-matlab-tools');
addpath(genpath(openephys_fld));

%% set paths
%data_fld = '/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_collection';
%data_fld = '/media/chris/CK4TB/Neuropixels_NHP/Data_collection';
data_fld = '/home/chris/Documents/CURRENT_PROJECTS/NEUROPIXELS/Sample_Data';
%data_fld = '/Users/chris/Dropbox/CURRENT_PROJECTS/NEUROPIXELS/Sample_Data';

subject = 'm032';
day = '20251022';
expt = 1; rec = 6;
ProbeLabel = 'Record Node 101';

%% re-order data for the tools to work
data_root_org = fullfile(data_fld,subject,day);
data_root_probe = fullfile(data_fld,subject,day,ProbeLabel);
% put experiments under probe lable
mkdir(data_root_probe);
movefile(fullfile(data_root_org,'experiment*'),data_root_probe);
% put it back
%movefile(fullfile(data_root_probe,'experiment*'),data_root_org);
%rmdir(data_root_probe)
% add the linked folder to the path so that the toolbox can search it
addpath(genpath(fullfile(data_root_probe,ProbeLabel)));

% cannot create a symlink on the server?
session = Session(fullfile(data_fld,subject,day));
% session now contain the recording info

%% node refers to recording node 
node = session.recordNodes{1};

%% each entry of node.recordings{} is a a run 
% note that runs get indexed in ascending order; use the filepath under
% node.recordings{i}.directory to identify the real run number
rec = node.recordings{2}; % refers to a run/recording
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
trange=(1:1e5);
achan = 5;
ts = 0:1/BNCstream.metadata.sampleRate:...
    (length(BNCstream.timestamps)/BNCstream.metadata.sampleRate)-1;
ts2 = BNCstream.timestamps - BNCstream.timestamps(1);
%plot(ts2(trange),BNCstream.samples(achan,trange))

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

%% check what info is available on external input
figure;
for i=1:8
    subplot(8,1,i);
    plot(BNCstream.samples(i,:))
end

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
b8 = bits{2}(bits{1}==8  & bits{4}==1); % correct
b1 = bits{2}(bits{1}==1  & bits{4}==1);

%TrC = bits{2}(bits{1}==6);

figure; hold on;
plot(ts3,dd,'Color',[0.2 0.2 0.2])
xline(ts3(idx),'r','LineWidth',3)
yline(6000,'g','LineWidth',1)
xline(TrOn,'c','LineWidth',3)
xline(b8,'m','LineWidth',1)
%xline(b1,'y','LineWidth',1) %timer

set(gca,'xlim',[0 30]);

%% get the Tracker log
load(fullfile(data_fld,subject,day,'run-006_T-143435',...
    'sub-m032_ses-20251022_task-ct_stim-fgdots_run-006'));
%%
TrLog = [];
for i=1:length(Log.events)
    if strcmp(Log.events(i).type,'trial_start')
        TrLog = [TrLog; Log.events(i).t];
    end
end

TrLog2 = TrLog-TrLog(1)+TrOn(1);
%plot(TrLog2,10000*ones(size(TrLog2)),'ow','MarkerSize',10)

% NB! onset times are not referenced to Experiment start
% check runstim and correct if desired

%% 
corrdur =[];
for i = 1:length(b8)
    tidx = find(TrOn<b8(i),1,'last');
    corrdur = [corrdur; b8(i)-TrOn(tidx)];
end


%% get the average response over all trials
ch = 1:10;
tw = [-0.200 0.800];

APt = single(APstream.timestamps);
APt = APt(1:4:end)-APt(1);

RAW = single(APstream.samples(ch,:));
FAP = zeros(size(RAW));

% get MUA envelope
Fs = APstream.metadata.sampleRate;
Fbp = [500, 5000];
N = 2;
Fn = Fs/2;
Fl = 200;
[B,A] = butter(N,[min(Fbp)/Fn max(Fbp)/Fn]);
[D,C] = butter(N,Fl/Fn,'low');

for i=1:size(RAW,1)
    disp(['Channel ' num2str(i)])
    S=RAW(i,:);
    dum1 = filtfilt(B,A,S);
    dum2 = abs(dum1);
    muafilt = filtfilt(D,C,dum2);
    

    % remove 50 Hz (line) and 60 Hz (monitor)
    buttLoop = muafilt;
    for i=1:5
        d = designfilt('bandstopiir','FilterOrder',2,...
            'HalfPowerFrequency1',50*i-5,...
            'HalfPowerFrequency2',50*i+5,...
            'DesignMethod','butter','SampleRate',Fs);
        buttLoop = filtfilt(d,buttLoop);
    end
    for i=1:5
        d = designfilt('bandstopiir','FilterOrder',2,...
            'HalfPowerFrequency1',60*i-5,...
            'HalfPowerFrequency2',60*i+5,...
            'DesignMethod','butter','SampleRate',Fs);
        buttLoop = filtfilt(d,buttLoop);
    end
    FAP(i,:) = buttLoop';
end

APd = FAP(:,1:4:end);
clear RAW muafilt S APstream dum1 dum2 FAP
ns = find(APt<=1,1,'last');
TRIALS=[];

for i=1:length(TrOn)
    tr0 = (i);
    si = find(APt>=tr0,1,'first');
    TRIALS = cat(3,TRIALS,APd(:,si:si+ns));
end

%size(TRIALS)
avgT = mean(TRIALS,3);

%%
figure; hold on;
for c= 1:length(ch)
    plot(APt(1:1+ns)+tw(1),avgT(c,:)-mean(avgT(c,:)));
end


%% put data back in the original file structure
rmpath(genpath(fullfile(data_root_probe,ProbeLabel)));
movefile(fullfile(data_root_probe,'experiment*'),data_root_org);
rmdir(data_root_probe)