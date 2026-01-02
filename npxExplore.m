%% add tools
scriptfld = '/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels';

openephys_fld = ['/media/DOCUMENTS/DOCUMENTS/EPHYS_ANALYSIS/NHP-Neuropixels/'...
    'OpenEphys/open-ephys-matlab-tools'];
addpath(genpath(openephys_fld));

%% set paths
data_fld = '/media/NETDISKS/VS03/VS03_6/Neuropixels_NHP/Data_collection';
%data_fld = '/media/chris/CK4TB/Neuropixels_NHP/Data_collection/';
subject = 'm032';
day = '20251022';
expt = 1;
rec = 5;

curr_data = fullfile(data_fld,subject,day);
oebin_path = fullfile(data_fld,subject,day,['experiment',num2str(expt)],['recording',num2str(rec)]);
recnode = fullfile(scriptfld,'Record Note 101');
recnode2 = fullfile(scriptfld,'Record\ Note\ 101'); % need to deal with spaces

%% load data [session method]
% temporarily softlink the content of the session to a folder named 
% 'Record Note 101' so that the SDK can read it.
[~,~] = mkdir(recnode);
system(sprintf('ln -s %s/* %s/', curr_data, recnode2));


% 
nodes = dir(fullfile(curr_data, 'Record Node*'));
sess = Session(curr_data) 


%% 
% path = folder containing 'structure.oebin'
data = load_open_ephys_binary(oebin_path, 'continuous', 1, 'mmap');

% Now you can access your streams directly by name
ap_data = rec.continuous('ProbeA-AP'); 
fprintf('Loaded %d channels\n', ap_data.info.num_channels);


%% 
% 1. Parse the JSON manually (bypass the broken library logic)
json_file = fileread(oebin_path);
metadata = jsondecode(json_file);

% 2. Get the specific folder name for your AP stream
% Usually 'Neuropix-PXI-100.ProbeA-AP' or similar
stream_folder = metadata.continuous(1).folder_name;
dat_path = fullfile(fileparts(oebin_path), 'continuous', stream_folder, 'continuous.dat');

% 3. Map the data directly
num_ch = metadata.continuous(1).num_channels;
file_info = dir(dat_path);
num_samples = file_info.bytes / (num_ch * 2); % 2 bytes per int16 sample

m = memmapfile(dat_path, 'Format', {'int16', [num_ch, num_samples], 'Data'});

% Access a 1-second snippet (assuming 30kHz)
snippet = m.Data.Data(:, 1:30000);