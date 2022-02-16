%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             PATHS & NAMES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

brainstorm3 = '/shares/tfawcett_group/amohammed2/ASSL_scripts/brainstorm3';
brainstorm_db = '/shares/tfawcett_group/amohammed2/ASSL_scripts/brainstorm_db';
raw_data_path = '/shares/tfawcett_group/amohammed2/ASSL_scripts/P01/test_db';
channel_path = '/shares/tfawcett_group/amohammed2/ASSL_scripts/P01/ANT_Channels.txt';

ProtocolName = 'test_db_16';

events = {'0001,0002'}; % Events to be imported Example: {'1,2'} NOT: {'1','2'}
new_names = {'Event_1','Event_2'}; % New event names Example: {'New1','New2'} NOT: {'New1, New2'}

image_folder_name = 'images'; % Name for folder where screenshots will be stored (it'll be created if it DNE)
image_folder_path = 1; % Value of 1 for default path (data folder of protocol)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empty arrays [] are the equivelant of choosing 'all' in BS gui

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Link raw data
link_raw_data.fcn = @link_raw_data_fcn;
link_raw_data.output = 'raw';
link_raw_data.file_type = 'EEG-ANT-CNT';
link_raw_data.image_preview = 1;
link_raw_data.image_file_tag = '1_raw_';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Add EEG positions from channel file
eeg_positions.fcn = @eeg_positions_fcn;
eeg_positions.channel_file_type = 'ASCII_NXYZ';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Process notch filter 
notch.fcn = @notch_fcn;
notch.input = 'raw';
notch.output = 'notch';
notch.sensor_types = 'EEG';
notch.freq_list = 60;
notch.cut_off_W = 1;
notch.use_old = 0;
notch.read_all = 0; 
notch.image_preview = 1;
notch.image_file_tag = '2_notch_';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Process band pass filter
band_pass.fcn = @band_pass_fcn;
band_pass.input = 'notch';
band_pass.output = 'bandpass';
band_pass.band_pass_sensor_types = 'EEG';
band_pass.lower_cutoff = 0.1; 
band_pass.upper_cutoff = 100; 
band_pass.tran_band = 0;
band_pass.attenuation = 'strict'; % 'strict' for 60dB 'relax' for 40dB
band_pass.version = '2019';
band_pass.process_entire_file_all_at_once = 0;
% Detect eye blinks
band_pass.detect_blinks_channel_name = 'Fp1,Fp2';
band_pass.detect_blinks_time_window = [];
band_pass.detect_blinks_event_name = 'blink';
% SSP EOG
band_pass.SSP_EOG_event_name = 'blink';
band_pass.SSP_EOG_sensor_types = 'EEG';
band_pass.SSP_EOG_use_existing_projectors = 1;
band_pass.SSP_EOG_select = 1;
% Re-reference EEG
band_pass.Re_ref_EEG_ref_channel = 'AVERAGE';
band_pass.Re_ref_sensor_types = 'EEG';
% Preview
band_pass.image_preview = 1;
band_pass.image_file_tag = '3_band_';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Import EEG events options
import_events.fcn = @import_events_fcn;
import_events.input = 'bandpass';
import_events.output = 'raw_events';
import_events.time_window = [];
import_events.epoch_time = [-0.1, 1];
import_events.create_cond = 1;
import_events.ignore_short = 1;
import_events.use_ctf_comp = 1;
import_events.use_ssp = 1;
% DC offset correction
import_events.DC_offset_freq = [];
import_events.DC_offset_baseline = [-0.1, -0.002]; 

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Process detrend
detrend.fcn = @detrend_fcn;
detrend.input = 'raw_events';
detrend.output = 'detrend_events';
detrend.time_window = [-0.1, 1];
detrend.sensor_types = 'EEG';
detrend.overwrite = 0;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Average event
avg_event.fcn = @avg_files_fcn;
avg_event.input = 'detrend_events';
avg_event.output = 'avg_event';
avg_event.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
avg_event.avg_type = 5;       
avg_event.avg_func = 1;     
avg_event.weighted = 0;      
avg_event.keep_events = 0;
avg_event.image_preview = 1;
avg_event.image_file_tag = '4_avg_event_';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Create/copy headmodel
head.fcn = @head_fcn;
head.method = 'isotropic';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Compute noise covarience
noise_average_event.fcn = @noise_fcn;
noise_average_event.identity = 0; % 0 for off, 1 for identity matrix (no noise)
noise_average_event.input = 'avg_event'; 
noise_average_event.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
noise_average_event.baseline = [-0.1, -0.002];
noise_average_event.time_window = [0, 1];
noise_average_event.remove_dcoffset = 1; % Block by block
noise_average_event.replace_file = 1; % Replace file if already exists

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Compute sources
sources.fcn = @sources_fcn;
sources.input = 'avg_event';
sources.output = 'sources';
sources.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
sources.process_option = 2; % 2 = Kernel only: one per file
sources.inverse_method = 'minnorm';
sources.inverse_measure = 'dspm2018';
sources.order = 0.5; % for depth weighting 
sources.maximal_amount = 10; % for depth weighting
sources.noise_cov_reg = 0.1;
sources.reg_par = 3;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

no_noise_avg_event.fcn = @noise_fcn;
no_noise_avg_event.identity = 1; % 0 for off, 1 for identity matrix (no noise)
no_noise_avg_event.input = 'avg_event'; 
no_noise_avg_event.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
no_noise_avg_event.baseline = [-0.1, -0.002];
no_noise_avg_event.time_window = [0, 1];
no_noise_avg_event.remove_dcoffset = 1; % Block by block
no_noise_avg_event.replace_file = 1; % Replace file if already exists

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Fast Fourier Transform
FFT_sen.fcn = @FFT_fcn;

% FFT for sensor-level spectral analysis
FFT_sen.sensor_source = 1; % 1 for sensor, 2 for source
FFT_sen.scouts = 0; % 1 for destrieux atlas, 'name of atlas' for custom atlas, 0 for off
FFT_sen.scout_list = {};
FFT_sen.input = 'avg_event'; % 'sources'
FFT_sen.output = 'FFT_sen';
FFT_sen.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
FFT_sen.avg_output = 1;
FFT_sen.image_preview = 1;
FFT_sen.image_file_tag = '5_FFT_sen_';
% For sensor FFT
FFT_sen.sensor_types = {'EEG'};
FFT_sen.time_window = [];
% For spectral FFT
FFT_sen.scout_func = 1;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

extract_sensor.fcn = @extract_fcn; 
extract_sensor.freq_range = [0, 100];
extract_sensor.input = 'FFT_sen';
extract_sensor.output = 'extract_FFT_sen';
extract_sensor.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
extract_sensor.time_window = [];
extract_sensor.rows = '';
extract_sensor.is_absolute = 1;
extract_sensor.avg_time = 0;
extract_sensor.avg_row = 0;
extract_sensor.avg_freq = 0;
extract_sensor.match_rows = 1;
extract_sensor.dim = 2;  % Concatenate time (dimension 2)
extract_sensor.image_preview = 1;
extract_sensor.image_file_tag = '6_SLSA_extract_';
extract_sensor.image_preview = 1;
extract_sensor.image_file_tag = '6_extract_sensor';



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

avg_sensor.fcn = @avg_signals_fcn;
avg_sensor.input = 'extract_FFT_sen';
avg_sensor.output = 'avg_sen';
avg_sensor.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
avg_sensor.avg_type = 1;       
avg_sensor.avg_func = 1;     
avg_sensor.overwrite = 0;
avg_sensor.image_preview = 1;
avg_sensor.image_file_tag = '7_avg_sen_';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Difference
difference.fcn = @difference_fcn; 
difference.A_dir = 'Event_1';
difference.A_file = 'extract_FFT_sen';
difference.B_dir = 'Event_2';
difference.B_file = 'extract_FFT_sen';
difference.output = 'difference_of_events';
difference.image_preview = 1;
difference.image_file_tag = '8_difference_';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

no_noise_notch.fcn = @noise_fcn;

no_noise_notch.identity = 1; % 0 for noise, 1 for identity matrix (no noise)
no_noise_notch.input = 'notch'; 
no_noise_notch.events = 0; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
no_noise_notch.baseline = [-0.1, -0.002];
no_noise_notch.time_window = [0, 1];
no_noise_notch.remove_dcoffset = 1; % Block by block
no_noise_notch.replace_file = 1; % Replace file if already exists

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

FFT_scout.fcn = @FFT_fcn;

% FFT for spectral output of each scout
FFT_scout.sensor_source = 2; % 1 for sensor, 2 for source
FFT_scout.scouts = 1; % 1 for destrieux atlas, 'name of atlas' for custom atlas, 0 for off
FFT_scout.scout_list = {'L HG_T_transv', 'R HG_T_transv', 'L STG_Lateral', 'R STG_Lateral',...
'L Plan_tempo','R Plan_tempo', 'L Transverse Sulcus','R Transverse Sulcus'};
FFT_scout.input = 'sources';
FFT_scout.output = 'FFT_scout';
FFT_scout.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
FFT_scout.avg_output = 1;
FFT_scout.image_preview = 1;
FFT_scout.image_file_tag = '9_FFT_scout_';
% For spectral FFT
FFT_scout.scout_func = 1;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

extract_scout.fcn = @extract_fcn;
extract_scout.freq_range = [0, 100];
extract_scout.input = 'FFT_scout';
extract_scout.output = 'extract_scout';
extract_scout.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
extract_scout.time_window = [];
extract_scout.rows = '';
extract_scout.is_absolute = 1;
extract_scout.avg_time = 0;
extract_scout.avg_row = 0;
extract_scout.avg_freq = 0;
extract_scout.match_rows = 1;
extract_scout.dim = 2;  % Concatenate time (dimension 2)
extract_scout.image_preview = 1;
extract_scout.image_file_tag = '10_extract_scout_';

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

FFT_source_space.fcn = @FFT_fcn;
% FFT for source space
FFT_source_space.sensor_source = 2; % 1 for sensor, 2 for source
FFT_source_space.scouts = 0; % 1 for destrieux atlas, 'name of atlas' for custom atlas, 0 for off
FFT_source_space.scout_list = {};
FFT_source_space.input = 'sources';
FFT_source_space.output = 'FFT_source';
FFT_source_space.events = 'all'; % 0 for none, 'all' for all, {'Event2','Event5','Event8'} to select specific events
FFT_source_space.avg_output = 1;
% For spectral FFT
FFT_source_space.scout_func = 1;
FFT_source_space.image_preview = 0;
FFT_source_space.image_file_tag = '11_FFT_source_space_';


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

group_tf_bands.fcn = @group_tf_bands_fcn;
group_tf_bands.input = 'FFT_source';
group_tf_bands.output = 'tf_bands';
group_tf_bands.events = 'all';
group_tf_bands.is_freq_bands = 1;
group_tf_bands.freq_bands = {'delta', '2, 4', 'mean'; 'theta', '5, 7', 'mean';...
 'alpha', '8, 12', 'mean'; 'beta', '15, 29', 'mean'; 'ASSR', '38, 40', 'mean';...
 'gamma1', '30, 59', 'mean'; 'gamma2', '60, 90', 'mean'};
group_tf_bands.is_time_bands = 0;
group_tf_bands.time_bands = '';
group_tf_bands.overwrite = 0;
group_tf_bands.image_preview = 0;
group_tf_bands.image_file_tag = '11_FFT_source_space_';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                PROCESSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose the functions you want to run (their order matters)

processes = {link_raw_data,eeg_positions,notch,band_pass,import_events,detrend,...
avg_event,head,noise_average_event,sources,no_noise_avg_event,FFT_sen,...
extract_sensor,avg_sensor,difference,no_noise_notch,FFT_scout,extract_scout,...
FFT_source_space,group_tf_bands};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           ***** Nothing else below is meant to be edited *****
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(brainstorm3);
addpath(brainstorm_db);
addpath(raw_data_path);

% Creates a structure called names which contains all paths and names
names.brainstorm_db = brainstorm_db;
names.ProtocolName = ProtocolName;
names.dir_path = raw_data_path;
names.channel_path = channel_path;
names.events = events;
names.old_event_names = split(events,',')';
names.new_event_names = new_names;
names.image_folder = image_folder_name;
names.db_data_path = fullfile(names.brainstorm_db,names.ProtocolName,'data');
if image_folder_path == 1
names.image_folder_path = names.db_data_path;
else
names.image_folder_path = image_folder_path;
end

% Finds the names of the subjects 
clear sub_name
sub_name = dir(names.dir_path);
sub_name = sub_name([sub_name(:).isdir]);
sub_name = sub_name(~ismember({sub_name(:).name},{'.','..','images'}));
sub_name = struct2cell(sub_name);
sub_name = sub_name(1,:);

% SubjectName = [91 83 117 98 106 101 99 116 78 97 109 101 93];

% Finds the raw data files

for i = 1:length(sub_name)
cnt_dir = dir(fullfile(names.dir_path,sub_name{i},'*.cnt'));
cnt_dir = struct2cell(cnt_dir);
cnt.path{i} = cnt_dir(2,:)';
cnt.file_name{i} = cnt_dir(1,:)';
end

% Converts sub_name into string array with one row
for i = 1:length(sub_name)
sub(i) = string(sub_name(i));
end

for i = 1:length(sub_name)

dates_s = [];
for j = 1:length(cnt.path{i})
id = fopen(fullfile(cnt.path{i}{j},cnt.file_name{i}{j}),'r');
fseek(id,-50000,'eof');
ASC = fread(id);

% StartDate = [91 83 116 97 114 116 68 97 116 101 93];
% StartFraction = [91 83 116 97 114 116 70 114 97 99 116 105 111 110 93];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ASC_string = convertStringsToChars(strjoin(string(ASC')));
SubjectName_string = "91 83 117 98 106 101 99 116 78 97 109 101 93*10 91";
SubjectName_string = regexptranslate('wildcard',SubjectName_string);

[start_loc,end_loc] = regexp(ASC_string, SubjectName_string);

ASC_string = ASC_string(start_loc:end_loc);

start_loc = regexp(ASC_string, "91 83 117 98 106 101 99 116 78 97 109 101 93",'end','once');
end_loc = regexp(ASC_string, "10 91",'start','once');

name_ascii = ASC_string((start_loc+5):(end_loc-2));

name_ascii = str2num(name_ascii);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ASC_string = convertStringsToChars(strjoin(string(ASC')));
StartDate_string = "91 83 116 97 114 116 68 97 116 101 93*10 91";
StartDate_string = regexptranslate('wildcard',StartDate_string);

[start_loc,end_loc] = regexp(ASC_string, StartDate_string);

ASC_string = ASC_string(start_loc:end_loc);

start_loc = regexp(ASC_string, "91 83 116 97 114 116 68 97 116 101 93",'end','once');
end_loc = regexp(ASC_string, "10 91",'start','once');

date_ascii = ASC_string((start_loc+5):(end_loc-2));

date_ascii = str2num(date_ascii);



% loc_d = strfind(ASC',StartDate);
% loc_f = strfind(ASC',StartFraction);
% d_ascii = ASC((loc_d+12):(loc_f-2));
fclose(id);


d_serial = cast(date_ascii,'char');
d_serial = extractBefore(d_serial,'.');
dates_s{i}(j) = x2mdate(str2num(d_serial));
weeks{i}{j} = datestr(dates_s{i}(j));
end


dates{i} = weeks{i};
first = find(min(dates_s{i}));
weeks{i}{first} = 'Day_0';
first_week = week(datetime(dates{i}{first}));

for k = 1:length(weeks{i})

  if strcmp(weeks{i}{k},weeks{i}{first})
  continue
  end
 
  year_week{i}(k) = week(datetime(dates{i}{k}));
  week_num{i}(k) = year_week{i}(k)-first_week;

  if week_num{i}(k) < 0
  week_num{i}(k) = floor((dates_s{i}(k)- dates_s{i}(first))/7);
  end

  weeks{i}{k} = strcat('Week_',num2str(week_num{i}(k)));

end

if ~exist('subject_name','var')
subject_name = cell(1,length(sub_name));
end

for l = 1:length(cnt.path{i})

check_sn = strcat(sub{i}, '_', weeks{i}{l},'_',char(datetime(dates{i}{l},'format','MMM_d_yyyy')));

if ~any(strcmp(subject_name{i},check_sn))
subject_name{i}{l} = strcat(sub{i}, '_', weeks{i}{l},'_',char(datetime(dates{i}{l},'format','MMM_d_yyyy')));
else
t = 1;
w = 2;
while t > 0
new_attempt = strcat(check_sn,'_',num2str(w));
t = any(strcmp(subject_name{i},new_attempt));
w = w + 1;
end
subject_name{i}{l} = new_attempt;
end

end
end

cnt.path = vertcat(cnt.path{:});
cnt.file_name = vertcat(cnt.file_name{:});

names.sub_name = horzcat(subject_name{:})';
names.raw_path = fullfile(cnt.path,cnt.file_name);
names.file_name = cnt.file_name;
names.file_name = strrep(names.file_name,'.cnt','');

for i = 1:length(names.sub_name)
db.sub{i} = [];
for j = 1:length(names.old_event_names)
db.sub{i}.events{j} = [];
end
end

% load(strcat(ProtocolName,'_db.mat'))

% Start BS
 brainstorm % To run without gui, use: brainstorm sever
%%
% Create/load protocol
iProtocol = bst_get('Protocol', ProtocolName);
if ~isempty(iProtocol)
gui_brainstorm('SetCurrentProtocol', iProtocol);
elseif isempty(iProtocol)
gui_brainstorm('CreateProtocol', ProtocolName, 1, 0); 
end

% This is for the log file
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
d = datetime(clock, 'InputFormat', 'uuuu-DDD''T''HH:mm:ss.SSS');
fprintf(file_ID, '%s - Protocol: %s \n', d, names.ProtocolName);
fclose(file_ID);


%%%%%%%%%%%%%%%%%%%%%%%% Runs each chosen process %%%%%%%%%%%%%%%%%%%%%%%%%
%%



for p = 1:length(processes)

db = processes{p}.fcn(processes{p},names,db);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = link_raw_data_fcn(in, names, db)

for i = 1:length(names.sub_name)

tic

sFiles = [];
SubjectNames = {names.sub_name{i}};
RawFiles = {names.raw_path{i}};

% Process: Create link to raw file
sFiles = bst_process('CallProcess', 'process_import_data_raw', sFiles, [], ...
    'subjectname',    SubjectNames{1}, ...
    'datafile',       {RawFiles{1}, in.file_type}, ...
    'channelreplace', 1, ...
    'channelalign',   1, ...
    'evtmode',        'value');

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

dir_name = strcat('@raw',names.file_name{i});
rename_a = fullfile(names.sub_name{i}, dir_name);
new_name = strcat(names.sub_name{i},'_raw');
rename_b = fullfile(names.sub_name{i}, new_name);
db_rename_condition(rename_a, rename_b);

raw = sFiles;

raw.Condition = new_name;
splits = split(raw.FileName, '/');
splits{2} = new_name;
raw.FileName = fullfile(splits{1},splits{2},splits{3});
splits = split(raw.ChannelFile, '/');
splits{2} = new_name;
raw.ChannelFile = fullfile(splits{1},splits{2},splits{3});

db.sub{i}.(in.output) = raw;

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - Raw data linked\n',...
 floor(t/60), rem(t,60),names.sub_name{i});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - Raw data linked\n',...
 floor(t/60), rem(t,60),names.sub_name{i});
fclose(file_ID);
end

end

save(strcat(names.ProtocolName,'_db.mat'),'db')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = eeg_positions_fcn(in, names, db)

for i = 1:length(names.sub_name)

tic

  raw_dir = strcat(names.sub_name{i},'_raw');
  raw_link_path = fullfile(names.brainstorm_db,names.ProtocolName,'data',...
    names.sub_name{i},raw_dir,'data_0raw*');

Data_dir = dir(raw_link_path);
Data_name = {Data_dir.name};
path = fullfile(names.sub_name{i},raw_dir,Data_name);

sFiles = path;
RawFiles = {names.channel_path};

% Process: Add EEG positions
sFiles = bst_process('CallProcess', 'process_channel_addloc', sFiles, [], ...
    'channelfile', {RawFiles{1}, in.channel_file_type}, ...
    'usedefault',  1, ...  % 1 means dont use default 
    'fixunits',    0, ...
    'vox2ras',     0);

t = toc;

if i ==1 
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - EEG positions\n',...
 floor(t/60), rem(t,60),names.sub_name{i});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - EEG positions\n',...
 floor(t/60), rem(t,60),names.sub_name{i});
fclose(file_ID);
end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = notch_fcn(in, names, db)

for i = 1:length(names.sub_name)

tic

sFiles = db.sub{i}.(in.input);
        
% Process: Notch filter: 60Hz
sFiles = bst_process('CallProcess', 'process_notch', sFiles, [], ...
    'sensortypes', in.sensor_types, ...
    'freqlist',    in.freq_list, ...
    'cutoffW',     in.cut_off_W, ...
    'useold',      in.use_old, ...
    'read_all',    in.read_all);

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

dir_name = strcat('@raw',names.file_name{i},'_notch');
rename_a = fullfile(names.sub_name{i}, dir_name);
new_name = strcat(names.sub_name{i},'_notch');
rename_b = fullfile(names.sub_name{i}, new_name);

db_rename_condition(rename_a, rename_b);

notch = sFiles;

notch.Condition = new_name;
splits = split(notch.FileName, '/');
splits{2} = new_name;
notch.FileName = fullfile(splits{1},splits{2},splits{3});
splits = split(notch.ChannelFile, '/');
splits{2} = new_name;
notch.ChannelFile = fullfile(splits{1},splits{2},splits{3});

db.sub{i}.(in.output) = notch;

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - Notch filter processed\n',...
floor(t/60), rem(t,60), names.sub_name{i});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - Notch filter processed\n',...
floor(t/60), rem(t,60), names.sub_name{i});
fclose(file_ID);
end

end

save(strcat(names.ProtocolName,'_db.mat'),'db')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = band_pass_fcn(in, names, db)

for i = 1:length(names.sub_name)

tic

sFiles = db.sub{i}.(in.input);

% Process: Band-pass:0.1Hz-100Hz
sFiles = bst_process('CallProcess', 'process_bandpass', sFiles, [], ...
    'sensortypes', in.band_pass_sensor_types, ...
    'highpass',    in.lower_cutoff, ... % This is the lower cutoff frequency
    'lowpass',     in.upper_cutoff, ... % This is the upper cutoff frequency
    'tranband',    in.tran_band, ...
    'attenuation', in.attenuation, ...  
    'ver',         in.version, ...  
    'mirror',      0, ...
    'read_all',    in.process_entire_file_all_at_once);

% Process: Detect eye blinks
sFiles = bst_process('CallProcess', 'process_evt_detect_eog', sFiles, [], ...
    'channelname', in.detect_blinks_channel_name, ...
    'timewindow',  in.detect_blinks_time_window, ...
    'eventname',   in.detect_blinks_event_name);

% Process: SSP EOG: blink
sFiles = bst_process('CallProcess', 'process_ssp_eog', sFiles, [], ...
    'eventname',   in.SSP_EOG_event_name, ...
    'sensortypes', in.SSP_EOG_sensor_types, ...
    'usessp',      in.SSP_EOG_use_existing_projectors , ...
    'select',      in.SSP_EOG_select);

% Process: Re-reference EEG
sFiles = bst_process('CallProcess', 'process_eegref', sFiles, [], ...
    'eegref',      in.Re_ref_EEG_ref_channel, ...
    'sensortypes', in.Re_ref_sensor_types);

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

dir_name = strcat('@raw',names.file_name{i},'_notch_band');
rename_a = fullfile(names.sub_name{i}, dir_name);
new_name = strcat(names.sub_name{i},'_bandpass');
rename_b = fullfile(names.sub_name{i}, new_name);

db_rename_condition(rename_a, rename_b);

bandpass = sFiles;

bandpass.Condition = new_name;
splits = split(bandpass.FileName, '/');
splits{2} = new_name;
bandpass.FileName = fullfile(splits{1},splits{2},splits{3});
splits = split(bandpass.ChannelFile, '/');
splits{2} = new_name;
bandpass.ChannelFile = fullfile(splits{1},splits{2},splits{3});

db.sub{i}.(in.output) = bandpass;

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID,'\n%d min %.f sec - %s - Band-pass filter processed & eye blinks removed\n',...
floor(t/60), rem(t,60), names.sub_name{i});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID,'%d min %.f sec - %s - Band-pass filter processed & eye blinks removed\n',...
floor(t/60), rem(t,60), names.sub_name{i});
fclose(file_ID);
end

end

save(strcat(names.ProtocolName,'_db.mat'),'db')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = import_events_fcn(in, names, db)

for i = 1:length(names.sub_name)

tic  

sFiles = db.sub{i}.(in.input);

SubjectNames = {names.sub_name{i}};

% Import data event
% Process: Import MEG/EEG: Events
sFiles = bst_process('CallProcess', 'process_import_data_event', sFiles, [], ...
    'subjectname', SubjectNames{1}, ...
    'condition',   '', ...
    'eventname',   names.events{1}, ...
    'timewindow',  in.time_window, ...
    'epochtime',   in.epoch_time, ...
    'createcond',  in.create_cond, ...
    'ignoreshort', in.ignore_short, ...
    'usectfcomp',  in.use_ctf_comp, ...
    'usessp',      in.use_ssp, ...
    'freq',        in.DC_offset_freq, ... % for DC offset correction
    'baseline',    in.DC_offset_baseline); % DC offset correction

raw_events = sFiles;

for j = 1:length(names.old_event_names)

rename_a = fullfile(names.sub_name{i}, names.old_event_names{j});
rename_b = fullfile(names.sub_name{i}, names.new_event_names{j});

check_a = fullfile(names.db_data_path,rename_a);

if isempty(dir(check_a))
continue
end

db_rename_condition(rename_a, rename_b);

index = find(strcmp({raw_events.Condition}, names.old_event_names{j})==1);
first = index(1);
last = index(length(index));
db.sub{i}.events{j}.(in.output) = sFiles(first:last);


for y = 1:length(db.sub{i}.events{j}.(in.output))
cond = db.sub{i}.events{j}.(in.output)(y).Condition;
old = names.old_event_names{j};
if cond(1) == old(1)
db.sub{i}.events{j}.(in.output)(y).Condition = names.new_event_names{j};
splits = split(db.sub{i}.events{j}.(in.output)(y).FileName, '/');
splits{2} = names.new_event_names{j};
db.sub{i}.events{j}.(in.output)(y).FileName = fullfile(splits{1},splits{2},splits{3});
splits = split(db.sub{i}.events{j}.(in.output)(y).ChannelFile, '/');
splits{2} = names.new_event_names{j};
db.sub{i}.events{j}.(in.output)(y).ChannelFile = fullfile(splits{1},splits{2},splits{3});
else
continue
end
end
end
end



t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - All events imported\n',...
floor(t/60), rem(t,60), names.sub_name{i});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - All events imported\n',...
floor(t/60), rem(t,60), names.sub_name{i});
fclose(file_ID);
end

save(strcat(names.ProtocolName,'_db.mat'),'db')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = detrend_fcn(in, names,db)

for i = 1:length(names.sub_name)

for j = 1:length(names.new_event_names)

tic 

sFiles = db.sub{i}.events{j}.(in.input);

% Process detrend
% Process: Remove linear trend
sFiles = bst_process('CallProcess', 'process_detrend', sFiles, [], ...
    'timewindow',  in.time_window, ...
    'sensortypes', in.sensor_types, ...
    'overwrite',   in.overwrite);

db.sub{i}.events{j}.(in.output) = sFiles;

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Detrend Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Detrend Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},names.new_event_names{j});
fclose(file_ID);
end

end

end

save(strcat(names.ProtocolName,'_db.mat'),'db')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = avg_files_fcn(in, names,db)

if strcmp(in.events, 'all')
for i = 1:length(names.sub_name)

for j = 1:length(names.new_event_names)

tic 

sFiles = db.sub{i}.events{j}.(in.input);

% Process: Average: By trial group (folder average)
sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
    'avgtype',       in.avg_type, ...  
    'avg_func',      in.avg_func, ... 
    'weighted',      in.weighted, ...
    'keepevents',    in.keep_events);

db.sub{i}.events{j}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Average processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Average processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},names.new_event_names{j});
fclose(file_ID);
end

end

end

elseif iscell(in.events)
for i = 1:length(names.sub_name)

for j = 1:length(in.events)

tic

loc = find(strcmp(names.new_event_names, in.events{j}));

sFiles = db.sub{i}.events{loc}.(in.input);

% Process: Average: By trial group (folder average)
sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
    'avgtype',       in.avg_type, ...  % By trial group (folder average)
    'avg_func',      in.avg_func, ...  % Arithmetic average:  mean(x)
    'weighted',      in.weighted, ...
    'keepevents',    in.keep_events);

db.sub{i}.events{loc}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Average processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Average processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},in.events{j});
fclose(file_ID);
end

end

end

elseif in.events == 0
for i = 1:length(names.sub_name)

tic 

sFiles = db.sub{i}.(in.input);

% Process: Average: By trial group (folder average)
sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
    'avgtype',       in.avg_type, ...  % By trial group (folder average)
    'avg_func',      in.avg_func, ...  % Arithmetic average:  mean(x)
    'weighted',      in.weighted, ...
    'keepevents',    in.keep_events);

db.sub{i}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - Average processed for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.input);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - Average processed for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i},in.input);
fclose(file_ID);
end

end
end

save(strcat(names.ProtocolName,'_db.mat'),'db')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = avg_signals_fcn(in, names,db)

if strcmp(in.events, 'all')
for i = 1:length(names.sub_name)

for j = 1:length(names.new_event_names)

tic 

sFiles = db.sub{i}.events{j}.(in.input);

% Process: Average: All signals
sFiles = bst_process('CallProcess', 'process_average_rows', sFiles, [], ...
    'avgtype',   in.avg_type, ...  
    'avgfunc',   in.avg_func, ...  
    'overwrite', in.overwrite);

db.sub{i}.events{j}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Average processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Average processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},names.new_event_names{j});
fclose(file_ID);
end

end

end

elseif iscell(in.events)
for i = 1:length(names.sub_name)

for j = 1:length(in.events)

tic

loc = find(strcmp(names.new_event_names, in.events{j}));

sFiles = db.sub{i}.events{loc}.(in.input);

% Process: Average: All signals
sFiles = bst_process('CallProcess', 'process_average_rows', sFiles, [], ...
    'avgtype',   in.avg_type, ...  
    'avgfunc',   in.avg_func, ...  
    'overwrite', in.overwrite);

db.sub{i}.events{loc}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Average processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Average processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},in.events{j});
fclose(file_ID);
end

end

end

elseif in.events == 0
for i = 1:length(names.sub_name)

tic 

sFiles = db.sub{i}.(in.input);

% Process: Average: All signals
sFiles = bst_process('CallProcess', 'process_average_rows', sFiles, [], ...
    'avgtype',   in.avg_type, ...  
    'avgfunc',   in.avg_func, ...  
    'overwrite', in.overwrite);

db.sub{i}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - Average processed for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.input);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - Average processed for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i},in.input);
fclose(file_ID);
end

end
end

save(strcat(names.ProtocolName,'_db.mat'),'db')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = head_fcn(in, names, db)

tic

raw_dir = strcat(names.sub_name{1},'_raw');
raw_link_path = fullfile(names.brainstorm_db,names.ProtocolName,'data',...
    names.sub_name{1},raw_dir,'data_0raw*');

Data_dir = dir(raw_link_path);
Data_name = {Data_dir.name};
path = fullfile(names.sub_name{1},raw_dir,Data_name);

sFiles = path;
 
% Process: Compute head model
sFiles = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
    'Comment',     '', ...
    'sourcespace', 1, ...  % Cortex surface
    'volumegrid',  struct(...
         'Method',        in.method, ...
         'nLayers',       17, ...
         'Reduction',     3, ...
         'nVerticesInit', 4000, ...
         'Resolution',    0.005, ...
         'FileName',      ''), ...
    'meg',         3, ...  % Overlapping spheres
    'eeg',         3, ...  % OpenMEEG BEM
    'ecog',        2, ...  % OpenMEEG BEM
    'seeg',        2, ...  % OpenMEEG BEM
    'openmeeg',    struct(...
         'BemFiles',     {{}}, ...
         'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
         'BemCond',      [1, 0.0125, 1], ...
         'BemSelect',    [1, 1, 1], ...
         'isAdjoint',    0, ...
         'isAdaptative', 1, ...
         'isSplit',      0, ...
         'SplitLength',  4000));

t = toc;

log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - Head model created for first subject\n',...
floor(t/60), rem(t,60));
fclose(file_ID);

tic

headmodel_path = fullfile(names.brainstorm_db,names.ProtocolName,'data',...
    names.sub_name{1},raw_dir,'headmodel*');

headmodel_dir = dir(headmodel_path);
headmodel_name = {headmodel_dir.name}; 

first_headmodel_path = fullfile(names.sub_name{1},raw_dir,headmodel_name);

db_set_headmodel(first_headmodel_path{1}, 'AllSubjects');

t = toc; 

log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - Head model copied to all subjects\n',...
floor(t/60), rem(t,60));
fclose(file_ID);
db=db;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = noise_fcn(in, names, db)

if strcmp(in.events, 'all') 
for i = 1:length(names.sub_name)

for j = 1:length(names.new_event_names)

tic

sFiles = db.sub{i}.events{j}.(in.input);

 % Process: Compute covariance (noise or data)
sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
    'baseline',       in.baseline, ...
    'datatimewindow', in.time_window, ...
    'sensortypes',    'EEG', ...
    'target',         1, ...  % Noise covariance     (covariance over baseline time window)
    'dcoffset',       in.remove_dcoffset, ...  % Block by block, to avoid effects of slow shifts in data
    'identity',       in.identity, ...
    'copycond',       0, ...
    'copysubj',       0, ...
    'copymatch',      0, ...
    'replacefile',    in.replace_file);  % Replace

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
if in.identity == 0
fprintf(file_ID, '\n%d min %.f sec - %s - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j}, in.input);
fclose(file_ID);
elseif in.identity == 1
fprintf(file_ID, '\n%d min %.f sec - %s - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j}, in.input);
fclose(file_ID);
end
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');

if in.identity == 0
fprintf(file_ID, '%d min %.f sec - %s - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j}, in.input);
fclose(file_ID);
elseif in.identity == 1
fprintf(file_ID, '%d min %.f sec - %s - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j}, in.input);
fclose(file_ID);
end
end
end
end

elseif iscell(in.events)
for i = 1:length(names.sub_name)

for j = 1:length(in.events)

tic

loc = find(strcmp(names.new_event_names, in.events{j}));

sFiles = db.sub{i}.events{loc}.(in.input);

 % Process: Compute covariance (noise or data)
sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
    'baseline',       in.baseline, ...
    'datatimewindow', in.time_window, ...
    'sensortypes',    'EEG', ...
    'target',         1, ...  % Noise covariance     (covariance over baseline time window)
    'dcoffset',       in.remove_dcoffset, ...  % Block by block, to avoid effects of slow shifts in data
    'identity',       in.identity, ...
    'copycond',       0, ...
    'copysubj',       0, ...
    'copymatch',      0, ...
    'replacefile',    in.replace_file);  % Replace

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
if in.identity == 0
fprintf(file_ID, '\n%d min %.f sec - %s - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.input);
fclose(file_ID);
elseif in.identity == 1
fprintf(file_ID, '\n%d min %.f sec - %s - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.input);
fclose(file_ID);
end
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');

if in.identity == 0
fprintf(file_ID, '%d min %.f sec - %s - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.input);
fclose(file_ID);
elseif in.identity == 1
fprintf(file_ID, '%d min %.f sec - %s - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.input);
fclose(file_ID);
end
end
end
end

elseif in.events == 0
for i = 1: length(names.sub_name)

tic

sFiles = db.sub{i}.(in.input);

 % Process: Compute covariance (noise or data)
sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
    'baseline',       in.baseline, ...
    'datatimewindow', in.time_window, ...
    'sensortypes',    'EEG', ...
    'target',         1, ...  % Noise covariance     (covariance over baseline time window)
    'dcoffset',       in.remove_dcoffset, ...  % Block by block, to avoid effects of slow shifts in data
    'identity',       in.identity, ...
    'copycond',       0, ...
    'copysubj',       0, ...
    'copymatch',      0, ...
    'replacefile',    in.replace_file);  % Replace

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
if in.identity == 0
fprintf(file_ID, '\n%d min %.f sec - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.input);
fclose(file_ID);
elseif in.identity == 1
fprintf(file_ID, '\n%d min %.f sec - %s - No noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.input);
fclose(file_ID);
end

else

log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
if in.identity == 0
fprintf(file_ID, '%d min %.f sec - %s - Noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.input);
fclose(file_ID);
elseif in.identity == 1
fprintf(file_ID, '%d min %.f sec - %s - No noise set for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.input);
fclose(file_ID);
end
end
end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = sources_fcn(in, names, db)

if strcmp(in.events, 'all') 
for i = 1:length(names.sub_name)

for j = 1:length(names.new_event_names)

tic

sFiles = db.sub{i}.events{j}.(in.input);

% Process: Compute sources [2018]
sFiles = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
    'output',  in.process_option, ...  % Kernel only: one per file
    'inverse', struct(...
         'Comment',        'dSPM-unscaled: EEG', ...
         'InverseMethod',  in.inverse_method, ...
         'InverseMeasure', in.inverse_measure, ...
         'SourceOrient',   {{'fixed'}}, ...
         'Loose',          0.2, ...
         'UseDepth',       1, ...
         'WeightExp',      in.order, ...
         'WeightLimit',    in.maximal_amount, ...
         'NoiseMethod',    'reg', ...
         'NoiseReg',       in.noise_cov_reg, ...
         'SnrMethod',      'fixed', ...
         'SnrRms',         0, ...
         'SnrFixed',       in.reg_par, ...
         'ComputeKernel',  1, ...
         'DataTypes',      {{'EEG'}}));

db.sub{i}.events{j}.(in.output) = sFiles;

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Sources Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Sources Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},names.new_event_names{j});
fclose(file_ID);
end

end

end

elseif iscell(in.events)
for i = 1:length(names.sub_name)

for j = 1:length(in.events)

tic

loc = find(strcmp(names.new_event_names, in.events{j}));

sFiles = db.sub{i}.events{loc}.(in.input);

% Process: Compute sources [2018]
sFiles = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
    'output',  in.process_option, ...  % Kernel only: one per file
    'inverse', struct(...
         'Comment',        'dSPM-unscaled: EEG', ...
         'InverseMethod',  in.inverse_method, ...
         'InverseMeasure', in.inverse_measure, ...
         'SourceOrient',   {{'fixed'}}, ...
         'Loose',          0.2, ...
         'UseDepth',       1, ...
         'WeightExp',      in.order, ...
         'WeightLimit',    in.maximal_amount, ...
         'NoiseMethod',    'reg', ...
         'NoiseReg',       in.noise_cov_reg, ...
         'SnrMethod',      'fixed', ...
         'SnrRms',         0, ...
         'SnrFixed',       in.reg_par, ...
         'ComputeKernel',  1, ...
         'DataTypes',      {{'EEG'}}));

db.sub{i}.events{loc}.(in.output) = sFiles;

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Sources Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Sources Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},in.events{j});
fclose(file_ID);
end

end

end

elseif in.events == 0
for i = 1:length(names.sub_name)

tic

sFiles = db.sub{i}.(in.input);

% Process: Compute sources [2018]
sFiles = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
    'output',  in.process_option, ...  % Kernel only: one per file
    'inverse', struct(...
         'Comment',        'dSPM-unscaled: EEG', ...
         'InverseMethod',  in.inverse_method, ...
         'InverseMeasure', in.inverse_measure, ...
         'SourceOrient',   {{'fixed'}}, ...
         'Loose',          0.2, ...
         'UseDepth',       1, ...
         'WeightExp',      in.order, ...
         'WeightLimit',    in.maximal_amount, ...
         'NoiseMethod',    'reg', ...
         'NoiseReg',       in.noise_cov_reg, ...
         'SnrMethod',      'fixed', ...
         'SnrRms',         0, ...
         'SnrFixed',       in.reg_par, ...
         'ComputeKernel',  1, ...
         'DataTypes',      {{'EEG'}}));

db.sub{i}.(in.output) = sFiles;

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - Sources Processed for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.input);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - Sources Processed for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i},in.input);
fclose(file_ID);
end

end

end
save(strcat(names.ProtocolName,'_db.mat'),'db')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = FFT_fcn(in, names, db)
tic 

if in.sensor_source == 1

if strcmp(in.events, 'all') 
for i = 1:length(names.sub_name)

for j = 1:length(names.new_event_names)

sFiles = db.sub{i}.events{j}.(in.input);

% Process: Fourier transform (FFT)
sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'timewindow',  in.time_window, ...
    'sensortypes', in.sensor_types, ...
    'avgoutput',   in.avg_output);

db.sub{i}.events{j}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},names.new_event_names{j}, in.output);
fclose(file_ID);
end

end

end
elseif iscell(in.events)
for i = 1:length(names.sub_name)

for j = 1:length(in.events)

tic

loc = find(strcmp(names.new_event_names, in.events{j}));

sFiles = db.sub{i}.events{loc}.(in.input);

% Process: Fourier transform (FFT)
sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'timewindow',  in.time_window, ...
    'sensortypes', in.sensor_types, ...
    'avgoutput',   in.avg_output);

db.sub{i}.events{loc}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.output);
fclose(file_ID);
end

end

end
elseif in.events == 0
for i = 1:length(names.sub_name)

tic

sFiles = db.sub{i}.(in.input);

% Process: Fourier transform (FFT)
sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'timewindow',  in.time_window, ...
    'sensortypes', in.sensor_types, ...
    'avgoutput',   in.avg_output);

db.sub{i}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.output);
fclose(file_ID);
end

end
end

elseif in.sensor_source == 2 

if strcmp(in.events, 'all') 
for i = 1:length(names.sub_name)

for j = 1:length(names.new_event_names)

sFiles = db.sub{i}.events{j}.(in.input);

if ischar(in.scouts)

% User scouts (from the scout_Destrieux_Auditory file)
sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'clusters',  {in.scouts, in.scout_list}, ...
    'scoutfunc', in.scout_func, ...  % Mean
    'avgoutput', in.avg_output);

elseif in.scouts == 1

% Process: Fourier transform (FFT)
sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'clusters',  {'Destrieux', {'G_Ins_lg_and_S_cent_ins L', 'G_Ins_lg_and_S_cent_ins R',...
 'G_and_S_cingul-Ant L', 'G_and_S_cingul-Ant R', 'G_and_S_cingul-Mid-Ant L',...
 'G_and_S_cingul-Mid-Ant R', 'G_and_S_cingul-Mid-Post L', 'G_and_S_cingul-Mid-Post R',...
 'G_and_S_frontomargin L', 'G_and_S_frontomargin R', 'G_and_S_occipital_inf L',...
 'G_and_S_occipital_inf R', 'G_and_S_paracentral L', 'G_and_S_paracentral R',...
 'G_and_S_subcentral L', 'G_and_S_subcentral R', 'G_and_S_transv_frontopol L',...
 'G_and_S_transv_frontopol R', 'G_cingul-Post-dorsal L', 'G_cingul-Post-dorsal R',...
 'G_cingul-Post-ventral L', 'G_cingul-Post-ventral R', 'G_cuneus L', 'G_cuneus R',...
 'G_front_inf-Opercular L', 'G_front_inf-Opercular R', 'G_front_inf-Orbital L',...
 'G_front_inf-Orbital R', 'G_front_inf-Triangul L', 'G_front_inf-Triangul R',...
 'G_front_middle L', 'G_front_middle R', 'G_front_sup L', 'G_front_sup R',...
 'G_insular_short L', 'G_insular_short R', 'G_oc-temp_lat-fusifor L',...
 'G_oc-temp_lat-fusifor R', 'G_oc-temp_med-Lingual L', 'G_oc-temp_med-Lingual R',...
 'G_oc-temp_med-Parahip L', 'G_oc-temp_med-Parahip R', 'G_occipital_middle L',...
 'G_occipital_middle R', 'G_occipital_sup L', 'G_occipital_sup R', 'G_orbital L',...
 'G_orbital R', 'G_pariet_inf-Angular L', 'G_pariet_inf-Angular R',...
 'G_pariet_inf-Supramar L', 'G_pariet_inf-Supramar R', 'G_parietal_sup L',...
 'G_parietal_sup R', 'G_postcentral L', 'G_postcentral R', 'G_precentral L',...
 'G_precentral R', 'G_precuneus L', 'G_precuneus R', 'G_rectus L', 'G_rectus R',...
 'G_subcallosal L', 'G_subcallosal R', 'G_temp_sup-G_T_transv L', 'G_temp_sup-G_T_transv R',...
 'G_temp_sup-Lateral L', 'G_temp_sup-Lateral R', 'G_temp_sup-Plan_polar L',...
 'G_temp_sup-Plan_polar R', 'G_temp_sup-Plan_tempo L', 'G_temp_sup-Plan_tempo R',...
 'G_temporal_inf L', 'G_temporal_inf R', 'G_temporal_middle L', 'G_temporal_middle R',...
 'Lat_Fis-ant-Horizont L', 'Lat_Fis-ant-Horizont R', 'Lat_Fis-ant-Vertical L',...
 'Lat_Fis-ant-Vertical R', 'Lat_Fis-post L', 'Lat_Fis-post R', 'Pole_occipital L',...
 'Pole_occipital R', 'Pole_temporal L', 'Pole_temporal R', 'S_calcarine L', 'S_calcarine R',...
 'S_central L', 'S_central R', 'S_cingul-Marginalis L', 'S_cingul-Marginalis R',...
 'S_circular_insula_ant L', 'S_circular_insula_ant R', 'S_circular_insula_inf L',...
 'S_circular_insula_inf R', 'S_circular_insula_sup L', 'S_circular_insula_sup R',...
 'S_collat_transv_ant L', 'S_collat_transv_ant R', 'S_collat_transv_post L',...
 'S_collat_transv_post R', 'S_front_inf L', 'S_front_inf R', 'S_front_middle L',...
 'S_front_middle R', 'S_front_sup L', 'S_front_sup R', 'S_interm_prim-Jensen L',...
 'S_interm_prim-Jensen R', 'S_intrapariet_and_P_trans L', 'S_intrapariet_and_P_trans R',...
 'S_oc-temp_lat L', 'S_oc-temp_lat R', 'S_oc-temp_med_and_Lingual L',...
 'S_oc-temp_med_and_Lingual R', 'S_oc_middle_and_Lunatus L', 'S_oc_middle_and_Lunatus R',...
 'S_oc_sup_and_transversal L', 'S_oc_sup_and_transversal R', 'S_occipital_ant L',...
 'S_occipital_ant R', 'S_orbital-H_Shaped L', 'S_orbital-H_Shaped R', 'S_orbital_lateral L',...
 'S_orbital_lateral R', 'S_orbital_med-olfact L', 'S_orbital_med-olfact R',...
 'S_parieto_occipital L', 'S_parieto_occipital R', 'S_pericallosal L', 'S_pericallosal R',...
 'S_postcentral L', 'S_postcentral R', 'S_precentral-inf-part L', 'S_precentral-inf-part R',...
 'S_precentral-sup-part L', 'S_precentral-sup-part R', 'S_suborbital L', 'S_suborbital R',...
 'S_subparietal L', 'S_subparietal R', 'S_temporal_inf L', 'S_temporal_inf R',...
 'S_temporal_sup L', 'S_temporal_sup R', 'S_temporal_transverse L', 'S_temporal_transverse R'}}, ...
    'scoutfunc', in.scout_func, ...  % Mean
    'avgoutput', in.avg_output);

elseif in.scouts == 0

sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'clusters',  {}, ...
    'scoutfunc', in.scout_func, ...  % Mean
    'avgoutput', in.avg_output);


end

db.sub{i}.events{j}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},names.new_event_names{j}, in.output);
fclose(file_ID);
end

end

end
elseif iscell(in.events)
for i = 1:length(names.sub_name)

for j = 1:length(in.events)

tic

loc = find(strcmp(names.new_event_names, in.events{j}));

sFiles = db.sub{i}.events{loc}.(in.input);

if ischar(in.scouts)

% User scouts (from the scout_Destrieux_Auditory file)
sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'clusters',  {in.scouts, in.scout_list}, ...
    'scoutfunc', in.scout_func, ...  % Mean
    'avgoutput', in.avg_output);

elseif in.scouts == 1

% Process: Fourier transform (FFT)
sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'clusters',  {'Destrieux', {'G_Ins_lg_and_S_cent_ins L', 'G_Ins_lg_and_S_cent_ins R',...
 'G_and_S_cingul-Ant L', 'G_and_S_cingul-Ant R', 'G_and_S_cingul-Mid-Ant L',...
 'G_and_S_cingul-Mid-Ant R', 'G_and_S_cingul-Mid-Post L', 'G_and_S_cingul-Mid-Post R',...
 'G_and_S_frontomargin L', 'G_and_S_frontomargin R', 'G_and_S_occipital_inf L',...
 'G_and_S_occipital_inf R', 'G_and_S_paracentral L', 'G_and_S_paracentral R',...
 'G_and_S_subcentral L', 'G_and_S_subcentral R', 'G_and_S_transv_frontopol L',...
 'G_and_S_transv_frontopol R', 'G_cingul-Post-dorsal L', 'G_cingul-Post-dorsal R',...
 'G_cingul-Post-ventral L', 'G_cingul-Post-ventral R', 'G_cuneus L', 'G_cuneus R',...
 'G_front_inf-Opercular L', 'G_front_inf-Opercular R', 'G_front_inf-Orbital L',...
 'G_front_inf-Orbital R', 'G_front_inf-Triangul L', 'G_front_inf-Triangul R',...
 'G_front_middle L', 'G_front_middle R', 'G_front_sup L', 'G_front_sup R',...
 'G_insular_short L', 'G_insular_short R', 'G_oc-temp_lat-fusifor L',...
 'G_oc-temp_lat-fusifor R', 'G_oc-temp_med-Lingual L', 'G_oc-temp_med-Lingual R',...
 'G_oc-temp_med-Parahip L', 'G_oc-temp_med-Parahip R', 'G_occipital_middle L',...
 'G_occipital_middle R', 'G_occipital_sup L', 'G_occipital_sup R', 'G_orbital L',...
 'G_orbital R', 'G_pariet_inf-Angular L', 'G_pariet_inf-Angular R',...
 'G_pariet_inf-Supramar L', 'G_pariet_inf-Supramar R', 'G_parietal_sup L',...
 'G_parietal_sup R', 'G_postcentral L', 'G_postcentral R', 'G_precentral L',...
 'G_precentral R', 'G_precuneus L', 'G_precuneus R', 'G_rectus L', 'G_rectus R',...
 'G_subcallosal L', 'G_subcallosal R', 'G_temp_sup-G_T_transv L', 'G_temp_sup-G_T_transv R',...
 'G_temp_sup-Lateral L', 'G_temp_sup-Lateral R', 'G_temp_sup-Plan_polar L',...
 'G_temp_sup-Plan_polar R', 'G_temp_sup-Plan_tempo L', 'G_temp_sup-Plan_tempo R',...
 'G_temporal_inf L', 'G_temporal_inf R', 'G_temporal_middle L', 'G_temporal_middle R',...
 'Lat_Fis-ant-Horizont L', 'Lat_Fis-ant-Horizont R', 'Lat_Fis-ant-Vertical L',...
 'Lat_Fis-ant-Vertical R', 'Lat_Fis-post L', 'Lat_Fis-post R', 'Pole_occipital L',...
 'Pole_occipital R', 'Pole_temporal L', 'Pole_temporal R', 'S_calcarine L', 'S_calcarine R',...
 'S_central L', 'S_central R', 'S_cingul-Marginalis L', 'S_cingul-Marginalis R',...
 'S_circular_insula_ant L', 'S_circular_insula_ant R', 'S_circular_insula_inf L',...
 'S_circular_insula_inf R', 'S_circular_insula_sup L', 'S_circular_insula_sup R',...
 'S_collat_transv_ant L', 'S_collat_transv_ant R', 'S_collat_transv_post L',...
 'S_collat_transv_post R', 'S_front_inf L', 'S_front_inf R', 'S_front_middle L',...
 'S_front_middle R', 'S_front_sup L', 'S_front_sup R', 'S_interm_prim-Jensen L',...
 'S_interm_prim-Jensen R', 'S_intrapariet_and_P_trans L', 'S_intrapariet_and_P_trans R',...
 'S_oc-temp_lat L', 'S_oc-temp_lat R', 'S_oc-temp_med_and_Lingual L',...
 'S_oc-temp_med_and_Lingual R', 'S_oc_middle_and_Lunatus L', 'S_oc_middle_and_Lunatus R',...
 'S_oc_sup_and_transversal L', 'S_oc_sup_and_transversal R', 'S_occipital_ant L',...
 'S_occipital_ant R', 'S_orbital-H_Shaped L', 'S_orbital-H_Shaped R', 'S_orbital_lateral L',...
 'S_orbital_lateral R', 'S_orbital_med-olfact L', 'S_orbital_med-olfact R',...
 'S_parieto_occipital L', 'S_parieto_occipital R', 'S_pericallosal L', 'S_pericallosal R',...
 'S_postcentral L', 'S_postcentral R', 'S_precentral-inf-part L', 'S_precentral-inf-part R',...
 'S_precentral-sup-part L', 'S_precentral-sup-part R', 'S_suborbital L', 'S_suborbital R',...
 'S_subparietal L', 'S_subparietal R', 'S_temporal_inf L', 'S_temporal_inf R',...
 'S_temporal_sup L', 'S_temporal_sup R', 'S_temporal_transverse L', 'S_temporal_transverse R'}}, ...
    'scoutfunc', in.scout_func, ...  % Mean
    'avgoutput', in.avg_output);

elseif in.scouts == 0

sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'clusters',  {}, ...
    'scoutfunc', in.scout_func, ...  % Mean
    'avgoutput', in.avg_output);

end

db.sub{i}.events{loc}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.output);
fclose(file_ID);
end

end

end
elseif in.events == 0
for i = 1:length(names.sub_name)

tic

sFiles = db.sub{i}.(in.input);

if ischar(in.scouts)

% User scouts (from the scout_Destrieux_Auditory file)
sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'clusters',  {in.scouts, in.scout_list}, ...
    'scoutfunc', in.scout_func, ...  % Mean
    'avgoutput', in.avg_output);

elseif in.scouts == 1

% Process: Fourier transform (FFT)
sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'clusters',  {'Destrieux', {'G_Ins_lg_and_S_cent_ins L', 'G_Ins_lg_and_S_cent_ins R',...
 'G_and_S_cingul-Ant L', 'G_and_S_cingul-Ant R', 'G_and_S_cingul-Mid-Ant L',...
 'G_and_S_cingul-Mid-Ant R', 'G_and_S_cingul-Mid-Post L', 'G_and_S_cingul-Mid-Post R',...
 'G_and_S_frontomargin L', 'G_and_S_frontomargin R', 'G_and_S_occipital_inf L',...
 'G_and_S_occipital_inf R', 'G_and_S_paracentral L', 'G_and_S_paracentral R',...
 'G_and_S_subcentral L', 'G_and_S_subcentral R', 'G_and_S_transv_frontopol L',...
 'G_and_S_transv_frontopol R', 'G_cingul-Post-dorsal L', 'G_cingul-Post-dorsal R',...
 'G_cingul-Post-ventral L', 'G_cingul-Post-ventral R', 'G_cuneus L', 'G_cuneus R',...
 'G_front_inf-Opercular L', 'G_front_inf-Opercular R', 'G_front_inf-Orbital L',...
 'G_front_inf-Orbital R', 'G_front_inf-Triangul L', 'G_front_inf-Triangul R',...
 'G_front_middle L', 'G_front_middle R', 'G_front_sup L', 'G_front_sup R',...
 'G_insular_short L', 'G_insular_short R', 'G_oc-temp_lat-fusifor L',...
 'G_oc-temp_lat-fusifor R', 'G_oc-temp_med-Lingual L', 'G_oc-temp_med-Lingual R',...
 'G_oc-temp_med-Parahip L', 'G_oc-temp_med-Parahip R', 'G_occipital_middle L',...
 'G_occipital_middle R', 'G_occipital_sup L', 'G_occipital_sup R', 'G_orbital L',...
 'G_orbital R', 'G_pariet_inf-Angular L', 'G_pariet_inf-Angular R',...
 'G_pariet_inf-Supramar L', 'G_pariet_inf-Supramar R', 'G_parietal_sup L',...
 'G_parietal_sup R', 'G_postcentral L', 'G_postcentral R', 'G_precentral L',...
 'G_precentral R', 'G_precuneus L', 'G_precuneus R', 'G_rectus L', 'G_rectus R',...
 'G_subcallosal L', 'G_subcallosal R', 'G_temp_sup-G_T_transv L', 'G_temp_sup-G_T_transv R',...
 'G_temp_sup-Lateral L', 'G_temp_sup-Lateral R', 'G_temp_sup-Plan_polar L',...
 'G_temp_sup-Plan_polar R', 'G_temp_sup-Plan_tempo L', 'G_temp_sup-Plan_tempo R',...
 'G_temporal_inf L', 'G_temporal_inf R', 'G_temporal_middle L', 'G_temporal_middle R',...
 'Lat_Fis-ant-Horizont L', 'Lat_Fis-ant-Horizont R', 'Lat_Fis-ant-Vertical L',...
 'Lat_Fis-ant-Vertical R', 'Lat_Fis-post L', 'Lat_Fis-post R', 'Pole_occipital L',...
 'Pole_occipital R', 'Pole_temporal L', 'Pole_temporal R', 'S_calcarine L', 'S_calcarine R',...
 'S_central L', 'S_central R', 'S_cingul-Marginalis L', 'S_cingul-Marginalis R',...
 'S_circular_insula_ant L', 'S_circular_insula_ant R', 'S_circular_insula_inf L',...
 'S_circular_insula_inf R', 'S_circular_insula_sup L', 'S_circular_insula_sup R',...
 'S_collat_transv_ant L', 'S_collat_transv_ant R', 'S_collat_transv_post L',...
 'S_collat_transv_post R', 'S_front_inf L', 'S_front_inf R', 'S_front_middle L',...
 'S_front_middle R', 'S_front_sup L', 'S_front_sup R', 'S_interm_prim-Jensen L',...
 'S_interm_prim-Jensen R', 'S_intrapariet_and_P_trans L', 'S_intrapariet_and_P_trans R',...
 'S_oc-temp_lat L', 'S_oc-temp_lat R', 'S_oc-temp_med_and_Lingual L',...
 'S_oc-temp_med_and_Lingual R', 'S_oc_middle_and_Lunatus L', 'S_oc_middle_and_Lunatus R',...
 'S_oc_sup_and_transversal L', 'S_oc_sup_and_transversal R', 'S_occipital_ant L',...
 'S_occipital_ant R', 'S_orbital-H_Shaped L', 'S_orbital-H_Shaped R', 'S_orbital_lateral L',...
 'S_orbital_lateral R', 'S_orbital_med-olfact L', 'S_orbital_med-olfact R',...
 'S_parieto_occipital L', 'S_parieto_occipital R', 'S_pericallosal L', 'S_pericallosal R',...
 'S_postcentral L', 'S_postcentral R', 'S_precentral-inf-part L', 'S_precentral-inf-part R',...
 'S_precentral-sup-part L', 'S_precentral-sup-part R', 'S_suborbital L', 'S_suborbital R',...
 'S_subparietal L', 'S_subparietal R', 'S_temporal_inf L', 'S_temporal_inf R',...
 'S_temporal_sup L', 'S_temporal_sup R', 'S_temporal_transverse L', 'S_temporal_transverse R'}}, ...
    'scoutfunc', in.scout_func, ...  % Mean
    'avgoutput', in.avg_output);

elseif in.scouts == 0

sFiles = bst_process('CallProcess', 'process_fft', sFiles, [], ...
    'clusters',  {}, ...
    'scoutfunc', in.scout_func, ...  % Mean
    'avgoutput', in.avg_output);


end

db.sub{i}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.output);
fclose(file_ID);
end

end
end

end
bst_memory('UnloadAll', 'Forced');
save(strcat(names.ProtocolName,'_db.mat'),'db')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = extract_fcn(in, names, db)

if strcmp(in.events, 'all') 
for i = 1:length(names.sub_name)

for j = 1:length(names.new_event_names)

tic

sFiles = db.sub{i}.events{j}.(in.input);

% Extract values
% Process: Extract values: [0ms,1000ms] 0-100Hz abs
sFiles = bst_process('CallProcess', 'process_extract_values', sFiles, [], ...
    'timewindow', in.time_window, ...
    'freqrange',  in.freq_range, ...
    'rows',       in.rows, ...
    'isabs',      in.is_absolute, ...
    'avgtime',    in.avg_time, ...
    'avgrow',     in.avg_row, ...
    'avgfreq',    in.avg_freq, ...
    'matchrows',  in.match_rows, ...
    'dim',        in.dim, ...  % Concatenate time (dimension 2)
    'Comment',    '');

db.sub{i}.events{j}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Frequencies extracted\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Frequencies extracted\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j});
fclose(file_ID);
end

end

end

elseif iscell(in.events)
for i = 1:length(names.sub_name)

for j = 1:length(in.events)

tic

loc = find(strcmp(names.new_event_names, in.events{j}));

sFiles = db.sub{i}.events{loc}.(in.input);

% Extract values
% Process: Extract values: [0ms,1000ms] 0-100Hz abs
sFiles = bst_process('CallProcess', 'process_extract_values', sFiles, [], ...
    'timewindow', in.time_window, ...
    'freqrange',  in.freq_range, ...
    'rows',       in.rows, ...
    'isabs',      in.is_absolute, ...
    'avgtime',    in.avg_time, ...
    'avgrow',     in.avg_row, ...
    'avgfreq',    in.avg_freq, ...
    'matchrows',  in.match_rows, ...
    'dim',        in.dim, ...  % Concatenate time (dimension 2)
    'Comment',    '');

db.sub{i}.events{loc}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Frequencies extracted\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j});
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Frequencies extracted\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j});
fclose(file_ID);
end

end

end

elseif in.events == 0
for i = 1:length(names.sub_name)

tic

sFiles = db.sub{i}.(in.input);

% Extract values
% Process: Extract values: [0ms,1000ms] 0-100Hz abs
sFiles = bst_process('CallProcess', 'process_extract_values', sFiles, [], ...
    'timewindow', in.time_window, ...
    'freqrange',  in.freq_range, ...
    'rows',       in.rows, ...
    'isabs',      in.is_absolute, ...
    'avgtime',    in.avg_time, ...
    'avgrow',     in.avg_row, ...
    'avgfreq',    in.avg_freq, ...
    'matchrows',  in.match_rows, ...
    'dim',        in.dim, ...  % Concatenate time (dimension 2)
    'Comment',    '');

db.sub{i}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - Frequencies extracted for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.input);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - Frequencies extracted for %s\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.input);
fclose(file_ID);
end

end

end
bst_memory('UnloadAll', 'Forced');
save(strcat(names.ProtocolName,'_db.mat'),'db')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = difference_fcn(in,names, db)


for i = 1:length(names.sub_name)

tic

if any(strcmp(names.new_event_names,in.A_dir))
loc = find(strcmp(names.new_event_names, in.A_dir));
sFiles = db.sub{i}.events{loc}.(in.A_file); 
elseif ~any(strcmp(names.new_event_names,in.A_dir))
sFiles = db.sub{i}.(in.A_dir).(in.A_file); 
end

if any(strcmp(names.new_event_names,in.B_dir))
loc = find(strcmp(names.new_event_names, in.B_dir));
sFiles2 = db.sub{i}.events{loc}.(in.B_file); 
elseif ~any(strcmp(names.new_event_names,in.B_dir))
sFiles2 = db.sub{i}.(in.B_dir).(in.B_file); 
end

% Difference 
% Process: Difference: A-B
sFiles = bst_process('CallProcess', 'process_diff_ab', sFiles, sFiles2);

db.sub{i}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.output);
fclose(file_ID);
end

end
bst_memory('UnloadAll', 'Forced');
save(strcat(names.ProtocolName,'_db.mat'),'db')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [db] = group_tf_bands_fcn(in, names, db)

if strcmp(in.events, 'all') 
for i = 1:length(names.sub_name)

for j = 1:length(names.new_event_names)

sFiles = db.sub{i}.events{j}.(in.input);

% Process: Group in time or frequency bands
sFiles = bst_process('CallProcess', 'process_tf_bands', sFiles, [], ...
    'isfreqbands', in.is_freq_bands, ...
    'freqbands',   in.freq_bands, ...
    'istimebands', in.is_time_bands , ...
    'timebands',   in.time_bands, ...
    'overwrite',   in.overwrite);

db.sub{i}.events{j}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, names.new_event_names{j}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i},names.new_event_names{j}, in.output);
fclose(file_ID);
end

end

end
elseif iscell(in.events)
for i = 1:length(names.sub_name)

for j = 1:length(in.events)

tic

loc = find(strcmp(names.new_event_names, in.events{j}));

sFiles = db.sub{i}.events{loc}.(in.input);

% Process: Group in time or frequency bands
sFiles = bst_process('CallProcess', 'process_tf_bands', sFiles, [], ...
    'isfreqbands', in.is_freq_bands, ...
    'freqbands',   in.freq_bands, ...
    'istimebands', in.is_time_bands , ...
    'timebands',   in.time_bands, ...
    'overwrite',   in.overwrite);

db.sub{i}.events{loc}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1 && j == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.events{j}, in.output);
fclose(file_ID);
end

end

end
elseif in.events == 0
for i = 1:length(names.sub_name)

tic

sFiles = db.sub{i}.(in.input);

% Process: Group in time or frequency bands
sFiles = bst_process('CallProcess', 'process_tf_bands', sFiles, [], ...
    'isfreqbands', in.is_freq_bands, ...
    'freqbands',   in.freq_bands, ...
    'istimebands', in.is_time_bands , ...
    'timebands',   in.time_bands, ...
    'overwrite',   in.overwrite);

db.sub{i}.(in.output) = sFiles;

datafile = sFiles.FileName;

if strcmp(in.image_preview,'all')
save_image_tf(datafile,names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
elseif ismember(i,in.image_preview) 
save_image_tf(datafile, names.file_name{i}, names.image_folder_path, names.image_folder, in.image_file_tag)
end

t = toc;

if i == 1
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '\n%d min %.f sec - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.output);
fclose(file_ID);
else
log = strcat(names.ProtocolName,'_log','.txt');
file_ID = fopen(log, 'a');
fprintf(file_ID, '%d min %.f sec - %s - %s Processed\n',...
floor(t/60), rem(t,60), names.sub_name{i}, in.output);
fclose(file_ID);
end

end
end
bst_memory('UnloadAll', 'Forced');
save(strcat(names.ProtocolName,'_db.mat'),'db')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_image(data_file_name, subject_name, directory_to_image_folder, image_folder_name, image_file_tag)

bst_report('Snapshot', 'data', data_file_name,'',[]);
report_file_name = 'report_name';
ReportFile = bst_report('Save', data_file_name, report_file_name);
image_dir = fullfile(directory_to_image_folder,image_folder_name);

if ~exist(image_dir)
mkdir(image_dir)
end

image_name = strcat(image_file_tag,subject_name,'.png');
image_path = fullfile(image_dir, image_name);
bst_report('Export', ReportFile, image_path,'PNG');
bst_report('Close');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_image_tf(data_file_name, subject_name, directory_to_image_folder, image_folder_name, image_file_tag)
hFig = view_spectrum(data_file_name, 'Spectrum');
bst_report('Snapshot', hFig, [], '', [300 100 600 400]);
report_file_name = 'report_name';
ReportFile = bst_report('Save', data_file_name, report_file_name);
image_dir = fullfile(directory_to_image_folder,image_folder_name);

if ~exist(image_dir)
mkdir(image_dir)
end

image_name = strcat(image_file_tag,subject_name,'.png');
image_path = fullfile(image_dir, image_name);
bst_report('Export', ReportFile, image_path,'PNG');
bst_report('Close');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%