    %% Visual Auditory Steady State Response (VASSR) script and Nback task
% Initially written by Mehdi Senoussi
% Readapted by Bertille Somon, 2020

close all; clear all; clc

% boolean to use/discard visual/auditory/Switch stimuli and points for
% visual N-back task
VISUAL=1; AUDIO=1; SWITCH=1;

% boolean to send triggers to the biosemi system and to display fullscreen
fullscreen=1;  LSL=1;
% Screen('Preference', 'SkipSyncTests', 1);
% Screen('Preference','VisualDebugLevel', 0);
Screen('Preference','TextRenderer', 0); %Remove PTB bug - BS
addpath(genpath('G:\Expé\Toolboxes\LiveAmp'))
% Instantiate the LSL library -BS
if LSL
    disp('Loading library...');
    lib = lsl_loadlib();
    % make a new stream outlet
    disp('Creating a new streaminfo...');
    info = lsl_streaminfo(lib,'Triggers','Markers',1,0,'cf_int32',...
		'ssvap_markers');
    disp('Opening an outlet...');
    outlet = lsl_outlet(info);
    disp('Do not forget to launch LabRecorder... Press any key to continue.')
    pause;
	mrk.outlet = outlet;
end


% observer number
obs_ind = input('Participant number: ', 's');

% Trial parameters
trialLength=180; % Time of each trial (sec)
n_trialsCombi = 5; % Number of trials per combination condition/task (ex. visual condition with 1-back task)
conditions={'Visual' 'Audio' 'Audio & visual'}; % Instructions/conditions to display before each trial

% timings of the parts of each trial
timing.preInstruction=1; % time before instruction display
timing.timeInstruction=5; % time for instruction display
timing.postInstruction=1; % time after instruction display

% Visual stimulus (rectangle)
visFreq=48; % flickering frequencies of the checkers
RatioX=1; % X dimension = Screen X size/Ratio
RatioY=1; % Y dimension = Screen Y size/Ratio

% Auditory stimulus
fs=44100; % sample frequency of sound signal
h=1; % modulation amplitude index (0<h<1)
Ac=1;  % Amplitude of carrier signal
Am=1;  % Amplitude modulation of carrier signal
fc=500; % Frequency of carrier signal
audFreq=40; % Frequency modulation of carrier signal amplitude
Nbcycles=2; % bip length in number of periods in terms of frequency modulation
Famp=2; % Increased factor of signal amplitude during bips 

% N-back task settings
trainSizes=[2 3 4]; % random number of points/bips in each visual/audio train 
Nback0=2; % number to detect for 0-back task
% Inter Point Interval: possible delays between successive points/bips (sec)
timing.IPI=0.7:0.05:.9;
% Inter Train Interval: possible delays between successive points/bips train (sec)
timing.ITrI=2:0.1:3;

% number of periods (of stimulus carry freq) while points are displayed
timing.cyclePointOn=2;%2
ProbaTarget=0.3; % probability of target per trial

% response parameters
respKeyVisual={'p'}; % visual key to press when visual target 
respKeyAudio={'a'}; % audio key to press when audio target
timing.respWin=2; % Time delay to respond to visual/audio targets (sec)
timing.respInt=1; % Refractory time between two successive responses (sec)


% Triggers index
mrk.PreInstruction=10; % mrk pre-instruction
vConditions = [11 12 13]; %14 for audio start
mrk.Conditions=vConditions(1:length(conditions)); % conditions mrk, same length as 'conditions'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mrk.PostInstruction=15; % mrk post-instruction
mrk.StimOn=20; % mrk stimulation onset
mrk.TrialEnd=40; % mrk trial end

mrk.keyRespAudio=50; % mrk keyboard auditory target response
mrk.keyRespVisual=52; % mrk keyboard visual target response
mrk.AudioTarget=70; % mrk audio target
mrk.AudioTargetTrue=71; % audio target correct response 
mrk.AudioTargetFalse=72; % audio target incorrect response 
mrk.VisualTarget=60; % mrk visual target
mrk.VisualTargetTrue=61; % Visual target correct response
mrk.VisualTargetFalse=62; % Visual target incorrect response

mrk.VisualPoints=101; % mrk visual points
% mrk start visual train (one mrk for each train length)
mrk.StartVisualTrain=100+trainSizes;
% mrk end visual train (one mrk for each train length)
mrk.EndVisualTrain=106+trainSizes;

mrk.AudioBip = [112]; % audio bip
% mrk start audio train (one mrk for each train length)
mrk.StartAudioTrain = 112+trainSizes;
% mrk end audio train (one mrk for each train length)
mrk.EndAudioTrain = 117+trainSizes;


%% Visual Auditory Steady State Response (VASSR) and N-back task function calling

[visStimVals,visStimTimes] = FuncSsvapNback_flicker(obs_ind, ...
     VISUAL, AUDIO, SWITCH, trialLength, LSL, fullscreen, RatioX, RatioY,...
    trainSizes, fs, h, Ac, Am, fc, audFreq, Nbcycles, Famp, n_trialsCombi, visFreq,...
    mrk, timing, respKeyVisual, respKeyAudio, conditions,...
    Nback0, ProbaTarget);
Screen('CloseAll')