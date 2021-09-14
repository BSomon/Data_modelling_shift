function [visStimVals,visStimTimes,NbBipsMemo] = FuncSsvapNback_flicker(...
    obs_ind, VISUAL, AUDIO, SWITCH, trialLength, LSL, fullscreen, RatioX, RatioY,...
    trainSizes, fs, h, Ac, Am, fc, audFreq, Nbcycles, Famp, n_trialsCombi, visFreq,...
    mrk, timing, respKeyVisual, respKeyAudio, conditions,...
    Nback0, ProbaTarget)

%% Visual Auditory Steady State Response (VASSR) script
% presents 1 flickering square and 1 flickering sound
% the task is to attend to the visual, auditory or both modalities
% while doing a 0-back task on trains of points or auditory bips
% appearing in each modality

setpref('dsp','portaudioHostApi',3) ; % set ASIO driver for low latency

%% Setting up log file
dateLaunch=datestr(now, 30);
filename=sprintf('./results/Subj-%s-%s.txt', obs_ind, dateLaunch);
fid=fopen(filename,'w');
fprintf(fid,sprintf('obs_ind\ttrial_n\tcondition\ttask\tinfo\ttime\n'));

format long g
ListenChar(2);
rng('shuffle');
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up the screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting info
[windw, pixindeg, diffWandG, grey, xCenter, yCenter, ifi,...
    screenXpixels, screenYpixels] = screen_init(0, fullscreen);

% Rect stimulation settings
% For adapting the size of the rectangle to the size of the screen
% dimX=screenXpixels/RatioX; % dimension in x-axis proportional to screen size
% dimY=screenYpixels/RatioY; % dimension in y-axis proportional to screen size
% baseRect = [0 0 dimX dimY];

% For adapting the size of the rectangle to the visual degrees it should
% cover
dimInDeg=2;
dim = round(dimInDeg/2/pixindeg);
baseRect = [0 0 dim dim];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up the grid of square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xPos, yPos] = meshgrid(-3:1:3, -3:1:3);
% Calculate number of squares and reshape matrices of coordinates into a vector
[s1, s2] = size(xPos);
numSquares = s1 * s2;
xPos = reshape(xPos, 1, numSquares);
yPos = reshape(yPos, 1, numSquares);
xPos = xPos .* dim + xCenter;
yPos = yPos .* dim + yCenter;

% define colors of the checkers for black and white
% gridCols=[0 1];
% % gridCols = [0 0]; %for a flicker instead of checkers -- BS
% randBandW=repmat(repmat(gridCols,3,1),1,ceil(numSquares/2));
% randBandW = randBandW(:,1:numSquares);

% define the colors of the checkers for sine wave
gridCols=[0 1];
randBandW=repmat(repmat(gridCols,3,1),1,ceil(numSquares/2));
randBandW=randBandW(:,1:numSquares);
n_oddSquares=size(randBandW(:,1:2:end),2);
n_evenSquares=size(randBandW(:,2:2:end),2);

allRects = nan(4, numSquares);
for j = 1:numSquares
    allRects(:, j) = CenterRectOnPointd(baseRect, xPos(j), yPos(j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up the fixation cross and point stim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define central cross (displaying before trial start)
fixSizeDeg=.3; fixCrossDimPix = round(fixSizeDeg/2/pixindeg);
xCoords = [-fixCrossDimPix fixCrossDimPix-1 fixCrossDimPix-1 -fixCrossDimPix];
yCoords = [-fixCrossDimPix fixCrossDimPix -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords]; lineWidthPix = 3;

% Define visual points
sizecue=round(.7/pixindeg); halfsc=round(sizecue/2.);
[x, y] = meshgrid(-halfsc:halfsc, -halfsc:halfsc);
circleMat = sqrt((x .^ 2) + (y .^ 2)); outerVal=circleMat(halfsc,end);
circleMat(circleMat<=outerVal)=0; circleMat(circleMat>outerVal)=grey;
circleMat=cat(3,circleMat,circleMat,circleMat); a=circleMat(:,:,1);
amask=a==0; a(amask)=1; circleMat(:,:,1)=a;

alph=ones(size(circleMat,1),size(circleMat,1),1); alph(circleMat(:,:,1)==.5)=0;
circleMat=cat(3,circleMat,alph);

targetTexture=Screen('MakeTexture', windw, circleMat);

dimInDeg=1; dim = round(dimInDeg/2/pixindeg); % visual point dimensions
baseTarget = [0 0 dim dim];
targetLoc = CenterRectOnPointd(baseTarget, xCenter, yCenter); % point localisation (screen center)

%% Combinations task/condition, to display, randomly distributed

% Full number of trials
n_trials=n_trialsCombi*(length(conditions)+1);
indexCondition=[];
for i=1:size(conditions,2)
    indexCondition(end+1)=i;
end
indexCondition(end+1) = i;

indexCondition=repmat(indexCondition,1,n_trialsCombi);

indexRandom=Shuffle(1:n_trials);
% Random distribution of conditions
indexCondition=indexCondition(indexRandom);
iStim = 1;
if SWITCH
    %Distribution of first (audio or visual) for audiovisual
    StartAV = randperm(length(find(indexCondition==3)));
end

%% Initalize Psychtoolbox sound
if AUDIO
    InitializePsychSound;
    %         pahandle = PsychPortAudio('Open', 1, [], 0, fs, 1);
    pahandle = PsychPortAudio('Open', [], [], 0, fs, 1);
end

% Send first LSL trigger
if LSL; outlet = mrk.outlet; outlet.push_sample(5);end
% Range of periods (1/frequency modulation) between successive bips
RangePeriodBip=timing.IPI*audFreq;
% Range of periods (1/frequency modulation) between successive bips train
RangePeriodBipTrain=timing.ITrI*audFreq;
Total = {};
%% trial loop
for triali=1:n_trials
    logtowrite=[];
    
%     TimeUsed=0; % number of period used
%     NbBipsMemo=[]; % record number of bips for all trains
%     % Start period for each bips train computed, according to the number
%     % of trains allowed by trial length
%     PeriodStartBipsTrain={};
%     PeriodBip={};% Period between first bip and the start of train
%     while TimeUsed<((trialLength-timing.respWin)*audFreq)
%         % while number of periods used is less than trial periods number
%         weightsNbBip=ones(size(trainSizes,2),1);
%         % Weights to be applied on 0-back trainsize
%         weightsNbBip(trainSizes==Nback0)=ProbaTarget;
%         % Weights to be applied on other trainsizes
%         weightsNbBip(trainSizes~=Nback0)=(1-ProbaTarget)/(size(trainSizes,2)-1);
%         NbBips=datasample(trainSizes,1,'Weight',weightsNbBip); % number of bips per audio train
%         NbBipsMemo=[NbBipsMemo, NbBips];
%         
%         timepassedAlready=0;
%         if TimeUsed~=0
%             timepassedAlready=PeriodStartBipsTrain{end}+sum(PeriodBip{end});
%         end
%         % Start period of each audio bips train
%         PeriodStartBipsTrain{end+1}=timepassedAlready+datasample(RangePeriodBipTrain,1);
%         % Period between successive bip for one train (first bip at the
%         % beginning of the train)
%         PeriodBip{end+1}=[0 datasample(RangePeriodBip,NbBips-1)];
%         
%         % Compute number of periods used
%         TimeUsed=PeriodStartBipsTrain{end}+sum(PeriodBip{end});
%         if TimeUsed>((trialLength-timing.respWin)*audFreq) % Enough periods for a new train? No: remove last train
%             PeriodBip(end)=[]; PeriodStartBipsTrain(end)=[]; NbBipsMemo(end)=[];
%         end
%     end
%     
%     % Test audio target (match)
%     NbAudioTarget=0; % Nb of audio targets per trial
%     PeriodAudioTarget=[]; % Period of each audio targets
%     PeriodStartBipsTrain=cell2mat(PeriodStartBipsTrain);
%     for i=1:size(NbBipsMemo,2)
%         A=PeriodBip{i};
%         % if number of bips match the 0-back requirement
%         if NbBipsMemo(i)==Nback0
%             % Period target is defined by the period of the last bip in the target train
%             PeriodAudioTarget=[PeriodAudioTarget,(PeriodStartBipsTrain(i)+...
%                 sum(PeriodBip{i}))];
%             NbAudioTarget=NbAudioTarget+1;
%         end
%     end
%     
%     % Define the sound signal
%     k=h/Am; % amplification factor of the amplitude modulation
%     t=0:1/fs:trialLength; % samples number in the sound signal
%     carrier=Ac*cos(2*pi*fc.*t); % signal carrier
%     modul=Am*cos(2*pi*audFreq.*t); % modulated signal
%     signal=carrier+k.*modul.*carrier; % amplitude modulation of signal 'carrier' by signal 'modul'
%     
%     t1min=1/(2*audFreq); % first time where amplitude of 'sound' signal is minimal (in sec)
%     mini=t1min:1/audFreq:trialLength; %  times where amplitude of 'sound' signal are minimal (in sec)
%     
%     PeriodBipTrial=[]; % Define period for each bip during the trial
%     for i=1:size(PeriodStartBipsTrain,2)
%         A=PeriodBip{i};
%         for j=1:size(A,2)
%             PeriodBipTrial=[PeriodBipTrial,(PeriodStartBipsTrain(i)+sum(A(1:j)))];
%         end
%     end
%     
%     PeriodEndBipsTrain=[]; % Define period for each end bips train during the trial
%     for i=1:size(PeriodStartBipsTrain,2)
%         A=PeriodBip{i};
%         PeriodEndBipsTrain=[PeriodEndBipsTrain,(PeriodStartBipsTrain(i)+sum(A))];
%     end
%     
%     porte=ones(length(signal),1); % Define 'porte' signal to increase amplitude of 'sound' signal during bips
%     for i=1:length(PeriodBipTrial)
%         porte(round(mini(PeriodBipTrial(i))*fs):round(mini(PeriodBipTrial(i))*fs)+Nbcycles*round(fs/audFreq))=Famp;
%     end
%     target=signal.*porte';
    [target,PeriodBipTrial,mini,PeriodEndBipsTrain,PeriodStartBipsTrain,PeriodAudioTarget,NbBipsMemo] = ...
    Create_Sounds(trialLength,timing,audFreq,trainSizes,Nback0,ProbaTarget,RangePeriodBipTrain,RangePeriodBip,...
    h,Am,fs,Ac,fc,Nbcycles,Famp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Delay pre-instruction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tempTStart=GetSecs; firstPass=1;
    while GetSecs < (tempTStart+timing.preInstruction -(ifi/2))
        DrawFormattedText(windw, '', 'center', 'center', [1 1 1]);
        Screen('Flip', windw);
        if firstPass && LSL
            firstPass=0;
            outlet.push_sample(mrk.PreInstruction);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Delay instruction (combination task/condition) displaying
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tempTCue=GetSecs; firstPass=1;
    % instruction to display
    
    switch indexCondition(triali)
        case 1
            StimType = ' points';
        case 2
            StimType = ' bips';
        case 3
            if mod(StartAV(iStim), 2)%If start by audio
                StimType = ' bips \n Sart by AUDIO targets';
            else
                StimType = ' points \n Start by VISUAL targets';
            end
    end
    cue=strcat(conditions{indexCondition(triali)}, ...
        sprintf(' - You have to detect trains of %i ', Nback0), StimType);
    while GetSecs < (tempTCue+timing.timeInstruction-(ifi/2))
        DrawFormattedText(windw, cue, 'center', 'center', [1 1 1]);
        Screen('Flip', windw);
        if firstPass && LSL
            firstPass=0;
            switch indexCondition(triali)
                case 3
                    outlet.push_sample(mrk.Conditions(indexCondition(triali))+mod(StartAV(iStim),2));% marker 14 if audio strart and 13 if visuo start
                    iStim = iStim + 1;
                otherwise
                    outlet.push_sample(mrk.Conditions(indexCondition(triali)));
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Delay post instruction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tempTDelayPreStim=GetSecs; firstPass=1;
    while GetSecs < (tempTDelayPreStim+timing.postInstruction -(ifi/2))
        % Draw fixation cross
        Screen('DrawLines', windw, allCoords,lineWidthPix,[0 1 0],...
            [xCenter yCenter],0);
        Screen('Flip', windw);
        if firstPass && LSL; firstPass=0; outlet.push_sample(mrk.PostInstruction); end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Stim presentation and settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    firstPass=1; % test if first pass
    VisualTargetOn=0;   % boolean if visual target is on (1) or off (0)
    TimeLastVisualTarget=0;  % time of the last visual target
    TempKeyPressVisualTarget=GetSecs; % Time from the last visual key press
    TempKeyPressAudioTarget=GetSecs; % Time from the last audio key press
    indicAudioTarget=1; % increment for each audio target passed
    indicAudioBip=1;    % increment for each audio bip passed
    indicAudioStartTrain=1; % increment for each audio start train passed
    indicAudioEndTrain=1;  % increment for each audio end train passed
    AudioTargetOn=0;  % boolean if audio target is on (1) or off (1)
    VisualTrainMemory=[]; % store number of visual points for each train
    timeForNextVisTrain=1; % boolean whether there is time to present the next visual train
    visPeriodCount=0; % count the visual stimulation period passed
    checkColorSwitch=[]; % check if color of visual stimulus is changing
    timeUpdvisStim=[]; % times when the visual stimulation was computed/updated
    firstPassTrain=1; % boolean for start visual train mrk sending
    firstPassPoint=1; % boolean for visual point mrk sending
    
    timesWhenPointsWereDrawn=[]; % saves times when points were drawn to check the length of display
    
    % settings of each visual train
    % define probality for each available number of points
    weightsNbPoints=ones(size(trainSizes,2),1);
    weightsNbPoints(trainSizes==Nback0)=ProbaTarget;
    weightsNbPoints(trainSizes~=Nback0)=(1-ProbaTarget)/(size(trainSizes,2)-1);
    % define number of points for the next train
    NbPoints=datasample(trainSizes,1,'Weight',weightsNbPoints);
    NbPointsMemo = [];
    % number of available periods between successive points, according
    % to the visual stimulation frequency
    Apoint=ceil(timing.IPI*visFreq);
    % define, randomly, number of periods between successive points for the next train
    PeriodPoints=datasample(Apoint,NbPoints);
    PeriodPoints(1)=0; % first point is displayed at the start of the train (no periods delay)
    
    % number of available periods between successive point trains,
    % according to the visual stimulation frequency
    Atrain=ceil(timing.ITrI*visFreq);
    % define, randomly, start period of the next visual train
    PeriodVisualStartTrain = datasample(Atrain,1);
    
    PeriodEndLastVisualTrain=0; % define end period of the last visual train
    VisualTrainEnd=0; % boolean if visual train is finished (1) or not (0)
    countDisplayPoints=1; % count number of displayed points per train
    
    %Send a bip for switch trials whenever there is a N-back target
    if strcmp(conditions{indexCondition(triali)}, 'Audio & visual')
        if firstPass
            nAV = 0;
            idxVisuoTarget_changed  = zeros(length(NbBipsMemo)-1);
        end
        nAV = nAV+1;
        tic;
        while length(find(sum(idxVisuoTarget_changed,2)))<7 
            [NbPointsMemo, PeriodVisualStartTrain, PeriodEndLastVisualTrain,...
                PeriodPoints, NbVisualTarget, PeriodVisualTarget, miniV] = Create_Points(trialLength,...
                timing, visFreq, weightsNbPoints, trainSizes, Atrain, Apoint, Nback0,ifi);
            idxAudioTarget = find(NbBipsMemo == Nback0);
            idxVisuoTarget = find(NbPointsMemo == Nback0);
            AudioTargetTimeSwitch = mini(PeriodStartBipsTrain(idxAudioTarget));
            VisuoTargetTimeSwitch = miniV(PeriodVisualStartTrain(idxVisuoTarget));
            idxVisuoTarget_changed  = zeros(length(AudioTargetTimeSwitch)-1,length(VisuoTargetTimeSwitch));
            for i = 1:length(AudioTargetTimeSwitch)-1
                for j = 1:length(VisuoTargetTimeSwitch)
                    if VisuoTargetTimeSwitch(j)>=AudioTargetTimeSwitch(i)+2 && AudioTargetTimeSwitch(i+1)-2>=VisuoTargetTimeSwitch(j)
                        idxVisuoTarget_changed(i,j) = 1;
                    end
                end
            end
            if toc>5 %if no solution
                [target,PeriodBipTrial,mini,PeriodEndBipsTrain,PeriodStartBipsTrain,PeriodAudioTarget,NbBipsMemo] = ...
                Create_Sounds(trialLength,timing,audFreq,trainSizes,Nback0,ProbaTarget,RangePeriodBipTrain,RangePeriodBip,...
                h,Am,fs,Ac,fc,Nbcycles,Famp);
            end
        end
        
        %%
        %         idxAudioTarget = find(NbBipsMemo == Nback0);
        %         idxVisuoTarget = find(NbPointsMemo == Nback0);
        %         AudioTargetTimeSwitch = mini(PeriodStartBipsTrain(idxAudioTarget));
        %         VisuoTargetTimeSwitch = miniV(PeriodVisualStartTrain(idxVisuoTarget));
        % 		idxVisuoTarget_changed = zeros(1, length(idxVisuoTarget));
        % 		while idxVisuoTarget_changed ~= idxVisuoTarget
        % 			try
        % 				[idxAudioTarget, idxVisuoTarget_changed, AudioTargetTimeSwitch,...
        % 					VisuoTargetTimeSwitch, NbPointsMemo, PeriodVisualStartTrain,...
        % 					PeriodEndLastVisualTrain, PeriodPoints]=CheckAVidx(idxAudioTarget,...
        % 					idxVisuoTarget, AudioTargetTimeSwitch, VisuoTargetTimeSwitch,...
        % 					NbPointsMemo, NbBipsMemo, Nback0, mini, miniV, PeriodStartBipsTrain, ...
        % 					PeriodVisualStartTrain, PeriodEndLastVisualTrain, PeriodPoints);
        
        %             catch
        %         [NbPointsMemo, PeriodVisualStartTrain, PeriodEndLastVisualTrain,...
        %             PeriodPoints, NbVisualTarget, PeriodVisualTarget, miniV] = Create_Points(trialLength,...
        %             timing, visFreq, weightsNbPoints, trainSizes, Atrain, Apoint, Nback0);
        %         idxAudioTarget = find(NbBipsMemo == Nback0);
        %         idxVisuoTarget = find(NbPointsMemo == Nback0);
        %         AudioTargetTimeSwitch = mini(PeriodStartBipsTrain(idxAudioTarget));
        %         VisuoTargetTimeSwitch = miniV(PeriodVisualStartTrain(idxVisuoTarget));
        %         idxVisuoTarget_changed = zeros(1, length(idxVisuoTarget));
        % 			end
        % 		end
        
        tmpPeriodEndLastVisualTrain = [0,PeriodEndLastVisualTrain];
        tmpPeriodVisualStartTrain = [PeriodVisualStartTrain, 180*visFreq];
        tmpPeriodPoints = [PeriodPoints, [0 0 0]];%Adding a last value for programming issues
        tmpNbPointsMemo = [NbPointsMemo, 0];%Adding a last value for programing issues
        
    end
    if AUDIO; PsychPortAudio('FillBuffer', pahandle, target); end
    
    
    if LSL; outlet.push_sample(mrk.StimOn);end
    
    StartStim=GetSecs; % Start stimulation time
    tempTFlick2=StartStim; % time for visual stimulation
    %% Start stimulation
    while GetSecs < (StartStim+trialLength-(ifi/2))
        FlushEvents('keyDown');
        
        % Visual stimulation
        if VISUAL
            % store at which moment the new color was computed (to
            % verify how often it is updated)
            timeUpdvisStim=[timeUpdvisStim GetSecs-StartStim];
            % Test with a sin wave and check board
            % compute the time at which each of the odd and even
            % squares are to have an opposite phase between them
            xs=[timeUpdvisStim(end).*(pi*2*visFreq), timeUpdvisStim(end).*(pi*2*visFreq)+pi];
            %             xs=[timeUpdvisStim(end).*(pi*2*visFreq), timeUpdvisStim(end).*(pi*2*visFreq)];
            
            % compute the color for odd and even squares
            visStimCol=(repmat(sin(xs),3,1)+1)/2;
            % store color for all odd squares (one half of the squares)
            randBandW(:,1:2:end)=repmat(visStimCol(:,1),1,n_oddSquares);
            % store color for all even squares (one half of the squares)
            randBandW(:,2:2:end)=repmat(visStimCol(:,2),1,n_evenSquares);
            
            Screen('FillRect', windw, randBandW, allRects);
            checkColorSwitch=[checkColorSwitch visStimCol(1,1)];
            if size(checkColorSwitch,2)>1
                if checkColorSwitch(end-1) > .5 && checkColorSwitch(end) < .5
                    visPeriodCount=visPeriodCount+0.5;
                elseif checkColorSwitch(end-1) < .5 && checkColorSwitch(end) > .5
                    visPeriodCount=visPeriodCount+0.5;
                end
            end
            % Draw fixation cross
            Screen('DrawLines', windw, allCoords,lineWidthPix,[0 1 0],[xCenter yCenter],0);
        end
        
        % Audio stimulation
        if firstPass
            if AUDIO
                PsychPortAudio('Start', pahandle, 1, 0, 1);
                % Define audio targets time in sec
                TimeAudioTarget=GetSecs+mini(PeriodAudioTarget);
                % Define audio bips time in sec
                TimeBipTrial=GetSecs+mini(PeriodBipTrial);
                % Define time, in sec, of audio start trains
                TimeStartBipsTrain=GetSecs+mini(PeriodStartBipsTrain);
                % Define time, in sec, of audio end trains
                TimeEndBipsTrain=GetSecs+mini(PeriodEndBipsTrain);
            end
            if strcmp(conditions{indexCondition(triali)}, 'Audio & visual')
                PeriodEndLastVisualTrain = tmpPeriodEndLastVisualTrain(1);
                PeriodVisualStartTrain = tmpPeriodVisualStartTrain(1);
                PeriodPoints = tmpPeriodPoints{1};
                NbPoints = tmpNbPointsMemo(1);
            end
            firstPass=0;
        end
        
        % Visual points display
        if timeForNextVisTrain %&& ~strcmp(conditions{indexCondition(triali)}, 'Audio & visual')% check if visual train is okay
            
            if ~VisualTrainEnd
                % compute the first period of the current visual train
                startPeriodPointOn=PeriodEndLastVisualTrain+PeriodVisualStartTrain+...
                    sum(PeriodPoints(1:countDisplayPoints));
                
                % compute the last period of the current visual train
                endPeriodPointOn=PeriodEndLastVisualTrain+PeriodVisualStartTrain+...
                    sum(PeriodPoints(1:countDisplayPoints))+timing.cyclePointOn;
                
                if visPeriodCount >= startPeriodPointOn && visPeriodCount <= endPeriodPointOn
                    % display visual points between periods defined previously
                    timesWhenPointsWereDrawn=[timesWhenPointsWereDrawn GetSecs-StartStim];
                    Screen('DrawTextures', windw, targetTexture, [], targetLoc);
                    
                    if firstPassTrain % if this is the first point of the train
                        disp 'Start visual Train'; firstPassTrain=0;
                        if countDisplayPoints==1 && LSL
                            outlet.push_sample(mrk.StartVisualTrain...
                                (trainSizes==NbPoints));
                        end
                    end
                    
                    
                    if firstPassPoint % if a point was just displayed
                        disp 'Visual point'; firstPassPoint=0;
                        % Write in log file
                        logtowrite=[logtowrite sprintf('%i\t%i\t%s\t%.5f\n',obs_ind,triali,'visual',GetSecs()-StartStim)];
                        if LSL; outlet.push_sample(mrk.VisualPoints); end
                    end
                    
                    % if this is the last point of the current train
                    if countDisplayPoints==NbPoints && visPeriodCount >= endPeriodPointOn
                        disp 'End visual Train';
                        VisualTrainEnd=1; countDisplayPoints=0;
                        if LSL; outlet.push_sample(mrk.EndVisualTrain(trainSizes==NbPoints)); end
                    end
                    
                    % increment when a point was displayed
                    if visPeriodCount == endPeriodPointOn
                        countDisplayPoints=countDisplayPoints+1;
                        firstPassPoint=1;
                    end
                end
            end
            
        end
        
        if VisualTrainEnd % test if visual train has ended
            firstPassTrain=1;
            % store number of visual points for each train
            VisualTrainMemory=[VisualTrainMemory, NbPoints];
            
            % test if there is a visual target, according to the last visual train
            if VisualTrainMemory(end)==Nback0
                disp 'Visual target';
                VisualTargetOn=1;
                TimeLastVisualTarget=GetSecs;
                if LSL; outlet.push_sample(mrk.VisualTarget); end
            end
            
            if strcmp(conditions{indexCondition(triali)}, 'Audio & visual')
                tmpPeriodEndLastVisualTrain = tmpPeriodEndLastVisualTrain(2:end);
                tmpPeriodVisualStartTrain = tmpPeriodVisualStartTrain(2:end);
                try
                    tmpPeriodPoints = tmpPeriodPoints(2:end);
                    tmpNbPointsMemo = tmpNbPointsMemo(2:end);
                end
                PeriodEndLastVisualTrain = tmpPeriodEndLastVisualTrain(1);
                PeriodVisualStartTrain = tmpPeriodVisualStartTrain(1)-PeriodEndLastVisualTrain;
                PeriodPoints = tmpPeriodPoints{1};
                NbPoints = tmpNbPointsMemo(1);
            else
                % Reinitialization of points settings (number, display time)
                % for the next visual train
                weightsNbPoints=ones(size(trainSizes,2),1); % define probality for each available number of points
                weightsNbPoints(trainSizes==Nback0)=ProbaTarget;
                weightsNbPoints(trainSizes~=Nback0)=(1-ProbaTarget)/(size(trainSizes,2)-1);
                
                NbPoints=datasample(trainSizes,1,'Weight',weightsNbPoints);
                
                % number of available periods between successive points,
                % according to the visual stimulation frequency
                Apoint=ceil(timing.IPI*visFreq);
                % define, randomly, number of periods between successive points for the next train
                PeriodPoints=datasample(Apoint,NbPoints);
                % first point is displayed at the start of the train (no periods delay)
                PeriodPoints(1)=0;
                
                % number of available periods between successive points
                % trains, according to the visual stimulation frequency
                Atrain=ceil(timing.ITrI*visFreq);
                % define, randomly, start period of the next visual train
                PeriodVisualStartTrain=datasample(Atrain,1);
                
                % define end period of the last visual train
                PeriodEndLastVisualTrain=visPeriodCount;
            end
            VisualTrainEnd=0; % boolean if visual train is end (1) or not (0)
            
            if (trialLength*visFreq) < (PeriodEndLastVisualTrain+PeriodVisualStartTrain+sum(PeriodPoints) + timing.respWin * visFreq)
                % test if there is time to display the next visual train
                timeForNextVisTrain=0;
            end
        end
        
        Screen('Flip', windw);
        
        
        
        % Test for start audio train
        if AUDIO && indicAudioStartTrain<=length(TimeStartBipsTrain) &&...
                round(TimeStartBipsTrain(indicAudioStartTrain),1)==...
                round(GetSecs,1)
            %             if AUDIO && indicAudioStartTrain<=length(PeriodStartBipsTrain) &&...
            %                     round(PeriodStartBipsTrain(indicAudioStartTrain),1)==round(GetSecs,1)
            disp('Start audio train')
            if LSL; outlet.push_sample(mrk.StartAudioTrain...
                    (trainSizes==NbBipsMemo(indicAudioStartTrain))); end
            indicAudioStartTrain=indicAudioStartTrain+1;
        end
        
        % Test for audio bip
        if AUDIO && indicAudioBip<=length(TimeBipTrial) &&...
                round(TimeBipTrial(indicAudioBip),1)==round(GetSecs,1)
            disp('Audio Bip')
            
            % Write in log file
            logtowrite=[logtowrite sprintf('%i\t%i\t%s\t%.5f\n',obs_ind,triali,'audio',GetSecs()-StartStim)];
            
            indicAudioBip=indicAudioBip+1;
            if LSL;outlet.push_sample(mrk.AudioBip); end
        end
        
        % Test for audio target
        if AUDIO && indicAudioTarget<=length(TimeAudioTarget) &&...
                round(TimeAudioTarget(indicAudioTarget),1)==round(GetSecs,1)
            disp('Audio target')
            indicAudioTarget=indicAudioTarget+1;
            AudioTargetOn=1;
            if LSL;outlet.push_sample(mrk.AudioTarget); end
        end
        
        
        % Test for end audio train
        if AUDIO && indicAudioEndTrain<=length(TimeEndBipsTrain) &&...
                round(TimeEndBipsTrain(indicAudioEndTrain),1)==...
                round(GetSecs,1)
            disp('End audio train')
            if LSL; outlet.push_sample(mrk.EndAudioTrain(trainSizes...
                    ==NbBipsMemo(indicAudioEndTrain)));  end
            indicAudioEndTrain=indicAudioEndTrain+1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        % Catch response
        [keyIsDown, keyTime, keyCode] = KbCheck;
        keyTemp=KbName(keyCode);
        RefractoryKeyTime = toc;
        
        if keyIsDown && sum(strcmp(keyTemp, [respKeyVisual, respKeyAudio]))...
                && RefractoryKeyTime>1.5
%             event=sprintf('resp_%s',keyTemp);
%             logtowrite=[logtowrite sprintf('%i\t%i\t%s\t%.5f\n',obs_ind,triali,event,GetSecs()-StartStim)];
            
            if LSL && strcmp(keyTemp,respKeyVisual)
                outlet.push_sample(mrk.keyRespVisual);
            elseif LSL && strcmp(keyTemp,respKeyAudio)
                outlet.push_sample(mrk.keyRespAudio);
            end
            tic;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%
        
    end % end stimulation
    Total = [Total, {VisualTrainMemory}];
    if LSL; outlet.push_sample(mrk.TrialEnd); end
    
    NextTrial = 0;
    
    DrawFormattedText(windw, 'Press space bar for next trial', 'center', 'center', [1 1 1]);
    Screen('Flip', windw);
    while ~NextTrial % wait for the next trial while 'space bar' is not pressed
        [keyIsDown, keyTime, keyCode] = KbCheck;
        if keyIsDown==1; NextTrial=1; end
    end
    
    visStimVals{triali}=checkColorSwitch;
    visStimTimes{triali}=timeUpdvisStim;
    
    % finally write the trial events in a log file
    fprintf(fid, logtowrite);
    save([num2str(obs_ind),'.mat'], 'Total')
end % end trial

sca;
fclose('all'); % close log file
ListenChar(0);
end


function [NbPointsMemo, PeriodVisualStartTrain, PeriodEndLastVisualTrain,...
    PeriodPoints, NbVisualTarget, PeriodVisualTarget, miniV] = Create_Points(trialLength,...
    timing, visFreq, weightsNbPoints, trainSizes, Atrain, Apoint, Nback0,ifi)

%% ADD POINTS
TimeUsed=0; % number of period used
NbPointsMemo=[]; % record number of bips for all trains
% Start period for each bips train computed, according to the number
% of trains allowed by trial length
PeriodVisualStartTrain={};
PeriodPoints={};% Period between first bip and the start of train
while TimeUsed<((trialLength-ifi/2)*visFreq)
    NbPoints=datasample(trainSizes,1,'Weight',weightsNbPoints); % number of bips per audio train
    NbPointsMemo=[NbPointsMemo, NbPoints];
    timepassedAlready=0;
    if TimeUsed~=0
        timepassedAlready=PeriodVisualStartTrain{end}+sum(PeriodPoints{end});
    end
    % Start period of each visual point train
    PeriodVisualStartTrain{end+1}=timepassedAlready+datasample(Atrain,1);
    % Period between successive points for one train (first point at the
    % beginning of the train)
    PeriodPoints{end+1}=[0 datasample(Apoint,NbPoints-1)];
    
    % Compute number of periods used
    TimeUsed=PeriodVisualStartTrain{end}+sum(PeriodPoints{end});
    if TimeUsed>((trialLength-ifi/2)*visFreq) % Enough periods for a new train? No: remove last train
        PeriodPoints(end)=[]; PeriodVisualStartTrain(end)=[]; NbPointsMemo(end)=[];
    end
end

% Test audio target (match)
NbVisualTarget=0; % Nb of visual targets per trial
PeriodVisualTarget=[]; % Period of each audio targets
PeriodVisualStartTrain=cell2mat(PeriodVisualStartTrain);
for i=1:size(NbPointsMemo,2)
    A=PeriodPoints{i};
    % if number of bips match the 0-back requirement
    if NbPointsMemo(i)==Nback0
        % Period target is defined by the period of the last bip in the target train
        PeriodVisualTarget=[PeriodVisualTarget,(PeriodVisualStartTrain(i)+...
            sum(PeriodPoints{i}))];
        NbVisualTarget=NbVisualTarget+1;
    end
end

% 		t1minV=1/(2*visFreq); % first time where amplitude of 'sound' signal is minimal (in sec)
miniV=(1/(2*visFreq)):(1/visFreq):trialLength; %  times where amplitude of 'sound' signal are minimal (in sec)

PeriodEndLastVisualTrain=[]; % Define period for each end bips train during the trial
for i=1:size(PeriodVisualStartTrain,2)
    A=PeriodPoints{i};
    PeriodEndLastVisualTrain=[PeriodEndLastVisualTrain,(PeriodVisualStartTrain(i)+sum(A))];
end

disp('Visual Points created!')

end

function [target,PeriodBipTrial,mini,PeriodEndBipsTrain,PeriodStartBipsTrain,PeriodAudioTarget,NbBipsMemo] = ...
    Create_Sounds(trialLength,timing,audFreq,trainSizes,Nback0,ProbaTarget,RangePeriodBipTrain,RangePeriodBip,...
    h,Am,fs,Ac,fc,Nbcycles,Famp)

    TimeUsed=0; % number of period used
    NbBipsMemo=[]; % record number of bips for all trains
    % Start period for each bips train computed, according to the number
    % of trains allowed by trial length
    PeriodStartBipsTrain={};
    PeriodBip={};% Period between first bip and the start of train
    while TimeUsed<((trialLength-timing.respWin)*audFreq)
        % while number of periods used is less than trial periods number
        weightsNbBip=ones(size(trainSizes,2),1);
        % Weights to be applied on 0-back trainsize
        weightsNbBip(trainSizes==Nback0)=ProbaTarget;
        % Weights to be applied on other trainsizes
        weightsNbBip(trainSizes~=Nback0)=(1-ProbaTarget)/(size(trainSizes,2)-1);
        NbBips=datasample(trainSizes,1,'Weight',weightsNbBip); % number of bips per audio train
        NbBipsMemo=[NbBipsMemo, NbBips];
        
        timepassedAlready=0;
        if TimeUsed~=0
            timepassedAlready=PeriodStartBipsTrain{end}+sum(PeriodBip{end});
        end
        % Start period of each audio bips train
        PeriodStartBipsTrain{end+1}=timepassedAlready+datasample(RangePeriodBipTrain,1);
        % Period between successive bip for one train (first bip at the
        % beginning of the train)
        PeriodBip{end+1}=[0 datasample(RangePeriodBip,NbBips-1)];
        
        % Compute number of periods used
        TimeUsed=PeriodStartBipsTrain{end}+sum(PeriodBip{end});
        if TimeUsed>((trialLength-timing.respWin)*audFreq) % Enough periods for a new train? No: remove last train
            PeriodBip(end)=[]; PeriodStartBipsTrain(end)=[]; NbBipsMemo(end)=[];
        end
    end
    
    % Test audio target (match)
    NbAudioTarget=0; % Nb of audio targets per trial
    PeriodAudioTarget=[]; % Period of each audio targets
    PeriodStartBipsTrain=cell2mat(PeriodStartBipsTrain);
    for i=1:size(NbBipsMemo,2)
        A=PeriodBip{i};
        % if number of bips match the 0-back requirement
        if NbBipsMemo(i)==Nback0
            % Period target is defined by the period of the last bip in the target train
            PeriodAudioTarget=[PeriodAudioTarget,(PeriodStartBipsTrain(i)+...
                sum(PeriodBip{i}))];
            NbAudioTarget=NbAudioTarget+1;
        end
    end
    
    % Define the sound signal
    k=h/Am; % amplification factor of the amplitude modulation
    t=0:1/fs:trialLength; % samples number in the sound signal
    carrier=Ac*cos(2*pi*fc.*t); % signal carrier
    modul=Am*cos(2*pi*audFreq.*t); % modulated signal
    signal=carrier+k.*modul.*carrier; % amplitude modulation of signal 'carrier' by signal 'modul'
    
    t1min=1/(2*audFreq); % first time where amplitude of 'sound' signal is minimal (in sec)
    mini=t1min:1/audFreq:trialLength; %  times where amplitude of 'sound' signal are minimal (in sec)
    
    PeriodBipTrial=[]; % Define period for each bip during the trial
    for i=1:size(PeriodStartBipsTrain,2)
        A=PeriodBip{i};
        for j=1:size(A,2)
            PeriodBipTrial=[PeriodBipTrial,(PeriodStartBipsTrain(i)+sum(A(1:j)))];
        end
    end
    
    PeriodEndBipsTrain=[]; % Define period for each end bips train during the trial
    for i=1:size(PeriodStartBipsTrain,2)
        A=PeriodBip{i};
        PeriodEndBipsTrain=[PeriodEndBipsTrain,(PeriodStartBipsTrain(i)+sum(A))];
    end
    
    porte=ones(length(signal),1); % Define 'porte' signal to increase amplitude of 'sound' signal during bips
    for i=1:length(PeriodBipTrial)
        porte(round(mini(PeriodBipTrial(i))*fs):round(mini(PeriodBipTrial(i))*fs)+Nbcycles*round(fs/audFreq))=Famp;
    end
    target=signal.*porte';
end
