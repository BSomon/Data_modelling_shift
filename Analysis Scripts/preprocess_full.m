%% EEG data preprocessing for VASSR analyses with trials extraction and computation of perf/RT etc for each condition (A, V, VA, AV)
SUBJ_tot = {'PCB0805' 'PEJ0511' 'PMFH0507' 'PND1607' 'PSV2801' 'PJFS0801' 'PMF0711' ...
    'PSL2804' 'PEF0102' 'PJE0610' 'PLD2501' 'PZG0801' 'PZG2305' 'PBJ0703' 'PAI0912'...
    'PAS0505' 'PJJTT0601' 'POL3007' 'PVB2304' 'PMR1501' 'PCH0107' 'PNP1504' 'PTH1804'...
    'PCD0306'};
%Subject to analyse
SUBJ = {'PCD0306'};
SUBJ_rej = {'PRS1907'};
%Select the processing step you want to create 
Preprocess = 0;
ICA = 1;
for s = 1:length(SUBJ)
    if Preprocess
        EEG = pop_loadset('filename',['sub-',SUBJ{s},'_ses-S001_task-Default_run-001_eeg.set'],'filepath',['G:/Expé/ANITI/Data/sub-',SUBJ{s},'/ses-S001/eeg/']);
        EEG = pop_headplot(EEG, 1); %Creating the spline file required further
        EEG = pop_eegfiltnew(EEG,1,[],[]);
        [zapdata,EEG.artifacts]=nt_zapline(EEG.data',50/500,1);
        EEG.data = zapdata';
        [~,filename] = fileparts(EEG.filename);
        EEG = pop_saveset(EEG,'filename', strcat([filename,'_HP']), 'filepath', EEG.filepath);
        EEG = pop_reref( EEG, [10 21]);
        [~,filename] = fileparts(EEG.filename);
        EEG = pop_saveset(EEG,'filename', strcat([filename,'_reref']), 'filepath', EEG.filepath);
        
        %ASR
        originalEEG = EEG;
        EEG=clean_asr(EEG,5,[],[],[],'off');
        %    Interpolate all the removed channels
        EEG = pop_interp(EEG, originalEEG.chanlocs, 'spherical');
        [~,filename] = fileparts(EEG.filename);
        EEG = pop_saveset(EEG,'filename', strcat([filename,'_ASR']), 'filepath', EEG.filepath);
        
        %ICA
        EEG = pop_runica(EEG, 'icatype','picard', 'maxiter', 500);
        EEG = eeg_checkset(EEG);
        [~,filename] = fileparts(EEG.filename);
        EEG = pop_saveset(EEG,'filename', strcat([filename,'_ICA']), 'filepath', EEG.filepath);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
    elseif ICA
        %Requires to remove bad ICs beforehand
        EEG = pop_loadset('filename',['sub-',SUBJ{s},'_ses-S001_task-Default_run-001_eeg_HP_reref_ASR_ICAp.set'],'filepath',['G:/Expé/ANITI/Data/sub-',SUBJ{s},'/ses-S001/eeg/']);
        EEG = pop_epoch( EEG, {  '11'  '12'  '13'  '14'  }, [0  180], 'newname', 'AllTrials', 'epochinfo', 'yes'); %5 seconds at the begining for instructions
        
        %     Epoching
        Table_stim(SUBJ(s))
        load('Trials_All.mat')
        idx_trials = find(contains(Trials_All(:,1), SUBJ{s}));
        [~,filename] = fileparts(EEG.filename);
        EEG1 = pop_selectevent( EEG, 'event',Trials_All{idx_trials, 2}' ,'deleteevents','on');
        EEG1 = pop_epoch(EEG1, {'60'}, [-2 2], 'newname','VisualTarget', 'epochinfo', 'yes');
        pop_saveset(EEG1, 'filename', strcat([filename,'_VisualTargetResp.set']), 'filepath', EEG.filepath);
        EEG2 = pop_selectevent( EEG, 'event',Trials_All{idx_trials, 3}' ,'deleteevents','on');
        EEG2 = pop_epoch(EEG2, {'70'}, [-2 2], 'newname','AudioTarget', 'epochinfo', 'yes');
        pop_saveset(EEG2, 'filename', strcat([filename,'_AudioTargetResp.set']), 'filepath', EEG.filepath);
        EEG3 = pop_selectevent( EEG, 'event',Trials_All{idx_trials, 4}' ,'deleteevents','on');
        EEG3 = pop_epoch(EEG3, {'60'}, [-2 2], 'newname','VASwitch', 'epochinfo', 'yes');
        pop_saveset(EEG3, 'filename', strcat([filename,'_AV-VisualTargetSwitch.set']), 'filepath', EEG.filepath);
        EEG4 = pop_selectevent( EEG, 'event',Trials_All{idx_trials, 5}' ,'deleteevents','on');
        EEG4 = pop_epoch(EEG4, {'70'}, [-2 2], 'newname','AVSwitch', 'epochinfo', 'yes');
        pop_saveset(EEG4, 'filename', strcat([filename,'_AV-AudioTargetSwitch.set']), 'filepath', EEG.filepath);
        
        
    end
end


