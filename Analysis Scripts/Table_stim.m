%% Get the number of events in each block and returns and save a file for
%each participant separately, but also a general file with the overall
%%number of hits, misses, and FP for each type of bloc and each participant
%%(Total)
function Table_stim(SUBJ)
load('G:\Expé\ANITI\Data\Scripts analyse\Total.mat')
load('G:\Expé\ANITI\Data\Scripts analyse\RT_All.mat')
load('G:\Expé\ANITI\Data\Scripts analyse\Trials_All.mat')
if isempty(find(contains(Trials_All(:,1), SUBJ)))
    if isempty(Total)||(~exist('Total'))
        init = 0;
    else
        init = Total(end,1);
    end
    for s = 1:length(SUBJ)
        EEG = pop_loadset('filename',['sub-',SUBJ{s},'_ses-S001_task-Default_run-001_eeg_HP_reref_ASR_ICAp.set'],...
            'filepath',['G:/Expé/ANITI/Data/sub-',SUBJ{s},'/ses-S001/eeg/']);
        EEG = pop_epoch( EEG, {  '11'  '12'  '13'  '14'  }, [0  180], 'newname', 'AllTrials', 'epochinfo', 'yes'); %5 seconds at the begining for instructions
        
        event = struct2cell(struct('type', {EEG.event(1:end).type}, 'epoch', {EEG.event(1:end).epoch}, 'latency', {EEG.event(1:end).latency}));
        event(1,:,:) = cellfun(@(x)str2num(x), event(1,:,:), 'UniformOutput', false);
        event = squeeze(cell2mat(event))';
        
        Event = cell(1,max(unique(event(:,2))));
        
        for iBloc = 1:max(unique(event(:,2)))
            idx= find(event(:,2)==iBloc);
            Event(1,iBloc) = {event(idx,1)};
        end
        
        [nTrig, Trig] = cellfun(@(x)groupcounts(x), Event, 'UniformOutput', false);
        
        Tableau = [];
        for iBloc = 1:length(Trig)
            tmpEvent = [];
            for iTarg = [60,52,70,50,108,119]
                if isempty(find(Trig{iBloc}==iTarg))
                    nEvent = 0;
                else
                    nEvent = nTrig{iBloc}(find(Trig{iBloc}==iTarg));
                end
                tmpEvent = cat(2,tmpEvent, nEvent);
            end
            Tableau(iBloc,:) = cat(2,iBloc, Trig{iBloc}(1), tmpEvent);
        end
        Tableau = sortrows(Tableau, 2);
        Tableau = array2table(Tableau);
        Tableau.Properties.VariableNames = {'Bloc', 'Type', 'VisTarg', 'VisResp',...
            'AudioTarg','AudioResp','VisEnd','AudioEnd'};
        save(['G:/Expé/ANITI/Data/sub-',SUBJ{s},'/ses-S001/behavior/Trials_num_',...
            SUBJ{s}, '.mat'], 'Tableau')
        
        [event, Total, RT, Trials] = Trial_count(init+s, event, Total);
        
        RT_All = cat(1, RT_All, RT);
        Trials_All = cat(1, Trials_All, [SUBJ{s}, Trials]);
    end
    save 'RT_All' 'RT_All'
    save 'Total' 'Total'
    save 'Trials_All' 'Trials_All'
end
RT_ALL = [];
Acc_ALL = [];
for ncol = 1:size(RT_All,2)
    for nlin = 1:size(RT_All,1)
        RT_ALL =cat(1,RT_ALL, cat(2, RT_All{nlin,ncol}(:,1), nlin*ones(size(RT_All{nlin,ncol},1),1),...
            ncol*ones(size(RT_All{nlin,ncol},1),1)));
    end
end
Acc_ALL = [Total(:,[8,1]),repmat([1:4]',[length(SUBJ),1])];

csvwrite('RT_All.csv', RT_ALL)
csvwrite('Acc_All.csv', Acc_ALL)
%% Average and boxplots
row1 = {'Unimodal' 'Unimodal' 'Switch' 'Switch'};
row2 = {'Visual' 'Auditory' 'Auditory' 'Visual'};
labelArray = [row1; row2];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
figure
subplot(1,2,1)
boxplot(Total(:, 8)*100, Total(:,2))
ax = gca();
ax.XTickLabel = tickLabels;
ax.TickLabelInterpreter = 'tex';
xlabel('Condition');ylabel('Hit rate (%)');
subplot(1,2,2)
boxplot(cell2mat(cellfun(@mean, RT_All, 'UniformOutput', false)))
ax = gca();
ax.XTickLabel = tickLabels;
ax.TickLabelInterpreter = 'tex';
xlabel('Condition');ylabel('Reaction time (s)');

%Scatter subject by subject for RT
s = cell2mat(cellfun(@size, RT_All, 'UniformOutput', false));
figure
for iSubj = 1:size(RT_All,1)
    subplot(4,6,iSubj)
    hold on
    scatter(ones(s(iSubj, 1:2)),RT_All{iSubj,1})
    scatter(2*ones(s(iSubj, 3:4)),RT_All{iSubj,2})
    scatter(3*ones(s(iSubj, 5:6)),RT_All{iSubj,3})
    scatter(4*ones(s(iSubj, 7:8)),RT_All{iSubj,4})
    hold off
end
end

%% Function to count the number of hits misses and late responses
function [event, Total, RT, Trials] = Trial_count(subj,event,Total)
RTAV = [];
RTVA = [];
SwitchAV = 0;
SwitchVA = 0;
NoSwitch = 0;
TotSwitchAV = 0;
TotSwitchVA = 0;
Trials_to_keep_AV = [];
Trials_to_keep_VA = [];
TrueNegAV = 0;
TrueNegVA = 0;
RT_FP_AV = [];
RT_FP_VA = [];
FalsePosAV = 0;
FalsePosVA = 0;
MissSwitchAV = 0;
MissSwitchVA = 0;
SwitchAV = 0;
SwitchVA = 0;
for iCond = [11,12,13,14] %Moyennes
    idx_event = find(ismember(event(:,1), [11,12,13,14,50,52,60,70]));
    tmpEvent = event(idx_event,:);
    idx_Bloc = tmpEvent(find(tmpEvent(:,1)==iCond),2);
    switch iCond
        case 11 %If Visual target
            BlocEvent = tmpEvent(find(and(ismember(tmpEvent(:,2), idx_Bloc),...
                ismember(tmpEvent(:,1), [60,50,52]))), :);
            idx60 = find(BlocEvent(:,1)==60); %Visual target
            idx52 = find(BlocEvent(:,1)==52); %Visual resp
            idx50 = find(BlocEvent(:,1)==50); %Audio resp
            idxHits = idx60(ismember(idx60+1, idx52));
            Late = sum(and([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500>2,[BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500<=8));
            Late_FP = sum([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500>8);
            Hits = sum(ismember(idx60+1, idx52))-(Late+Late_FP);
            Miss = sum(~ismember(idx60+1, idx52))+Late_FP;
            FalsePos60 = sum(~ismember(idx52-1, idx60))+Late_FP;
            FalsePos70 = length(idx50);
            Hit_rate = Hits/length(idx60);
            Miss_rate = Miss/length(idx60);
            Trials_to_keep_V = BlocEvent(idxHits([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500<=2),3);
            RT_Hits_V = {(BlocEvent(idxHits([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500<=2)+1,3)...
                -BlocEvent(idxHits([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500<=2),3))/500};
            TrueNeg = sum(~ismember(idx70+1, idx50));
        case 12 %If Audio target
            BlocEvent = tmpEvent(find(and(ismember(tmpEvent(:,2), idx_Bloc), ismember(tmpEvent(:,1), [70,50,52]))), :);
            idx70 = find(BlocEvent(:,1)==70); %Audio target
            idx52 = find(BlocEvent(:,1)==52); %Visual resp
            idx50 = find(BlocEvent(:,1)==50); %Audio resp
            idxHits = idx70(ismember(idx70+1, idx50));
            Late = sum(and([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500>2,[BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500<=8));
            Late_FP = sum([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500>8);
            Hits = sum(ismember(idx70+1, idx50))-(Late+Late_FP);
            Miss = sum(~ismember(idx70+1, idx50))+Late_FP;
            FalsePos70 = sum(~ismember(idx50-1, idx70))+Late_FP;
            FalsePos60 = length(idx52);
            Hit_rate = Hits/length(idx70);
            Miss_rate = Miss/length(idx70);
            Trials_to_keep_A = BlocEvent(idxHits([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500<=2),3);
            RT_Hits_A = {(BlocEvent(idxHits([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500<=2)+1,3)...
                -BlocEvent(idxHits([BlocEvent(idxHits+1,3)-BlocEvent(idxHits,3)]./500<=2),3))/500};
            TrueNeg = sum(~ismember(idx60+1, idx52));
        case {13,14} %
            BlocEvent60 = tmpEvent(find(and(ismember(tmpEvent(:,2), idx_Bloc), ismember(tmpEvent(:,1), [60,50,52]))), :);
            BlocEvent70 = tmpEvent(find(and(ismember(tmpEvent(:,2), idx_Bloc), ismember(tmpEvent(:,1), [70,50,52]))), :);
            
            idx60 = find(BlocEvent60(:,1)==60); %Visual target
            idx70 = find(BlocEvent70(:,1)==70); %Audio target
            idx52_60 = find(BlocEvent60(:,1)==52); %Visual resp
            idx50_60 = find(BlocEvent60(:,1)==50); %Audio resp
            idx52_70 = find(BlocEvent70(:,1)==52); %Visual resp
            idx50_70 = find(BlocEvent70(:,1)==50); %Audio resp
            Resp60 = BlocEvent60(idx60(ismember(idx60+1, idx52_60)),[1,3]);
            Resp52_60 = BlocEvent60(idx60(ismember(idx60+1, idx52_60))+1,[1,3]); %Visual hits
            RT60 = Resp52_60(:,2)-Resp60(:,2);
            Resp70 = BlocEvent70(idx70(ismember(idx70+1, idx50_70)),[1,3]);
            Resp50_70 = BlocEvent70(idx70(ismember(idx70+1, idx50_70))+1,[1,3]); %Audio hits
            RT70 = Resp50_70(:,2)-Resp70(:,2);
            Resp_AV = sortrows(cat(2,cat(1,Resp60,Resp70),cat(1,RT60/500,RT70/500)),2); %Trigger|tTrigger(nb points)|RT
            tmpResp = Resp_AV;
            Resp_AV = Resp_AV(find(Resp_AV(:,3)<=2),:);
            Late = size(tmpResp(find(and(tmpResp(:,3)>2, tmpResp(:,3)<=6))),1);
            if sum(~diff(tmpEvent(find(ismember(tmpEvent(:,3),tmpResp(tmpResp(:,3)>6, 2)))+1,2)))~=0
                MissSwitchVA = MissSwitchVA+sum(tmpResp(~diff(tmpEvent(find(ismember(tmpEvent(:,3),...
                    tmpResp(tmpResp(:,3)>6, 2)))+1,2)))==60);
                MissSwitchAV = MissSwitchAV+sum(tmpResp(~diff(tmpEvent(find(ismember(tmpEvent(:,3),...
                tmpResp(tmpResp(:,3)>6, 2)))+1,2)))==70);
            end
                    
            iFalsePos60 = BlocEvent60(idx52_60(~ismember(idx52_60-1, idx60)),[1,3]); %Visual FP
            FalsePos60 = size(iFalsePos60,1);
            iFalsePos70 = BlocEvent70(idx50_70(~ismember(idx50_70-1, idx70)),[1,3]); %Audio FP
            FalsePos70 = size(iFalsePos70,1);
            Altern = tmpEvent(find(and(ismember(tmpEvent(:,2), idx_Bloc), ismember(tmpEvent(:,1), [60,70]))), :);
            
            %Bining into True positive (hits/Switch), True Neg (No resp),
            %False positive (Resp wrong) and False Neg (No resp wrong) for
            %audio to visual (AV) and visual to audio (VA)
            for i = 1:length(idx_Bloc)
                if iCond == 13
                    first = min(find(and(Altern(:,1)==60, Altern(:,2)==idx_Bloc(i))));
                    last = max(find(Altern(:,2)==idx_Bloc(i)));
                    for iSwitch = first:last-1
                        switch Altern(iSwitch,1)
                            case 60
                                switch Altern(iSwitch+1,1)
                                    case 70
                                        TotSwitchVA = TotSwitchVA+1;
                                        if ismember(Altern(iSwitch,3), Resp_AV(:,2))
                                            SwitchVA = SwitchVA+1;
                                            RTVA = cat(1, RTVA, Resp_AV(Altern(iSwitch,3)==Resp_AV(:,2),:));
                                        else
                                            MissSwitchVA = MissSwitchVA+1;
                                        end
                                    case 60
                                        if ~ismember(Altern(iSwitch,3), Resp_AV(:,2))
                                            TrueNegVA = TrueNegVA+1;
                                        else
                                            FalsePosVA = FalsePosVA+1;
                                            RT_FP_VA = cat(1, RT_FP_VA, Resp_AV(Altern(iSwitch,3)==Resp_AV(:,2),:));
                                        end
                                end
                            case 70
                                switch Altern(iSwitch+1,1)
                                    case 60
                                        TotSwitchAV = TotSwitchAV+1;
                                        if ismember(Altern(iSwitch,3), Resp_AV(:,2))
                                            SwitchAV = SwitchAV+1;
                                            RTAV = cat(1, RTAV, Resp_AV(Altern(iSwitch,3)==Resp_AV(:,2),:));
                                        else
                                            MissSwitchAV = MissSwitchAV+1;
                                        end
                                    case 70
                                        if ~ismember(Altern(iSwitch,3), Resp_AV(:,2))
                                            TrueNegAV = TrueNegAV+1;
                                        else
                                            FalsePosAV = FalsePosAV+1;
                                            RT_FP_AV = cat(1, RT_FP_AV, Resp_AV(Altern(iSwitch,3)==Resp_AV(:,2),:));
                                        end
                                end
                        end
                    end
                    
                elseif iCond == 14
                    first = min(find(and(Altern(:,1)==70, Altern(:,2)==idx_Bloc(i))));
                    last = max(find(Altern(:,2)==idx_Bloc(i)));
                    for iSwitch = first:last-1
                        switch Altern(iSwitch,1)
                            case 60
                                switch Altern(iSwitch+1,1)
                                    case 70
                                        TotSwitchVA = TotSwitchVA+1;
                                        if ismember(Altern(iSwitch,3), Resp_AV(:,2))
                                            SwitchVA = SwitchVA+1;
                                            RTVA = cat(1, RTVA, Resp_AV(Altern(iSwitch,3)==Resp_AV(:,2),:));
                                        else
                                            MissSwitchVA = MissSwitchVA+1;
                                        end
                                    case 60
                                        if ~ismember(Altern(iSwitch,3), Resp_AV(:,2))
                                            TrueNegVA = TrueNegVA+1;
                                        else
                                            FalsePosVA = FalsePosVA+1;
                                            RT_FP_VA = cat(1, RT_FP_VA, Resp_AV(Altern(iSwitch,3)==Resp_AV(:,2),:));
                                        end
                                end
                            case 70
                                switch Altern(iSwitch+1,1)
                                    case 60
                                        TotSwitchAV = TotSwitchAV+1;
                                        if ismember(Altern(iSwitch,3), Resp_AV(:,2))
                                            SwitchAV = SwitchAV+1;
                                            RTAV = cat(1, RTAV, Resp_AV(Altern(iSwitch,3)==Resp_AV(:,2),:));
                                        else
                                            MissSwitchAV = MissSwitchAV+1;
                                        end
                                    case 70
                                        if ~ismember(Altern(iSwitch,3), Resp_AV(:,2))
                                            TrueNegAV = TrueNegAV+1;
                                        else
                                            FalsePosAV = FalsePosAV+1;
                                            RT_FP_AV = cat(1, RT_FP_AV, Resp_AV(Altern(iSwitch,3)==Resp_AV(:,2),:));
                                        end
                                end
                        end
                    end
                end
            end
            
            

            Miss_rate = (MissSwitchVA+MissSwitchAV)/(TotSwitchAV+TotSwitchVA);
            Trials_to_keep_VA = cat(1,Trials_to_keep_VA, RTVA(:,2));
            Trials_to_keep_AV = cat(1,Trials_to_keep_AV, RTAV(:,2));
    end
    Total = cat(1, Total, cat(2,subj, iCond, Hits, Miss, FalsePos60, FalsePos70, Late, Hit_rate, Miss_rate));
end

Total(and(Total(:,2)==13, Total(:,1)==subj), [3:9]) =  [size(RTVA,1), MissSwitchVA, size(RT_FP_VA,1),...
    size(RT_FP_AV,1), Late, SwitchVA/TotSwitchVA, MissSwitchVA/TotSwitchVA];%VisualToAudio
Total(and(Total(:,2)==14, Total(:,1)==subj), [3:9]) = [size(RTAV,1), MissSwitchAV, size(RT_FP_VA,1),...
    size(RT_FP_AV,1), Late, SwitchAV/TotSwitchAV, MissSwitchAV/TotSwitchAV]; %AudioToVisual
RT = [RT_Hits_V, RT_Hits_A, {RTVA(:,3)}, {RTAV(:,3)}];
Trials_to_keep_V = find(and(ismember(event(:,3), Trials_to_keep_V), ismember(event(:,1),60)));
Trials_to_keep_A = find(and(ismember(event(:,3), Trials_to_keep_A), ismember(event(:,1),70)));
Trials_to_keep_VA = find(and(ismember(event(:,3), Trials_to_keep_VA), ismember(event(:,1),60)));
Trials_to_keep_AV = find(and(ismember(event(:,3), Trials_to_keep_AV), ismember(event(:,1),70)));
Trials = {Trials_to_keep_V, Trials_to_keep_A, Trials_to_keep_VA, Trials_to_keep_AV};
end