% runs the stimlus sequence for classical audio oddball task
% it takes as arguments
% seq = the stimulus sequence based on categories created previously
% listReg = the list of regular stimlui
% listOdd = the list of odd stimuli
% the available list for this task are below, which are defined previously
% listHi = {'hi1' 'hi2'};
% listMe = {'me1' 'me2'};
% listLo = {'lo1' 'lo2'};
% date = date of today (yyyymmdd) specified previouisly

% the following arguments are general
% pl is a placeholder argument for getSequence() to proceed with operations
%   specific for classicalAud
%
% function [timingData,taskNames] = classicalAud(seq, listReg, listOdd, listOddEx, ...
%     dateString, window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames, pl)

function [timingData,taskNames] = classicalAud(seq, listReg, listOdd, listOddEx, ...
    dateString, screen, allCoords, lineWidthPix, taskNames, pl, devType)
window = screen.win;
white = screen.white;
xCenter = screen.xCenter;
yCenter = screen.yCenter;

% Draw loading instruction to wait for sequence creation
line = 'loading ...';

DrawFormattedText(window, line, 'center', 'center', white);
Screen('Flip', window);

% Create sequence
finalSequence = getSequence(seq, listReg, listOdd, listOddEx, pl);

% Instructions for the task
line1 = 'In this task, you will listen to a series of sound.';
line2 = '\n\n Please pay attention to the characteristics of the sound.';
line3 = '\n\n You may be asked to answer relevant questions at the end of this task.';
line4 = '\n\n Press [space] to continue.';

% Draw instructions
DrawFormattedText(window, [line1 line2 line3 line4], 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

% Draw fixation cross
Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter]); %, 2);
Screen('Flip', window);

% Determine regular stimlus by the input
if ismember('lo1', listReg)
    regFreq = 220;
elseif ismember('me1', listReg)
    regFreq = 440;
elseif ismember('hi1', listReg)
    regFreq = 880;
else
    error('listReg contains unexpected values');
end

% Determine odd stimulus by the output
if ismember('lo1', listOdd)
    oddFreq = 220;
elseif ismember('me1', listOdd)
    oddFreq = 440;
elseif ismember('hi1', listOdd)
    oddFreq = 880;
else
    error('listReg contains unexpected values');
end

% Determine oddEx stimulus by the output
if ismember('lo1', listOddEx)
    oddExFreq = 220;
elseif ismember('me1', listOddEx)
    oddExFreq = 440;
elseif ismember('hi1', listOddEx)
    oddExFreq = 880;
else
    error('listReg contains unexpected values');
end

% initialize data structure to store timing data
timingData = struct();

WaitSecs(2);

% -----------------!!!send trigger for starting!!!-----------------
if strcmp(devType, 'EEG')
    write(port, 1,"uint8");
elseif strcmp(devType, 'MEG')
    PTBSendTrigger(1,0);
else
    Beeper(2000)
end

% Record the start time of the experiment
startTime = GetSecs();

% running the stimlus sequence
for i = 1:size(finalSequence, 1)
    stiType = finalSequence{i, 1};
    stiName = finalSequence{i, 2};

    onsetTime = GetSecs() - startTime;

    if strcmp(stiType, 'A')
        current_code = 64;
        Beeper(regFreq,0.6,0.45);
    else
        if rand() < 0.5
            oddFreqPresent = oddFreq;
        else
            oddFreqPresent = oddExFreq;
        end
        current_code = 128;
        Beeper(oddFreqPresent,0.6,0.45);
    end

    offsetTime = GetSecs() - startTime;

    % -----------------!!!send trigger for ending!!!-----------------
    if strcmp(devType, 'EEG')
        write(port, current_code,"uint8");
    elseif strcmp(devType, 'MEG')
        PTBSendTrigger(current_code,0);
    else
        Beeper(2000)
    end
    % write(port, current_code,"uint8");
    % Beeper(2000)
    disp(stiName)
    % store
    timingData(i).stiType = stiType;
    timingData(i).stiName = stiName;
    timingData(i).onsetTime = onsetTime;
    timingData(i).offsetTime = offsetTime;

    % Randomize inter-stimulus interval between 750ms to 1250ms
    WaitSecs(0.75 + rand() * 0.5);

end

timingData(1).startTime = startTime;

regObject = listReg{1}(1:2);

filename = sprintf('%s_timingData_%s_%s.mat', dateString, 'caud', regObject);

taskNames{end+1} = filename;

save(filename, 'timingData');

line = 'This task is over. Please press [space] to proceed to the next task.';

DrawFormattedText(window, line, 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

end % end of function