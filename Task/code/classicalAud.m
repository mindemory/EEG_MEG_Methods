%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is created by Brian Yan (by2139@nyu.edu)
% And has been adapted for this course.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

function [timingData,taskNames] = classicalAud(cat, seq, labels, ...
    dateString, window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames, devType, port)

% Draw loading instruction to wait for sequence creation
line = 'loading ...';

DrawFormattedText(window, line, 'center', 'center', white);
Screen('Flip', window);


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
if strcmp(devType, 'EEG')
    Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter], 2);
elseif strcmp(devType, 'MEG')
    Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter], 2);
else
    Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter]);
end
Screen('Flip', window);


% initialize data structure to store timing data
timingData = struct();

KbStrokeWait;

% -----------------!!!send trigger for starting!!!-----------------
% write(port, 1,"uint8");
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
for i = 1:numel(seq)
    stim = seq{i};
    label = labels{i};



    if stim == "n200"
        stimulus = 200;
    elseif stim == "n300"
        stimulus = 300;
    elseif stim == "n400"
        stimulus = 400;
    elseif stim == "n500"
        stimulus = 500;
    elseif stim == "n600"
        stimulus = 600;
    elseif stim == "n700"
        stimulus = 700;
    elseif stim == "n800"
        stimulus = 800;
    elseif stim == "n900"
        stimulus = 900;
    end


    % -----------------!!!send trigger for stimulus onset!!!-----------------


    if label == "reg"
        current_code = 64;
    elseif label == "odd"
        current_code = 128;
    end

    if strcmp(devType, 'EEG')
        write(port, current_code,"uint8");
    elseif strcmp(devType, 'MEG')
        PTBSendTrigger(current_code,0);
    else
        Beeper(2000)
    end


    onsetTime = GetSecs() - startTime;

    if strcmp(devType, 'EEG')
        Beeper(stimulus, 0.45); % needs to be tested
    elseif strcmp(devType, 'MEG')
        Beeper(stimulus, 0.06, 0.45); % needs to be tested
    else
        Beeper(stimulus, 0.5, 0.45);
    end



    offsetTime = GetSecs() - startTime;

    % -----------------!!!send trigger for ending!!!-----------------
    % write(port, current_code,"uint8");
    % Beeper(2000)

    % store
    timingData(i).stiType = label;
    timingData(i).stiName = stim;
    timingData(i).onsetTime = onsetTime;
    timingData(i).offsetTime = offsetTime;

    % Randomize inter-stimulus interval between 750ms to 1250ms
    pause(0.75 + rand() * 0.5);

end

timingData(1).startTime = startTime;
dateStringBlah = datestr(now, 'yyyymmdd_HHMMSS');

filename = sprintf('%s_timingData_%s_%s.mat', dateStringBlah, 'caud', cat);

taskNames{end+1} = filename;
dirToSave = '../../../TaskTiming/';
if ~exist("dirToSave", 'dir')
    mkdir(dirToSave)
end
filename = [dirToSave filename];

save(filename, 'timingData');

line = 'This task is over. Please press [space] to proceed to the next task.';

DrawFormattedText(window, line, 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

end % end of function



