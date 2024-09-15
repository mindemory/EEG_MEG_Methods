% runs the stimlus sequence for semantic visual oddball task
% it takes as arguments
% seq = the stimulus sequence based on categories created previously
% listReg = the list of regular stimlui
% listOdd = the list of odd stimuli
% the available list for this task are below, which are defined previously
% animalsListV = [...]
% toolsListV = [...]
% bodyPartsListV = [...], n = 20
% date = date of today (yyyymmdd) specified previouisly

% the following arguments are general

% function [timingData, taskNames] = semanticVis(seq, listReg, listOdd, listOddEx, ...
%     dateString, window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames)

function [timingData, taskNames] = semanticVis(seq, listReg, listOdd, listOddEx, ...
    dateString, screen, allCoords, lineWidthPix, taskNames, devType)

window = screen.win;
white = screen.white;
xCenter = screen.xCenter;
yCenter = screen.yCenter;

% Draw loading instruction to wait for sequence creation
line = 'loading ...';

DrawFormattedText(window, line, 'center', 'center', white);
Screen('Flip', window);

% Create sequence
finalSequence = getSequence(seq, listReg, listOdd, listOddEx);

line1 = 'In this task, you will be presented with a set of words.';
line2 = '\n\n Please pay attention to the words presented to you, and try to';
line3 = '\n\n identify the categories the words fall under. Please bear';
line4 = '\n\n in mind that the words may come from more than one category.';
line5 = '\n\n After the presentation is finished, you will be asked to';
line6 = '\n\n report the categories you identified. Press [space] to continue.';

% Draw instructions
DrawFormattedText(window, [line1 line2 line3 line4 line5 line6], 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

% Draw fixation cross
Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter]); %, 2);
Screen('Flip', window);

WaitSecs(1);

% Empty window
Screen('FillRect', window, white / 2)
Screen('Flip', window);

WaitSecs(1);

% initialize data structure to store timing data
timingData = struct();

% -----------------!!!send trigger for starting!!!-----------------
if strcmp(devType, 'EEG')
    write(port, 2,"uint8");
elseif strcmp(devType, 'MEG')
    PTBSendTrigger(2,0);
else
    Beeper(2000)
end
% write(port, 2,"uint8");
% Beeper(2000)

% Record the start time of the experiment
startTime = GetSecs();

% running the stimlus sequence
for i = 1:size(finalSequence, 1)
    stiType = finalSequence{i, 1};
    stiName = finalSequence{i, 2};

    % onsetTime = GetSecs() - startTime;

    if strcmp(stiType, 'A')
        current_code = 64;
        stimulus = stiName;
    else
        current_code = 128;
        stimulus = stiName;
    end

    % Draw stimulus
    DrawFormattedText(window, stimulus, 'center', 'center', white);

    % Flip to the screen
    onsetTime = Screen('Flip', window);

    % Display word for 1 second and then offset word
    WaitSecs(1);

    Screen('FillRect', window, white / 2);

    offsetTime = Screen('Flip', window);

    onsetTime = onsetTime - startTime;
    offsetTime = offsetTime - startTime;

    disp(stimulus)

    % -----------------!!!send trigger for ending!!!-----------------
    if strcmp(devType, 'EEG')
        write(port, current_code,"uint8");
    elseif strcmp(devType, 'MEG')
        PTBSendTrigger(current_code,0);
    else
        Beeper(2000)
    end
    % write(port, current_code,"uint8");
%     Beeper(2000)

    % store
    timingData(i).stiType = stiType;
    timingData(i).stiName = stiName;
    timingData(i).onsetTime = onsetTime;
    timingData(i).offsetTime = offsetTime;

    % Randomize inter-stimulus interval between 750ms to 1250ms
    WaitSecs(0.5 + rand() * 0.5);

end

timingData(1).startTime = startTime;

regObject = inputname(2);

filename = sprintf('%s_timingData_%s_%s.mat', dateString, 'svis', regObject);

taskNames{end+1} = filename;

save(filename, 'timingData');


line1 = 'This is the end of the task. Please put down the categories';
line2 = '\n\n you identified from the task. Press [space] to continue when you finish.';

DrawFormattedText(window, [line1 line2], 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

end % end of function