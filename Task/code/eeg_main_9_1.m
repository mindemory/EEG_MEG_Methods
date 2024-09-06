% Initialize Psychtoolbox and set up the experiment environment
sca; % Close all screens
close all;
clear all;
clearvars;

Screen('Preference', 'SkipSyncTests', 1)

%---------------!!!!!!!!mods required!!!!!!----------------------%
% add "port" to function input for trigger sending
% figure out stimulus choice
%   odd words from two categories check----------------
% discuss stimulus presentation duration
%   try to unify the trial length (stimulus + gap) check--------
% instructions for recall?
%   no.
% make sure that a word does not repeat too soon
% check, now by 5--------------------
% for each task in each section, stimulus needs to be manually inputted
% the recorded file name is date+name+task+regular stimulus

% average display time for semantic audio stimulus is 0.457740545
% average gap time is 1 sec
% 

%---------trigger code-------------------------------------
% addpath('triggers');
% trigger code
% classicalau
% start = 1
% reg = 64, odd = 128

% visual semantic
% start = 2
% reg = 64, odd = 128

% audio semantic
% start = 4
% reg = 64, odd = 128

% story
% start = 8
% onset = 64, offset = 128
%-----------trigger code----------------------------

% get the date of the experiment to add to file name later
date = datestr(now, 'yyyymmdd_HHMMSS');
date = date(1:8);

% ---the following section imports all the audio files and store them
% under corresponding variables names to be used later

% curr_dir = pwd;

curr_dir = "D:\Desktop\MemoryLab";

% port = init_trigger;




%--auditory semantic oddball setup start------------------------------

listHi = {'hi1' 'hi2'};
listMe = {'me1' 'me2'};
listLo = {'lo1' 'lo2'};

%--auditory semantic oddball setup end------------------------------


disp('import list');

%--visual semantic oddball setup start------------------------------

% read the file that contains all the word
% wordTable = readtable([curr_dir filesep 'EventModel' filesep 'selected_categories.xlsx']);
% Define the file path
file_path = fullfile(curr_dir, 'stimuliSemantic', 'selected_categories.xlsx');

% Read the table
wordTable = readtable(file_path);


% slice the matrix into lists
animalsListV = wordTable.Animals;
toolsListV = wordTable.Tools;
bodyPartsListV = wordTable.BodyParts;

%-------visual semantic oddball setup end-----------------------------





disp('import audio files');

%-------auditory semantic oddball setup start---------------------------

% Define word sets directories
% animalsDir = [curr_dir filesep 'EventModel' filesep 'animals'];
% toolsDir = [curr_dir filesep 'EventModel' filesep 'tools'];
% bodyPartsDir = [curr_dir filesep 'EventModel' filesep 'bodyParts'];
% vehiclesDir = [curr_dir filesep 'EventModel' filesep 'vehicles'];

% animalsDir = [curr_dir filesep 'stimuliSemantic' filesep 'animals'];
% toolsDir = [curr_dir filesep 'stimuliSemantic' filesep 'tools'];
% bodyPartsDir = [curr_dir filesep 'stimuliSemantic' filesep 'bodyParts'];
% vehiclesDir = [curr_dir filesep 'stimuliSemantic' filesep 'vehicles'];

animalsDir = fullfile(curr_dir,'stimuliSemantic','animals');
toolsDir = fullfile(curr_dir,'stimuliSemantic','tools');
bodyPartsDir = fullfile(curr_dir,'stimuliSemantic','bodyParts');
vehiclesDir = fullfile(curr_dir,'stimuliSemantic','vehicles');

% Load audio files
% dir() shows all files in the folder ('*.mp3' restricts the search to .mp3 files)
% fullfile() makes sure that the format works across platforms

animalsAu = dir(fullfile(animalsDir, '*.mp3'));
toolsAu = dir(fullfile(toolsDir, '*.mp3'));
bodyPartsAu = dir(fullfile(bodyPartsDir, '*.mp3'));
vehiclesAu = dir(fullfile(vehiclesDir, '*.mp3'));

% Initialize structures to store audio data
animalsData = struct();
toolsData = struct();
bodyPartsData = struct();
vehiclesData = struct();

% Read audio files for animals and store in structure with file names as field names
for i = 1:numel(animalsAu)
    filePath = fullfile(animalsDir, animalsAu(i).name);
    [audioData, fs] = audioread(filePath);
    [~, fileName, ~] = fileparts(animalsAu(i).name);
    animalsData.(fileName) = audioData;
    %disp(['Loaded: ' animalsAu(i).name ' as ' fileName]);
end

% Read audio files for tools and store in structure with file names as field names
for i = 1:numel(toolsAu)
    filePath = fullfile(toolsDir, toolsAu(i).name);
    [audioData, fs] = audioread(filePath);
    [~, fileName, ~] = fileparts(toolsAu(i).name);
    toolsData.(fileName) = audioData;
    
end

% Read audio files for bodyParts and store in structure with file names as field names
for i = 1:numel(bodyPartsAu)
    filePath = fullfile(bodyPartsDir, bodyPartsAu(i).name);
    [audioData, fs] = audioread(filePath);
    [~, fileName, ~] = fileparts(bodyPartsAu(i).name);
    bodyPartsData.(fileName) = audioData;
    
end

% Read audio files for vehicles and store in structure with file names as field names
for i = 1:numel(vehiclesAu)
    filePath = fullfile(vehiclesDir, vehiclesAu(i).name);
    [audioData, fs] = audioread(filePath);
    [~, fileName, ~] = fileparts(vehiclesAu(i).name);
    vehiclesData.(fileName) = audioData;
    
end

% fieldnames gets all the names in the data structure
animalsListA = fieldnames(animalsData);
toolsListA = fieldnames(toolsData);
bodyPartsListA = fieldnames(bodyPartsData);
vehiclesListA = fieldnames(vehiclesData);

%-------auditory semantic oddball setup start---------------------------



%------story setup begin-----------------------------------------------

% Read the audio file for the story
% The story is sliced by the events, along with a cue "next passage" for
% the cued condition
% [storyAu, freq] = audioread([curr_dir filesep 'EventModel' filesep 'SecretLife.mp3']);

% [storyAu, freq] = audioread(fullfile(curr_dir, 'stimuliSemantic', 'test.mp3'));

storyDir = fullfile(curr_dir,'stimuliSemantic','story');
storyAu = dir(fullfile(storyDir, '*.mp3'));
storyData = struct();
for i = 1:numel(storyAu)
    filePath = fullfile(storyDir, storyAu(i).name);
    [audioData, fs] = audioread(filePath);
    [~, fileName, ~] = fileparts(storyAu(i).name);
    storyData.(fileName) = audioData;
    %disp(['Loaded: ' animalsAu(i).name ' as ' fileName]);
end
%-------story setup end----------------------------------------------

disp('import finished');


%-------visual display setup start--------------------------------------

% Setup Psychtoolbox
PsychDefaultSetup(2);
InitializePsychSound(1);

% Define screen numbers
screenNumbers = Screen('Screens');
secondaryScreen = max(screenNumbers);

% Define black and white
white = WhiteIndex(secondaryScreen);
black = BlackIndex(secondaryScreen);
grey = white / 2;

% Open an on screen window on the secondary monitor
[window, windowRect] = PsychImaging('OpenWindow', secondaryScreen, grey);

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Set waitframes to 1
waitframes = 1;

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Setup the text type for the window
%Screen('TextFont', window, 'Ariel');
Screen('TextSize', window, 48);

% Setup keyboard
KbName('UnifyKeyNames');

% Set up the fixation cross
fixCrossDimPix = 40;
xCoords = [-fixCrossDimPix fixCrossDimPix 0 0];
yCoords = [0 0 -fixCrossDimPix fixCrossDimPix];
allCoords = [xCoords; yCoords];
lineWidthPix = 4;

disp('visual setup finished')

%-----------visual display setup end---------------------------------

% Define audio device
audioDevice = PsychPortAudio('Open', [], 1, 1, [], 2);

% Set playback volume to 60%
PsychPortAudio('Volume', audioDevice, 0.6);

disp('audio setup finished')

% ----------basic audio setup end---------------------------------------


% Generate stimulus sequence
% The repmat function replicates and tiles an array. In this case, 
% it’s creating an array that consists of the string ‘A’ repeated numA 
% times. The {} around ‘A’ is used to create a cell array containing the string ‘A’. The 1 
% get the sequence representing two categories
numA = 5 * 4; % Each word in set A is played 4 times
numB = round(0.2 / 0.8 * numA); % Calculate number of B words to maintain 0.8 to 0.2 ratio
sequence = [repmat({'A'}, 1, numA) repmat({'B'}, 1, numB)];
sequence = sequence(randperm(length(sequence))); % Shuffle sequence

% define the function handles
tasks = {@classicalAud, @semanticVis, @semanticAud};

% initialize cell array to store the name and order of the task
taskNames = {};

disp('main setup finished')

%--------------------ALL SETUP FINISHED------------------------------

%-----------sequence testing-------------
% finalSequence = getSequence(sequence, animalsListA, bodyPartsListA, toolsListA);
















%--------------------START DISPLAY-----------------------------------


% Instructions
line = 'Press any key to initiate experiment. This may take a few seconds.';

% Draw instructions
DrawFormattedText(window, line, 'center', 'center', white);

% Flip to the screen
Screen('Flip', window);

KbStrokeWait;


% -----------display testing---------------------------------
% timingData = struct();
% 
% [timingData, taskNames] = semanticAud(sequence, animalsListA, bodyPartsListA, toolsListA, date, ...
%     window, white, allCoords, lineWidthPix, xCenter, yCenter, audioDevice, taskNames);
% [timingDataAll{i}, taskNames] = semanticVis(sequence, animalsListV, toolsListV, bodyPartsListV, date, ...
%     window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames);
% [timingDataAll{i}, taskNames] = classicalAud(sequence, listHi, listLo, listMe, date, ...
%     window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames, 1);
% % 
% sca;
% 










%-------------------Section I---------------------------------------

% generate order of tasks which returns e.g. [3 1 2]
order = randperm(length(tasks));

% intialize cell array that stores the order
taskOrder = cell(1,length(tasks));

% Initialize a cell array to store the timing data
timingDataAll = cell(1, length(tasks));

% try using struct
% timingDataAll = struct();

% the display of the three oddball tasks are organized as functions. 
% the display orders in a section is randomized
for i = 1:length(tasks)
    taskIndex = order(i);
    task = tasks{taskIndex};

    if taskIndex == 1
        [timingDataAll{i}, taskNames] = task(sequence, listHi, listMe, listLo, date, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames, 1);
    elseif taskIndex == 2
        [timingDataAll{i}, taskNames] = task(sequence, animalsListV, toolsListV, bodyPartsListV, date, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames);
    elseif taskIndex == 3
        [timingDataAll{i}, taskNames] = task(sequence, toolsListA, bodyPartsListA, animalsListA, date, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, audioDevice, taskNames);
    end

end    


%------------------end of Section I------------------------------






















% ----------Section II------------------------------------------------

line1 = 'In the following task, you will listen to a story.';
line2 = '\n\n Press [space] to continue whenever you are ready.';

% Draw instructions
DrawFormattedText(window, [line1 line2], 'center', 'center', white);

% Flip to the screen
Screen('Flip', window);


% Define audio device (no need since no clean up in the previous section)
% audioDevice = PsychPortAudio('Open', [], 1, 1, [], 2);

% Set playback volume to 60%
PsychPortAudio('Volume', audioDevice, 0.6);


% PsychPortAudio('FillBuffer', audioDevice, [storyAu']);

% Set the cell array that contains the sequence of audio to be played
%   (the field name for the cue audio is "Cue")
storySeq = {'Marine' 'Main1' 'Medical' 'Main2' 'Trial' 'Main3' 'War' 'Main4'};

% initialize ds to store timing structure
timingData2 = struct();

% Press any key to continue
KbStrokeWait;

% Draw fixation cross
Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter], 2);
Screen('Flip', window);

%--! send trigger ! --

% write(port, 8,"uint8");
Beeper(2000);

%--! send trigger ! --

startTime = GetSecs;

% Start audio playback
%[, repetitions=1] [, when=0] [, waitForStart=1] (so it records the actual
% onset time)

for i = 1:numel(storySeq)

% get event name
event = storySeq{i};

% get audio
audioData = storyData.(event);

PsychPortAudio('FillBuffer', audioDevice, audioData');
    
% write(port, current_code,"uint8");
    
% Start audio and record onset time
onsetTime = PsychPortAudio('Start', audioDevice, 1, 0, 1);
    

%--! send trigger ! --

% write(port, 64,"uint8");
Beeper(2000);

%--! send trigger ! --

% Wait for the audio to finish
PsychPortAudio('Stop', audioDevice, 1, 1);


% get offset time
offsetTime = GetSecs();

%--! send trigger ! --

% write(port, 128,"uint8");
Beeper(2000);

%--! send trigger ! --

onsetTime = onsetTime - startTime;
offsetTime = offsetTime - startTime;


% onsetTime = PsychPortAudio('Start', audioDevice, 1, 0, 1);
% 
% % Wait for the audio to finish
% offsetTime = PsychPortAudio('Stop', audioDevice, 1, 1);
% 

% Record duration time
duration = offsetTime-onsetTime;

timingData2(i).event = event;
timingData2(i).onsetTime = onsetTime;
timingData2(i).offsetTime = offsetTime;
timingData2(i).duration = duration;

end

filename = sprintf('%s_timingData_%s.mat', date, 'story1');

% Save timing data to a .mat file
save(filename, 'timingData2');

line = 'This task is over. Please press [space] to proceed to the next task.';

DrawFormattedText(window, line, 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

sca;

%------------end of Section II--------------------------------------



























%-------------------Section III----------------------------------

% generate order of tasks which returns e.g. [3 1 2]
order = randperm(length(tasks));

% intialize cell array that stores the order
taskOrder = cell(1,length(tasks));

% Initialize a cell array to store the timing data
timingDataAll = cell(1, length(tasks));

% try using struct
% timingDataAll = struct();


for i = 1:length(tasks)
    taskIndex = order(i);
    task = tasks{taskIndex};

    if taskIndex == 1
        [timingDataAll{i}, taskNames] = task(sequence, listMe, listLo, listHi, date, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames, 1);
    elseif taskIndex == 2
        [timingDataAll{i}, taskNames] = task(sequence, toolsListV, animalsListV, bodyPartsListV, date, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames);
    elseif taskIndex == 3
        [timingDataAll{i}, taskNames] = task(sequence, bodyPartsListA, animalsListA, toolsListA, date, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, audioDevice, taskNames);
    end
    
end    

%-----------------end of Section III---------------------------


















%------------------start of Section IV-----------------------------------

line1 = 'In the following task, you will listen to the story again.';
line2 = '\n\n Press [space] to continue whenever you are ready.';

% Draw instructions
DrawFormattedText(window, [line1 line2], 'center', 'center', white);

% Flip to the screen
Screen('Flip', window);




% Define audio device (no need since no clean up in the previous section)
% audioDevice = PsychPortAudio('Open', [], 1, 1, [], 2);

% Set playback volume to 60%
PsychPortAudio('Volume', audioDevice, 0.6);


% PsychPortAudio('FillBuffer', audioDevice, [storyAu']);

% Set the cell array that contains the sequence of audio to be played
%   (the field name for the cue audio is "Cue")
storySeq = {'Marine' 'Main1' 'Medical' 'Main2' 'Trial' 'Main3' 'War' 'Main4'};

% initialize ds to store timing structure
timingData4 = struct();

% Press any key to continue
KbStrokeWait;

% Draw fixation cross
Screen('DrawLines', window, allCoords, lineWidthPix, white, [xCenter yCenter], 2);
Screen('Flip', window);

%--! send trigger ! --

% write(port, 8,"uint8");
Beeper(2000);

%--! send trigger ! --

startTime = GetSecs;

% Start audio playback
%[, repetitions=1] [, when=0] [, waitForStart=1] (so it records the actual
% onset time)

for i = 1:numel(storySeq)

% get event name
event = storySeq{i};

% get audio
audioData = storyData.(event);

PsychPortAudio('FillBuffer', audioDevice, audioData');
    
% write(port, current_code,"uint8");
    
% Start audio and record onset time
onsetTime = PsychPortAudio('Start', audioDevice, 1, 0, 1);
    

%--! send trigger ! --

% write(port, 64,"uint8");
Beeper(2000);

%--! send trigger ! --

% Wait for the audio to finish
PsychPortAudio('Stop', audioDevice, 1, 1);


% get offset time
offsetTime = GetSecs();

%--! send trigger ! --

% write(port, 128,"uint8");
Beeper(2000);

%--! send trigger ! --

onsetTime = onsetTime - startTime;
offsetTime = offsetTime - startTime;


% onsetTime = PsychPortAudio('Start', audioDevice, 1, 0, 1);
% 
% % Wait for the audio to finish
% offsetTime = PsychPortAudio('Stop', audioDevice, 1, 1);
% 

% Record duration time
duration = offsetTime-onsetTime;

timingData4(i).event = event;
timingData4(i).onsetTime = onsetTime;
timingData4(i).offsetTime = offsetTime;
timingData4(i).duration = duration;

end


filename = sprintf('%s_timingData_%s.mat', date, 'story2');

% Save timing data to a .mat file
save(filename, 'timingData4');

line = 'This task is over. Please press [space] to proceed to the next task.';

DrawFormattedText(window, line, 'center', 'center', white);

Screen('Flip', window);

KbStrokeWait;

%-------------------end of Section IV--------------------------



















%--------------------Section V---------------------------------------


% generate order of tasks which returns e.g. [3 1 2]
order = randperm(length(tasks));

% intialize cell array that stores the order
taskOrder = cell(1,length(tasks));

% Initialize a cell array to store the timing data
timingDataAll = cell(1, length(tasks));

% try using struct
% timingDataAll = struct();


for i = 1:length(tasks)
    taskIndex = order(i);
    task = tasks{taskIndex};

    if taskIndex == 1
        [timingDataAll{i}, taskNames] = task(sequence, listLo, listHi, listMe, date, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames, 1);
    elseif taskIndex == 2
        [timingDataAll{i}, taskNames] = task(sequence, bodyPartsListV, animalsListV, toolsListV, date, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, taskNames);
    elseif taskIndex == 3
        [timingDataAll{i}, taskNames] = task(sequence, animalsListA, toolsListA, bodyPartsListA, date, ...
    window, white, allCoords, lineWidthPix, xCenter, yCenter, audioDevice, taskNames);
    end
    
end    

%----------------end of Section V-----------------------------------

filename = sprintf('%s_taskOrder.mat', date);

save(filename,'taskNames')

line2_1 = 'This is the end of this task experiment.';
line2_2 = ' Press [space] to exit.';

% Draw instructions         
DrawFormattedText(window, [line2_1 line2_2], 'center', screenYpixels * 0.25, white);

% Flip to the screen
Screen('Flip', window);

% Press any key to continue
KbStrokeWait;

sca;
