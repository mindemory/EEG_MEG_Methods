%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is created by Mrugank Dake (mrugank.dake@nyu.edu)
% And has been adapted for this course.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function screen = initScreen(parameters, devType)
% To make it transparent for working in demo mode
if parameters.isDemoMode
    PsychDebugWindowConfiguration(0, 0.8);
end

screen.id = max(Screen('Screens')); %get the screen
[screen.screenXpixels, screen.screenYpixels] = Screen('WindowSize', screen.id); % get x and y pixels of screen
[screen.screenWidth, screen.screenHeight] = Screen('DisplaySize', screen.id); % get screen width and height in mm
screen.screenWidth = screen.screenWidth/10; % mm to cm
screen.screenHeight = screen.screenHeight/10; % mm to cm
pixWidth = screen.screenWidth/screen.screenXpixels; % cm/pixel
pixHeight = screen.screenHeight/screen.screenYpixels; % cm/pixel
screen.pixSize = mean([pixWidth, pixHeight]); % cm/pixel
screen.xCenter = screen.screenXpixels/2; % pixels
screen.yCenter = screen.screenYpixels/2; % pixels

screen.deg_width = atand(screen.screenWidth/2 / parameters.viewingDistance) * 2; % width of screen in dva
screen.deg_height = atand(screen.screenHeight/2 / parameters.viewingDistance) * 2; % height of screen in dva

% Initialize colors
screen.white = WhiteIndex(screen.id);
screen.black = BlackIndex(screen.id);
screen.grey = screen.white*0.5;

% Initialize Screen Window
AssertOpenGL;
[screen.win, screen.screenRect] = PsychImaging('OpenWindow', screen.id,screen.grey, [], 32, 2, [], [], kPsychNeed32BPCFloat);
screen.ifi = Screen('GetFlipInterval', screen.win);
Screen('BlendFunction', screen.win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Screen('BlendFunction', screen.win, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
if strcmp(devType, 'MEG')
    Screen('TextSize', screen.win, 24);
else
    Screen('TextSize', screen.win, 48);
end
InitializePsychSound(1);

% Retreive the maximum priority number
topPriorityLevel = MaxPriority(screen.win);
Priority(topPriorityLevel);

KbName('UnifyKeyNames');

end