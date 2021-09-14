function [window, pixindeg, diffWandG, grey, xCenter, yCenter, ifi,...
    screenXpixels, screenYpixels] = screen_init(useScreenCalib, fullscreen)
    
%     Screen('Preference','SkipSyncTests', 0)
    % Here we call some default settings for setting up Psychtoolbox
    PsychDefaultSetup(2);

    % Get the screen numbers
    screens = Screen('Screens'); screenNumber = max(screens);

    % Define black and white
    white = WhiteIndex(screenNumber); black = BlackIndex(screenNumber);
    grey = white / 2; diffWandG = abs(white - grey);

    % Load normalized gamma table
    if useScreenCalib
        calib_filename='gammaCalib_12122016.mat';
        load(calib_filename); % this and following code necessary to linearize the monitor was added on 12/12/2016
        Screen('LoadNormalizedGammaTable',screenNumber,CLUT,[]);
    end

    % Open an on screen window
    %use [] for full screen or [0,0,800,400] for a smaller screen
    if fullscreen; screensize=[]; else screensize=[0 0 800 400]; end;
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, screensize);
%     [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, screensize);
    % [window, windowRect] = Screen('OpenWindow',screenNumber, grey, screensize);

    % Get the size of the on screen window in pixels and cm
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
    [width, height]=Screen('DisplaySize', screenNumber);
    ifi = Screen('GetFlipInterval', window);
    [xCenter, yCenter] = RectCenter(windowRect);

    % Get size of pixel in degree. d is distance of observer to screen
    d=57; pixincm=width/10/screenXpixels;
    pixindeg=360/pi*tan(pixincm/(2*d));

    % Set up alpha-blending for smooth (anti-aliased) lines
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA'); 

    topPriorityLevel = MaxPriority(window);
    Priority(topPriorityLevel);
end
