%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Harrison Martin - River Avulsion Simulation                             %
% PhD Candidate, Department of Earth & Atmospheric Sciences               %
% Indiana University                                                      %
% hkmartin@iu.edu                                                         %
% October 2021                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% This code is intented to accompany, and reproduce the results of, a     %
% manuscript submitted to ESurf by myself and my advisor Dr. Douglas A.   %
% Edmonds (edmondsd@iu.edu), who also contributed some code.              %
%                                                                         %
% This release is intended to be as understandable as possible, but as    %
% with any large(r) model, there are some components that might seem      %
% opaque or idiosyncratic (because sometimes they are). I have tried to   %
% comment everything as clearly as possible, and usually used intuitive   %
% variable names, at least for stuff a user would usually interact with.  %
% I also tried to ensure that everything in the code works as advertised  %
% and is stable, but can't warrant perfect performance in all cases. As   %
% such, I removed some old, unusued components or setting options that    %
% are unnecessary or irrelevant to reproducing the manuscript's results,  %
% because I wasn't able to test old components sufficiently to have a     %
% stable release.                                                         %
%                                                                         %
% All said, I'm more than happy to provide whatever help I can if you     %
% want to use or understand this code! Feel free to reach out to me if    %
% you have any questions. I'm also always happy to chat if you have an    %
% interesting research problem that you think the model may be able to    %
% help solve. If you're interested in adding or modifying components,     %
% I'd love to chat! Collaboration is great. Also, if you reach out,       %
% there's a fair chance that I may have already tried it at some point    %
% and have results or code that I would be more than happy to share. I'm  %
% more than happy to collaborate!                                         %
%                                                                         %
% CC BY-NC-ND 4.0 ; Feel free to me with any questions.                   %
%                                                                         %
% Thanks for reading, and I hope you enjoy the code!                      %
%  - Harrison                                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% resetting and experimental variables

% gets rid of old variables and closes old windows, if there are any
clearvars
close all

% Example runs
TitlingForRuns = {'test1','test2'};	% this naming convention has attraction, repulsion, healing direction, and duration ...
% importantly, you have to make a folder with this exact title in exportDirectory
TriggeringForRuns = [30,30];	% characteristic avulsion trigger period (yrs)
AttractionForRuns = [0.25,0.25];	% attraction parameter value
ResistanceForRuns = [0.50,4.00];	% repulsion parameter value
HealingForRuns = [55000,55000];	% characteristic healing timescale (yrs)
InitializationLengthsForRuns = [122500,122500,122500,122500,122500];	% initialization lengths for floodplain slope (m)
HealDirectionForRuns = [3,3];	% choice in healing direction: 0 deposition-only, 1 erosion-only, 3 far-field directed

%% start the model

for runNum = [1:2]	% how many of your pre-set runs do you want to go through
    
    % reset all variables between runs except input experimental variables
    clearvars -except runNum TitlingForRuns TriggeringForRuns HealingForRuns ...
        InitializationLengthsForRuns AttractionForRuns ResistanceForRuns HealDirectionForRuns
    close all
    
    %% codified constants: changing these is anticipated
    tfinal = 1000000;   % total run time (yrs)
    dt = 10;	% timestep (yrs)
    timeSteps = tfinal/dt;  % implied number of timesteps
    numAvulsions = timeSteps/20;    % only used to set limits for colour bars and plots... kind of an initial guess at spectral resolution
    colorSetting = 5;   % 4 is detrended lowelevations, 5 is detrended highelevations
    useCustomColormap = 1;  % I made a custom colourmap, which is parula but with white (for contrast) for the first few lowest values ... 
                            % 1 uses it. 0 uses Matlab default instead (Parula).
    plotAvulsionLoci = 1;   % 0 to not plot avulsion loci, 1 to plot; for Fig.3 at final output
    animateCapturedPathfinding = 0;	% 0 to skip drawing frames for captured pathfinding, 1 otherwise. Skipping vastly speeds up the run.
    animateFloodplainPathfinding = 0;   % 0 to skip drawing frames for floodplain pathfinding, 1 otherwise. Same deal about speed-up.
    animateSubsidence = 1;  % if you have no avulsion on a frame, just subsidence going on, should you still render it in Fig1? 1 yes, 0 no.
    plotInitialEta = 1;	% do we want to update Fig2 during spin-up? Default is every 5,000 timesteps.
    plotOngoingEta = 1; % do we want to update Fig2 throughout the run? Following variable specifies how often to do so.
    plotEvery = 5000;   % every X frames, for plotting Fig2 throughout the run.
    plotSubEvery = 500;     % every X frames; if there's no avulsions going on, how often do we update Fig1?
    blackboxMode = 2;   % how often do we animate frames? speeds up runs. 0 is final output only, 1 is every avulsion, 2 is animateEvery only
    createGIF = 1;  % do you want to export a GIF from the run?
    animateEvery = 2000;    % every X frames; how often do you append a new frame to the GIF?
    toEqualibriumToggle = 0;    % how do we handle diffusion each timestep? 0 does a timestep worth of transient diffusion per timestep.
    insetFirstChannel = 1;  % does the first-ever channel get inset into the background floodplain? technically, we actually raise everything else 1 meandepth
    floodplainInitialLength = InitializationLengthsForRuns(runNum); % m, controls initial floodplain elev
    checkViableEveryYears = TriggeringForRuns(runNum);  % yr, average trigger period
    triggersAreProbabilistic = 1;   % 0 for a trigger occuring deterministically every X years, 1 for 1/X chance per year
    overbankMethod = 0;     % 0 to scale to channel depth
    healedChannelsBecomeFloodplain = 1;   % can abandoned channels eventually stop capturing flow?
    healingThreshold = AttractionForRuns(runNum);   % if yes to above, at which fraction of a mean channel depth in remnant relief?
    channelThreshold = 0.25;    % the fraction of a mean channel depth in remnant relief at which you turn abandoned channel cells into floodplain cells
    requireElevAdvantage = 1;   % do you require an elevation advantage (~gradient advantage) in addition to superelevation?
    useFactorForSE = 1;     % should you use a SE factor in units of mean depth (~beta - 1; instead of a fixed value in metres)
    SEFactor = 0.00;    % what fraction of a mean channel depth must an active channel bed rise above a neighboring low to be superelevated (~beta - 1)
    incisionRule = 1;   % when a new channel is created, it incises one channel depth down from 1) that cell's highelev
    exportOutput = 1;   % do you want to save screenshots of Fig1 every so often?
    exportEvery = 2000; % every X frames; if yes to above, how often?
    timeRuns = 0;   % if 1, times 5kyrs of simulation and then stops it, with a projection of how long 10Myr would take at that rate.
    trackLobeSwitching = 1;	% do you want to track channel position for every step? Fig8. Creates very large, slow .csv files.
    saveAvulsionLoci = 1;   % do you want to save where every avulsion occured? creates a pretty large .csv file
    diflIgnoresAlluvial = 3;    % overbank deposition applied to abandoned channel cells: 3 is same as row's neutral floodplain amount, ...
                                % (which causes no healing due to overbank deposition)
    hardcodeHealing = 1;    % do you want abandoned channels to heal (aside from overbank deposition)
    hardcodeHealingDirection = HealDirectionForRuns(runNum);	% 0 for 1-ended low-towards-high (infill only), 1 for 1-ended high-towards-low (erosion only), ...
                                                                % 3 for both going toward neutral floodplain
    hardcodeYearsToHeal = HealingForRuns(runNum);	% characteristic channel healing timescale, years, in multiples of timesteps, please
    disallowInheritance = 0;	% 0 checks relative superelevation, 1 requires a strict amount of aggradation per cell, 2 is the same as 1 but resets ... 
                                % aggradation the moment a cell isn't still occupied
    constrainDifl = 1;	% will you prevent overbank deposition from exceeding the magnitude of subsidence for a cell at each timestep?
    randomWalkMode = 1;	% 0 ignores slope, just using pre-determined weights; 1 is slope-ordered weights
    numSteps = 1000000;	% maximum timesteps per avulsion cycle before the simulation aborts; for troubleshooting, shouldn't matter as long as it's big

    
    % where does output get saved to?
%     exportDirectory = '/N/u/hkmartin/Carbonate/Documents/Revisions/Misc Runs/SwitchingTests/5Myrs';
    exportDirectory = 'C:\Users\Owner\Desktop\Programming\Matlab\Research\Misc Tests';
    exportFolder = TitlingForRuns{runNum};
    exportName = TitlingForRuns{runNum};
    GIFfilename = exportName;
    
    initialElevation = 1;	% do you want an initial floodplain elevation before initializing? colourbar settings assumes this value is >= numAvulsions
    
    gridWidth = 301;	% domain space size, in cells, horizontal/columns
    gridLength = 301;	% domain spaze size, in cells, vertical/rows
    cellSize = 500;	% m, dimension per side of the square cell
    
    % for randomWalkMode == 0; odds to move each direction. Enter as percentages (0 to 100, not decimals)
    % oddsUp = 0; % should remain zero; unexpected behaviour otherwise
    oddsR = 10;
    oddsDR = 20;
    oddsD = 40;
    oddsDL = 20;
    oddsL = 10;	% assumes that all five add up to 100!
    
    % the top and bottom betaSub values. Linearly interpolates to fill betaSubMtx.
    topBetaSub = 0.000010;	% how many metres to subside per year at the top/proximal end of the domain
    botBetaSub = 0.000005;	% same, but at bottom/distal end of the domain
    subsidence = topBetaSub;	% subsidence in meters per year (used in the 1D diffusion code)
    
    % in metres
    minSuperelevation = 0.625;	% m, approximately 0.5 mean-channel-depths; only used for first timestep and if useFactorForSE == 0
    
    % repulsion factor; how many times larger than flow depth does levee height have to be to repel an approaching avulsion
    resistanceFactor = ResistanceForRuns(runNum);
    
    % uniform overbank deposition constants. Linear interpolation between 'em.
    topDeposition = 0.0000001;	% m/yr, proximal
    botDeposition = 0.0000025;	% m/yr, distal
    
    % differential overbank deposition constants. Unit values, which are scaled based on cell elevation difference from max elevation of each row
    topDiffDeposition = 0.0000001;
    botDiffDeposition = 0.0000025;
    
    seCbarLow = 0;	% colourbar bounds for superelevation views, in metres
    seCbarHigh = 20;
    
    % 1D diffusion code vars (Paola et al. 1992):
    qo = 190000;    % initial specific discharge m2/yr at head of stream
    basinwidth = 50000;	% m, used for calculating discharge
    Acoeff = 1;
    cf = 0.01;	% nondimensional coefficient of friction = f/8
    C0 = 0.7;
    rohS = 2650;	% kg/m^3
    qsin = 400;	% incoming sediment supply, m^3/yr
    
    % stuff used to solve for depth
    rohW = 1000;	% kg/m^3, roh sub w, density of water
    g = 9.81;	% m^2/s, accelleration due to gravity 
    % via reference table at http://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm
    D50 = 0.05;	% m, grain size
    shields = 0.045;	% shields' parameter, via https://pubs.usgs.gov/sir/2008/5093/table7.html
    
    %% derivative variables (for calculations)
    
    centerCol = ceil(gridWidth/2); % finds the horizontal midpoint; used for initial starting point
    
    % errorchecking to ensure all four directional probabilities add up to 100%
    if oddsR + oddsDR + oddsD + oddsDL + oddsL ~= 100
        error ('Odds for all four movement directions must sum to 100%!');
    end
    
    % converting from percentage to decimal
    oddsR = oddsR/100;
    oddsDR = oddsDR/100;
    oddsD = oddsD/100;
    oddsDL = oddsDL/100;
    oddsL = oddsL/100;
    
    primeOddsR = oddsR; % for random walk mode == 0
    primeOddsDR = oddsDR;
    primeOddsD = oddsD;
    primeOddsDL = oddsDL;
    primeOddsL = oddsL;
    
    capturedOddsL = (oddsL) / (oddsL + oddsR); % for when you can only move L or R, scales the new odds to 0 - 1
    capturedOddsR = 1;
    
    capturedOddsDL = (oddsDL) / (oddsDL + oddsDR); % for when you can only move DL or DR, scales the new odds to 0 - 1
    capturedOddsDR = 1;
    
    topBetaSubPerStep = topBetaSub; % vestigial; we apply dt at the time of subsidence, so this can be a simple equality
    botBetaSubPerStep = botBetaSub;
    
    primeFirstStepOddsR = oddsR / (oddsR + oddsD + oddsL); % for random walk mode == 0, for first steps
    primeFirstStepOddsD = primeFirstStepOddsR + (oddsD / (oddsR + oddsD + oddsL));
    primeFirstStepOddsL = primeFirstStepOddsD + (oddsL / (oddsR + oddsD + oddsL));
    
    FSOwalkWeights = [0.50 0.35 0.15]; % for random walk mode == 1, first step weights
    walkWeights = [0.40 0.275 0.175 0.10 0.05]; % for rando mwalk mode == 1, five possible weights
    
    checkViableEvery = checkViableEveryYears / dt; % from years to steps
    hardcodeTimestepsToHeal = hardcodeYearsToHeal / dt; % from years to steps
    
    if randomWalkMode == 0 % just handles how we assign first-step pathfinding weights if we don't slope-weight the random walk
        
        firstStepOddsR = primeFirstStepOddsR;
        firstStepOddsD = primeFirstStepOddsD;
        firstStepOddsL = primeFirstStepOddsL;
    end
    
    %% initializing variables
    
    gridspace = zeros(gridLength,gridWidth); % creates a 2D matrix to work in; cell value is the cGen of the last channel to visit the cell
    % key: cell value 0 = never visited or fully healed to floodplain, j = visited during generation j
    
    % sets up a matrix to store active channel path history, used for finding a new avulsion starting point
    pathHist = zeros(4*gridLength,2);
    
    % keeps track of upstream vs downstream
    pathStep = zeros(gridLength,gridWidth);
    stepCounter = 1;
    
    % tracks elevation of each cell, colourbar assumes minimum is zero.
    lowelevation = zeros(gridLength,gridWidth) + initialElevation;
    % same as low, but with levee heights wherever a river is or has been
    % highelevation = elevation; % <- unusued, since it's copied over from 'elevation' later
    
    % tracks elevation of each cell relative to the neutral floodplain.
    superelevation = zeros(gridLength,gridWidth);
    leveesuperelevation = zeros(gridLength,gridWidth);
    
    % tracks 'background' elevation for the mtx.
    fpelevmtx = zeros(gridLength,gridWidth);
    
    % creates a matrix of beta subsidence values for each cell
    betaSubMtx = zeros(gridLength,gridWidth);
    slopeBSM = (topBetaSubPerStep - botBetaSubPerStep) / (gridLength-1);
    for gsr = 1:gridLength
        betaSubMtx(gsr,:) = botBetaSubPerStep + (slopeBSM * (gsr-1));
    end
    
    % creates a matrix of Uniform overbank deposition values for each cell
    depositionMtx = zeros(gridLength,gridWidth);
    slopeOBD = (topDeposition - botDeposition) / gridLength;
    for obdm = 1:gridLength
        depositionMtx(obdm,:) = botDeposition + (slopeOBD * obdm);
    end
    
    % creates a matrix of Differential overbank deposition values for each cell
    diffDepositionPreMtx = zeros(gridLength,gridWidth);
    slopeODBD = (topDiffDeposition - botDiffDeposition) / gridLength;
    for odb = 1:gridLength
        diffDepositionPreMtx(odb,:) = botDiffDeposition + (slopeODBD * odb);
    end
    
    % used for calculating and applying differential deposition
    rowminelev = zeros(gridLength,1);
    rowmaxelev = zeros(gridLength,1);
    lowdiffdeposcalar = zeros(gridLength,gridWidth); % same, but normalized to have each row sum to gridWidth, for low elevs
    highdiffdeposcalar = zeros(gridLength,gridWidth); % " for high (aka levee) elevations
    neutralElevMtx = zeros(gridLength,gridWidth); % a matrix with farfield elevations projected along-strike across the domain
    
    % needs to be one for pathHist to start at row 1 for the initial runthrough
    newRandEntry = 1;
    currentEntry = 1;
    
    % tracks, for each pathfinding step, whether the channel is currently captured by an abandoned channel
    captured = 0; % 0 = no, 1 = yes
    
    % tells the first-step-of-avulsion pathfinding when to break out of a loop
    newPathFound = 0;
    
    % a table that stores avulsion locus history/details
    avulsionLoci = zeros(timeSteps,5);
    
    % a list of all cells for which avulsions can occur
    viableLoci = zeros(4*gridLength,3);
    
    % did an avulsion attempt succeed?
    aTryTracker = zeros(timeSteps,1);
    
    % tracks what avulsion generation we're on
    cGen = 1;
    
    % tracks what time steps experienced avulsions
    didAvulse = zeros(timeSteps,1);
    
    fig1 = figure(1); % creates the figure that shows the basin and river
    set(fig1,'Units','normalized');
    set(fig1,'Position',[0 0 0.50 0.90]); % configured for my desktop monitor setup, but seems fairly flexible on other monitors
    
    % variables used in calculating depth
    depth = zeros(2*gridLength,1); % depth, m
    velocity = zeros(2*gridLength,1); % velocity, m/s
    To = zeros(2*gridLength,1); % Tau naught, bed-surface shear stress
    slope = zeros(2*gridLength,1); % slope, as a decimal/fraction
    chanwidth = zeros(2*gridLength,1); % channel width, in m. not used for much other than curiosity
    
    % variables used in recording OG (meaning initial runthrough) depth
    OGdepth = zeros(2*gridLength,1); % depth, m
    OGvelocity = zeros(2*gridLength,1); % velocity, m/s
    OGTo = zeros(2*gridLength,1); % Tau naught, bed-surface shear stress
    OGchanwidth = zeros(2*gridLength,1); % channel width, in m. not used for much other than curiosity
    
    stepsThisGen = 1; % counts the number of pathfinding steps this generation; can't be tied to i since some steps can fail
    
    stepDirections = zeros(2*gridLength,2); % tracks the direction each step took, used for total length; R/D/L are col1, DR/DL are col2
    
    orthoSize = cellSize; % metres per orthogonal step
    diagSize = sqrt(2) * cellSize; % metres per diagonal step
    
    Ltracker = zeros(timeSteps,1); % tracks river length at each step
    
    progradLength = zeros(timeSteps,3); % [deprecated] tracks how long each avulsion goes for before being captured or exiting
    progradLength(1,2) = 1; % step num
    progradLength(1,3) = 100; % row num
    
    triggerCount = 0; % how many triggers have occurred over the whole simulation?
    
    if trackLobeSwitching == 1
        lobeTracker = zeros(timeSteps,9); % this matrix will track lobe switching at up to 8 distances from the mountain-front
    end
    
    subOnlySEtracker = zeros(gridLength,1); % you want to know how superelevated you'd be from subsidence alone, along each row? this it is
    SEtracker = zeros(gridLength,1); % total superelevation experienced by each row
    aggradTracker = zeros(gridLength,1); % total aggradation experienced by each row
    repulsionTracker = zeros(gridLength,1);% how many repulsions occurred this step
    
    if useCustomColormap == 1 % this is where my custom colour map shows up
        load whitetip4 % it's parula but with the first few bars (I think 0.25 channel depths?) whited out, to make it easier to see
        customParula = whitetip4;
    else
        customParula = parula;
    end
    
    SETmtx = zeros(gridLength,gridWidth); % track superelevation experienced by each cell
    preSETmtx = zeros(gridLength,gridWidth); % stores the "before", so we can calculate absolute superelevation
    
    isChannel = zeros(gridLength,gridWidth); % is this cell an abandoned channel or is it floodplain?
    
    viableLociTracker = zeros(gridLength,gridWidth); % tracks how many times each cell has been a viable avulsion location
    
    %% MODEL BEGINS
    
    % Floodplain elevation initialization:
    
    % run the first 1D model
    % Set the number of grid points and build a cell-center grid
    N = gridLength-2; % number of cells in the vector
    L = floodplainInitialLength; % length of the domain, in cells
    dx = L/N;
    x = -.5*dx:dx:L+.5*dx; % the 1D set of spaced gridpoints in the 'x' direction
    x = x'; % Turn x into a column vector.
    
    if subsidence > 0
        sigma = (subsidence:-subsidence/(2*(N+1)):(subsidence/2))'; % takes subsidence at upstream boundary, interpolates to half-value at bottom
    else % catches for div by 0
        sigma = zeros(N+2,1);
    end
    
    q = zeros(length(x),1)+qo; % specific discharge of water, m^2/yr, i.e. rate*basin length. For L=100,000 m, and rain = 1 m/yr, q = 100,000 m2/yr
    D = (8.*q.*Acoeff.*sqrt(cf))/(C0*(2.65-1)); % fluvial diffusivity, m^2/yr
    
    % Load Dm with average values D(j-1/2) and Dp with D(j+1/2)
    Dm = zeros(N+2,1);
    Dp = zeros(N+2,1); % Make the column vectors
    Dm(2:N+1) = .5*(D(2:N+1)+D(1:N)); % average j and j-1
    Dp(2:N+1) = .5*(D(2:N+1)+D(3:N+2)); % average j and j+1
    
    % initialize eta
    eta = zeros(N+2,1);
    
    % time variables
    fdt = 10; % timestep (yrs)
    ftfinal = 1000000; % total runtime (years); make it long enough to reach eq'm
    fnsteps = ftfinal/fdt; % implied #steps
    
    % sediment variables
    const = 2*dx^2 / fdt;
    eq_slope = qsin/D(1);
    
    % Create the sparse matrices A and B
    aneg1 = [-Dm(2:N+1); 0;0]; % extra zeros are needed to pad row N+2 for boundary conditions
    a0 = [0;const + (Dm(2:N+1)+Dp(2:N+1));0]; % extra zeros are needed to pad row N+2 for boundary conditions
    a1 = [0;0; -Dp(2:N+1)];% extra zeros are needed to pad row N+2 for boundary conditions
    A = spdiags([aneg1 a0 a1],[-1 0 1],N+2,N+2);
    
    bneg1 = [Dm(2:N+1); 0;0]; % extra zeros are needed to pad row N+2 for boundary conditions
    b0 = [0;const - (Dm(2:N+1)+Dp(2:N+1));0]; % extra zeros are needed to pad row N+2 for boundary conditions
    b1 = [0;0; Dp(2:N+1)];% extra zeros are needed to pad row N+2 for boundary conditions
    B = spdiags([bneg1 b0 b1],[-1 0 1],N+2,N+2);
    
    % load the boundary conditions into A and B A(1,1)=1 and A(1,2)=1 makes first equation eta1-eta2=r(1), ...
    % this means that r(1) is the elevation difference between the cells which we set as the equilibrium slope ...
    % needed to transport sediment supplied.
    A(1,1) = 1;
    A(1,2) = -1;
    B(1,1) = 0;
    % T(0)=0
    A(N+2,N+1) = 0.5;
    A(N+2,N+2) = 0.5;
    B(N+2,N+2) = 0;
    % T(L)=0
    
    % This is the time advance loop for diffusion
    for mtime = 1:fnsteps
        
        % find the right-hand side for the solution at interior points
        r = B*(eta-sigma.*fdt);
        % apply the boundary conditions
        r(1) = eq_slope*dx; % T(0)=0
        r(N+2) = 0;
        % do the linear solve to update T
        eta = A\r;
        
        % for debugging, save eta0
        eta0 = eta;
        x0 = x;
        
        % Make a plot of T every once in a while (every 50k years).
        if plotInitialEta == 1 && (rem(mtime,5000) == 0)
            figure(2)
            plot(x,eta)
            %hold on
            title([num2str(mtime*fdt/ftfinal*100) ' percent complete'])
            pause(.01)
        end
        
    end
    
    % implementing new elevation along the floodplain
    for fpe = 1:gridLength
        lowelevation(fpe,:) = eta(gridLength+1-fpe);
    end
    
    % copying elevation to highelevation
    highelevation = lowelevation;
    
    % preserving ghost node elevation for eta(1)
    ghostElev = eta(1);
    ghostSigma = sigma(1);
    ghostSigma2 = sigma(N+2);
    
    % for Fig3, we track what's happening to entry and exit elevs over time
    elevTracker = zeros(timeSteps,4);
    
    % for bugfixing, checks what the topleft corner does for the whole run
    cornerTracker = zeros(timeSteps,1);
    
    % Determine channel slope
    for slp = 1:(N-1)
        slope(slp) = abs((eta(slp+1) - eta(slp)) / dx);
    end
    
    % sets the slope for the final (unused) cell equal to the penultimate one
    slope(N) = slope(N-1);
    
    % define a min & max slope based on the range of values observed in the first, unperturbed runthru
    maxslope = 2.0*slope(1);
    minslope = 0.5*slope(N);
    
    % saving OG parameters
    OGQ = zeros(length(eta),1)+(qo*basinwidth/3.154e7);
    OGslope = slope;
    
    for OGdpt = 1:length(eta)-2
        
        % METHOD 4: solving via shields
        OGTo(OGdpt) = shields * ((rohS - rohW) * g * D50);
        
        % calculate velocity via empirical bed shear stress eq'n
        OGvelocity(OGdpt) =  sqrt(OGTo(OGdpt)/(cf*rohW));
        
        % calculate depth via open channel analytical shear stress eq'n
        OGdepth(OGdpt) = OGTo(OGdpt) / (rohW * g * OGslope(OGdpt));
        
        % calculate depth via definition of discharge
        OGchanwidth(OGdpt) = OGQ(OGdpt)/(OGdepth(OGdpt)*OGvelocity(OGdpt));
        
    end
    
    %% MODEL BEGINS (for real)
    
    for j = 1:timeSteps %% Avulsion loop! Timesteps
        
        if j == 1 % gets an initial reading of elevations
            % tracks point source height
            elevTracker(j,1) = eta(1);
            elevTracker(j,2) = eta(2);
            % tracks output height
            elevTracker(j,3) = eta(N+1);
            elevTracker(j,4) = eta(N+2);
            % tracks corner
            cornerTracker(j,1) = lowelevation(gridLength,1);
        end
        
        if timeRuns == 1 % if we want to try a run, it does a tiny one and extrapolates to guess how long a big one would take
            
            if j == 100
                tic
            end
            
            if j == 600
                timedRun = toc;
                sprintf(['To run 5k years, it took ' num2str(timedRun) ' seconds.'])
                sprintf('That implies that a full 10Myr run would take an estimated:')
                sprintf([num2str(timedRun*2000/60) ' minutes...'])
                sprintf(['or ' num2str(timedRun*2000/3600) ' hours.'])
                error('Finished timing!') % stops the run so you can see
            end
            
        end
        
        if trackLobeSwitching == 1
            % tracks where the active channel is at different lengths downstream
            lobeTracker(j,1) = j*dt;
            %25th row
            for ltk = 1:gridWidth
                if gridspace(25,ltk) == cGen
                    lobeTracker(j,2) = ltk;
                end
            end
            %50th row
            for ltk = 1:gridWidth
                if gridspace(50,ltk) == cGen
                    lobeTracker(j,3) = ltk;
                end
            end
            %100th row
            for ltk = 1:gridWidth
                if gridspace(100,ltk) == cGen
                    lobeTracker(j,4) = ltk;
                end
            end
            %150th row
            for ltk = 1:gridWidth
                if gridspace(150,ltk) == cGen
                    lobeTracker(j,5) = ltk;
                end
            end
            %200th row
            for ltk = 1:gridWidth
                if gridspace(200,ltk) == cGen
                    lobeTracker(j,6) = ltk;
                end
            end
            %225th row
            for ltk = 1:gridWidth
                if gridspace(225,ltk) == cGen
                    lobeTracker(j,7) = ltk;
                end
            end
            %250th row
            for ltk = 1:gridWidth
                if gridspace(250,ltk) == cGen
                    lobeTracker(j,8) = ltk;
                end
            end
            %275th row
            for ltk = 1:gridWidth
                if gridspace(275,ltk) == cGen
                    lobeTracker(j,9) = ltk;
                end
            end
        end
        
        if j ~= 1 % don't subside or diffuse sediment on the first go, since you don't have a channel path to restore or update
            
            oldEta = eta; % for the "before" picture for debugging or calculating changes
            
            if disallowInheritance == 1 || disallowInheritance == 2 % for absolute aggradation for superelevation
                preSETmtx(gridspace == cGen) = lowelevation(gridspace == cGen);
            end
            
            if toEqualibriumToggle == 0
                % find the right-hand side for the solution at interior points
                r = B*(eta-sigma.*dt);
                % apply the boundary conditions
                r(1) = eq_slope*dx; % T(0)=0
                r(N+2) = 0;
                % do the linear solve to update T
                eta = A\r;
            else
                error('Invalid toEquilibriumToggle selection. Has to be 0, champ.')
            end
            
            % tracks point source height
            elevTracker(j,1) = eta(1);
            elevTracker(j,2) = eta(2);
            elevTracker(j,3) = eta(N+1);
            elevTracker(j,4) = eta(N+2);
            cornerTracker(j,1) = lowelevation(gridLength,1);
            
            % resets the aggradation tracker for this step
            aggradThisStep = zeros(gridLength,5);
            
            % loops along path to check how much aggradation
            for ats = 1:sum(pathHist(:,1)~=0)
                
                % tracks aggradation per row
                ATScol = 3;
                % will overestimate if there are multiple actives in a row but one had zero change
                if aggradThisStep(pathHist(ats,1),2) == 0 % odds of 'zero change' are expected to be low enough to not matter
                    aggradThisStep(pathHist(ats,1),2) = eta(ats) - oldEta(ats);
                else
                    aggradThisStep(pathHist(ats,1),ATScol) = eta(ats) - oldEta(ats);
                    ATScol = ATScol + 1;
                end
                
            end
            
            % corrects for aggradation per row when there's multiple cells per row part of the active channel (by averaging)
            sATS = size(aggradThisStep);
            for ata = 2:gridLength-1
                
                aggradThisStep(ata,1) = mean(aggradThisStep(ata,2:1+sum(aggradThisStep(ata,2:sATS(2))~=0)));
                
            end
            
            aggradTracker = aggradTracker + aggradThisStep(:,1);
            
            % once again takes a "before" snapshot of the active channel elevations
            eta0 = eta;
            
            % insets the first channel, depending on your setting
            if j == 2 && insetFirstChannel == 1
                
                meandepth = mean(depth(1:(sum(depth(:,1)~=0))));
                lowelevation(3:gridLength,:) = lowelevation(3:gridLength,:) + meandepth;
                highelevation(3:gridLength,:) = highelevation(3:gridLength,:) + meandepth;
                
            end
            
            % implementing new elevation along the active channel
            for ne = 1:sum(pathHist(:,1)~=0)
                lowelevation(pathHist(ne,1),pathHist(ne,2)) = eta(ne+1);
            end
            
            % adds levee height to highelevation mtx
            for hel = 1:sum(pathHist(:,1)~=0)
                if lowelevation(pathHist(hel,1),pathHist(hel,2)) + depth(hel) > highelevation(pathHist(hel,1),pathHist(hel,2))
                    highelevation(pathHist(hel,1),pathHist(hel,2)) = lowelevation(pathHist(hel,1),pathHist(hel,2)) + depth(hel);
                end
            end
            
            % if you don't allow relative superelevation, then calculate how much superelevation you earned per-cell this timestep
            if disallowInheritance == 1 || disallowInheritance == 2
                SETmtx(gridspace == cGen) = SETmtx(gridspace == cGen) + (lowelevation(gridspace == cGen) - preSETmtx(gridspace == cGen));
            end
            
            % if you selected this, deletes progress from abandoned channels
            if disallowInheritance == 2
                SETmtx(gridspace ~= cGen) = 0;
            end
            
            % checks farfield lowelevations
            preSEfloodplain = lowelevation(:,1);
            
            % SUBSIDENCE
            % elevation is lowered for whole gridspace
            lowelevation = lowelevation - dt*betaSubMtx;
            
            % elevation is restored for active channel (lowering is exactly undone on a per-cell basis)
            for ei = 1:(sum(pathHist(:,1)~=0))
                lowelevation(pathHist(ei,1),pathHist(ei,2)) = lowelevation(pathHist(ei,1),pathHist(ei,2)) + dt*betaSubMtx(pathHist(ei,1),pathHist(ei,2));
            end
            
            % applies subsidence to highelevation mtx, too
            highelevation = highelevation - dt*betaSubMtx;
            
            % (high)elevation is restored for active channel (lowering is exactly undone on a per-cell basis)
            for ei2 = 1:(sum(pathHist(:,1)~=0))
                highelevation(pathHist(ei2,1),pathHist(ei2,2)) = highelevation(pathHist(ei2,1),pathHist(ei2,2)) + dt*betaSubMtx(pathHist(ei2,1),pathHist(ei2,2));
            end
            
            % OVERBANK DEPOSITION
            
            % calculate differential/proportional deposition
            
            % calculate min and max for each row
            for rme = 1:gridLength
                rowminelev(rme) = min(lowelevation(rme,:));
                rowmaxelev(rme) = max(highelevation(rme,:));
            end
            
            % calculate mean channel depth and superelevation factor (to handle diff'l depos'n)
            meandepth = mean(depth(1:(sum(depth(:,1)~=0))));
            if useFactorForSE == 1
                minSuperelevation = meandepth*SEFactor;
            end
            
            if overbankMethod == 0 % unit-depth-scaling method
                
                if j ~= 1
                    % reinitialize deposcalars
                    lowdiffdeposcalar = zeros(gridLength,gridWidth);
                    highdiffdeposcalar = zeros(gridLength,gridWidth);
                end
                
                if diflIgnoresAlluvial == 3 % same as neutralelev cells for row                    
                    % create a matrix of 'neutral' elevations based on undeformed floodplain elevs
                    for nemx = 1:gridLength
                        neutralElevMtx(nemx,:) = lowelevation(nemx,1);
                    end                   
                    % create a new deposition matrix that accounts for elevation difference from
                    % max(row), for low elevations
                    for dmx = 1:gridLength                       
                        lowdiffdeposcalar(dmx,1) = ((rowmaxelev(dmx) - lowelevation(dmx,1))./meandepth);
                        lowdiffdeposcalar(dmx,:) = lowdiffdeposcalar(dmx,1); % if alluvial                        
                    end                    
                    % create a new deposition matrix that accounts for elevation difference from
                    % max(row), for high elevations
                    for dmx = 1:gridLength                        
                        highdiffdeposcalar(dmx,1) = ((rowmaxelev(dmx) - highelevation(dmx,1))./meandepth);
                        highdiffdeposcalar(dmx,:) = highdiffdeposcalar(dmx,1); % if alluvial                        
                    end                    
                else                    
                    error('invalid diflIgnoresAlluvial selection. Only 3''s allowed.')   
                end
            
            end
            
            % creates two matrices, by multiplying the diffdepositionPreMtx by the
            % two depositional scalars, one for low and one for high
            lowDiffDepoMtx = lowdiffdeposcalar.*diffDepositionPreMtx;
            highDiffDepoMtx = highdiffdeposcalar.*diffDepositionPreMtx;
            
            if constrainDifl == 1 % prevents aggradation from exceeding subsidence
                lowDiffDepoMtx(lowDiffDepoMtx>betaSubMtx) = betaSubMtx(lowDiffDepoMtx>betaSubMtx);
                highDiffDepoMtx(highDiffDepoMtx>betaSubMtx) = betaSubMtx(highDiffDepoMtx>betaSubMtx);
            end
            
            % elevation is adjusted by uniform deposition across the whole gridspace
            lowelevation = lowelevation + (dt*depositionMtx);
            
            % applies uniform deposition to highelevation mtx, too
            highelevation = highelevation + (dt*depositionMtx);
            
            % low elevation is adjusted by differential deposition
            lowelevation = lowelevation + (dt*lowDiffDepoMtx);
            
            % high elevation is adjusted by differential deposition
            highelevation = highelevation + (dt*highDiffDepoMtx);
            
            % lowelevation global deposition is undone along the active channel (deposition is exactly undone on a per-cell basis)
            for di = 1:(sum(pathHist(:,1)~=0))
                lowelevation(pathHist(di,1),pathHist(di,2)) = lowelevation(pathHist(di,1),pathHist(di,2)) - dt*depositionMtx(pathHist(di,1),pathHist(di,2));
            end
            
            % highelevation global deposition is undone along the active channel (deposition is exactly undone on a per-cell basis)
            for di2 = 1:(sum(pathHist(:,1)~=0))
                highelevation(pathHist(di2,1),pathHist(di2,2)) = highelevation(pathHist(di2,1),pathHist(di2,2)) - dt*depositionMtx(pathHist(di2,1),pathHist(di2,2));
            end
            
            % lowelevation differential deposition is undone along the active channel (deposition is exactly undone on a per-cell basis)
            for di3 = 1:(sum(pathHist(:,1)~=0))
                lowelevation(pathHist(di3,1),pathHist(di3,2)) = lowelevation(pathHist(di3,1),pathHist(di3,2)) - dt*lowDiffDepoMtx(pathHist(di3,1),pathHist(di3,2));
            end
            
            % highelevation differential deposition is undone along the active channel (deposition is exactly undone on a per-cell basis)
            for di4 = 1:(sum(pathHist(:,1)~=0))
                highelevation(pathHist(di4,1),pathHist(di4,2)) = highelevation(pathHist(di4,1),pathHist(di4,2)) - dt*highDiffDepoMtx(pathHist(di4,1),pathHist(di4,2));
            end
            
            % tracks what superelevation would have been from subsidence alone
            subOnlySEtracker = subOnlySEtracker + (preSEfloodplain - lowelevation(:,1));
            % tracks how much superelevation from overbank deposition occurred
            SEtracker = SEtracker + (preSEfloodplain - lowelevation(:,1)) + aggradThisStep(:,1);
            
            if hardcodeHealing == 1 % if you want to apply healing rules
                
                if hardcodeHealingDirection == 0 % one-ended low-to-high; deposition-only
                    
                    % converting from time to depth; how much to heal per timestep
                    HChealingAmount = meandepth / hardcodeTimestepsToHeal;
                    
                    % figure out which cells should be modified
                    HCnotCurrent = (gridspace~=cGen);
                    HCnotHealed = (lowelevation~=highelevation);
                    HCtoBeHealed = and(HCnotCurrent,HCnotHealed);
                    
                    % raises lows
                    lowelevation(logical(HCtoBeHealed)) = lowelevation(logical(HCtoBeHealed)) + HChealingAmount;
                    % corrects for overfill
                    HCtooHigh = and(lowelevation>highelevation,HCtoBeHealed);
                    lowelevation(HCtooHigh) = highelevation(HCtooHigh);
                    
                elseif hardcodeHealingDirection == 1 % one-ended high-to-low; erosion-only
                    
                    % converting from time to depth; how much to heal per timestep
                    HChealingAmount = meandepth / hardcodeTimestepsToHeal;
                    
                    % figure out which cells should be modified
                    HCnotCurrent = (gridspace~=cGen);
                    HCnotHealed = (lowelevation~=highelevation);
                    HCtoBeHealed = and(HCnotCurrent,HCnotHealed);
                                       
                    % lowers highs
                    highelevation(logical(HCtoBeHealed)) = highelevation(logical(HCtoBeHealed)) - HChealingAmount;
                    % corrects for overfill
                    HCtooLow = and(highelevation<lowelevation,HCtoBeHealed);
                    highelevation(HCtooLow) = lowelevation(HCtooLow);
                    
                elseif hardcodeHealingDirection == 3 % both ends trend towards neutralelev; farfield-directed
                    
                    % create a matrix of 'neutral' elevations based on undeformed floodplain elevs
                    for nemx = 1:gridLength
                        neutralElevMtx(nemx,:) = lowelevation(nemx,1);
                    end
                    
                    % converting from time to depth; how much to heal per timestep
                    HChealingAmount = meandepth / hardcodeTimestepsToHeal;
                    HChealingAmountPerEnd = HChealingAmount/2;
                    
                    % figure out which cells should be modified
                    HCnotCurrent = (gridspace~=cGen);
                    HClowElevtobeRaised = and(HCnotCurrent,lowelevation<neutralElevMtx);
                    HChighElevtobeRaised = and(HCnotCurrent,highelevation<neutralElevMtx);
                    HClowElevtobeLowered = and(HCnotCurrent,lowelevation>neutralElevMtx);
                    HChighElevtobeLowered = and(HCnotCurrent,highelevation>neutralElevMtx);
                    
                    % raise lowelevs
                    lowelevation(logical(HClowElevtobeRaised)) = lowelevation(logical(HClowElevtobeRaised)) + HChealingAmountPerEnd;
                    % check/correct for overfill
                    HClowElevTooRaised = and((lowelevation>neutralElevMtx),HClowElevtobeRaised);
                    lowelevation(logical(HClowElevTooRaised)) = neutralElevMtx(logical(HClowElevTooRaised));
                    
                    % raise highelevs
                    highelevation(logical(HChighElevtobeRaised)) = highelevation(logical(HChighElevtobeRaised)) + HChealingAmountPerEnd;
                    % check/correct for overfill
                    HChighElevTooRaised = and((highelevation>neutralElevMtx),HChighElevtobeRaised);
                    highelevation(logical(HChighElevTooRaised)) = neutralElevMtx(logical(HChighElevTooRaised));
                    
                    % lower lowelevs
                    lowelevation(logical(HClowElevtobeLowered)) = lowelevation(logical(HClowElevtobeLowered)) - HChealingAmountPerEnd;
                    % check/correct for lowering too much
                    HClowElevTooLowered = and(lowelevation<neutralElevMtx,HClowElevtobeLowered);
                    lowelevation(logical(HClowElevTooLowered)) = neutralElevMtx(logical(HClowElevTooLowered));
                    
                    % lower highelevs
                    highelevation(logical(HChighElevtobeLowered)) = highelevation(logical(HChighElevtobeLowered)) - HChealingAmountPerEnd;
                    % check/correct for lowering too much
                    HChighElevTooLowered = and(highelevation<neutralElevMtx,HChighElevtobeLowered);
                    highelevation(logical(HChighElevTooLowered)) = neutralElevMtx(logical(HChighElevTooLowered));
                    
                    % there shouldn't be a case, because of my implementation,where
                    % lowelev ends up above highelev. at least to my understanding.
                    % if this happens, I wanna flag it as an error so I know and
                    % can understand what cases result in this.
                    if any(any(lowelevation>highelevation))
                        error('Lowelev has some cell that''s above highelev')
                    end                    
                else                    
                    error('invalid hardcode healing direction selection, champ. 0, 1, or 3.')                    
                end               
            end            
        end
        
        % updates the matrix that tracks non-floodplain cells
        isChannel(gridspace~=0) = 1;
        
        if healedChannelsBecomeFloodplain == 1 && j~=1 % no need to use this double-loop scan if you aren't removing healed paleochannels
                                              % just kidding, I learned how matrix math works and recoded this
            gridspace(and(and((gridspace~=0),(gridspace~= cGen)),((highelevation-lowelevation) < (healingThreshold*meandepth)))) = 0;
            isChannel(and(gridspace~=cGen,(highelevation-lowelevation) < (channelThreshold*meandepth))) = 0; 
        end
        
        if colorSetting == 4 && (blackboxMode == 1 || blackboxMode == 2) % only spend the memory/time to make and calculate these matrices if you're gonna use 'em.
            for fpemr = 1:gridLength
                fpelevmtx(fpemr,:) = lowelevation(fpemr,1); % detrend lowelevations
            end
            superelevation = lowelevation - fpelevmtx;
        elseif colorSetting == 5 && (blackboxMode == 1 || blackboxMode == 2)
            for fpemr = 1:gridLength
                fpelevmtx(fpemr,:) = lowelevation(fpemr,1); % detrend highelevations
            end           
            leveesuperelevation = highelevation - fpelevmtx;
        end
        
        % Make a plot of eta (Fig2) every once in a while (once every 'plotEvery' # of frames).
        if plotOngoingEta == 1 && (rem(j,plotEvery) == 0)
            figure(2)
            plot(x,eta,'b')
            hold on
            plot(x,eta0,'k:')
            hold off
            title([num2str(j*dt) ' years.'])
            pause(.01)
        end
        
        % update Fig1 if our settings say you should this timestep
        if animateSubsidence == 1 && (rem(j,plotSubEvery) == 0) && blackboxMode == 1
            figure(1)
            if colorSetting == 4
                imagesc(superelevation);
            elseif colorSetting == 5
                imagesc(leveesuperelevation);
            else
                %             error('Why are you trying to animate subsidence without showing elevations?')
                %             % Whatever.
            end
            colormap(customParula) % sets the colour
            caxis([seCbarLow seCbarHigh]); % sets the colourbar limits (for elevation) at the start
            axis equal % ensures that cells plot as squares, not rectangles
            axis xy
            axis ([0 300 0 300])
            colorbar % displays the colourbar legend
            drawnow % updates the frame without making a new window
        end
        
        % resets the list of viable avulsion loci (makes bugtesting
        % easier to do it before the loop, rather than after)
        viableLoci = 0*viableLoci;

        % resets avulsion trigger flag
        avulseThisTimestep = 0;
        
        % check to see if an avulsion occurs this timestep
        if triggersAreProbabilistic == 0 
            if mod(j,checkViableEvery) == 0
                avulseThisTimestep = 1;
            end
        elseif triggersAreProbabilistic == 1    
            taplrand = rand; % pick a number, any number...
            if taplrand <= (1/checkViableEvery)
                avulseThisTimestep = 1;
            end   
        end
        
        % if you do have a trigger this timestep...
        if avulseThisTimestep == 1
            triggerCount = triggerCount + 1;
            nextViable = 1;
            for pj = 1:(sum(pathHist(:,1)~=0)) % let's check to see which cells along the active channel might be viable avulsion loci
                checkedRow = pathHist(pj,1); % iterate along active channel
                checkedCol = pathHist(pj,2);
                
                % reset flags for diff directions
                SEedR = 0;
                SEedD = 0;
                SEedL = 0;
                SEedC = 0; % center is used for the absolute-aggradation settings that disallowInheritance
                
                if disallowInheritance == 0 % check relative elevation above neighboring cells
                    if (lowelevation(checkedRow,checkedCol) - lowelevation(checkedRow,checkedCol + 1)) >= minSuperelevation
                        SEedR = 1;
                    end
                    if (lowelevation(checkedRow,checkedCol) - lowelevation((checkedRow - 1),checkedCol)) >= minSuperelevation
                        SEedD = 1;
                    end
                    if (lowelevation(checkedRow,checkedCol) - lowelevation(checkedRow,(checkedCol - 1))) >= minSuperelevation
                        SEedL = 1;
                    end
                elseif disallowInheritance == 1 || disallowInheritance == 2 % see if you've accumulated enough superelevation in isolation
                    if SETmtx(checkedRow,checkedCol) >= minSuperelevation
                        SEedC = 1;
                    end
                else
                    error('invalid disallowInheritance selection; must be 0, 1, or 2.')
                end
                
                % okay, let's check to see if we can safely avulse from this cell...
                if checkedRow <= 4
                    % do nothing. We don't want to allow avulsions downwards here.
                elseif (SEedR == 1 || SEedC == 1) && ...
                        gridspace(checkedRow, (checkedCol + 1)) ~= (cGen)
                    % moving R is valid (elevation and generation)
                    if checkedCol < (gridWidth - 2) % otherwise, goes out of bounds
                        if ((gridspace((checkedRow - 1), (checkedCol + 2)) ~= 0 || gridspace((checkedRow - 1), (checkedCol + 1)) ~= 0 || gridspace((checkedRow - 1), checkedCol) ~= 0) || ...
                                (gridspace(checkedRow, (checkedCol + 2)) ~= 0 && gridspace((checkedRow - 1), (checkedCol + 3)) ~= 0) || ...
                                (pathStep(checkedRow, (checkedCol + 1)) < pathStep(checkedRow,(checkedCol + 2)))) % moving R does not lead to a terminal cell
                            if requireElevAdvantage == 0 % if you don't care about elevation advantage
                                viableLoci(nextViable,1) = pj; % saves old pathHist row number
                                viableLoci(nextViable,2) = checkedRow; % adds coords of loci to viableLoci list
                                viableLoci(nextViable,3) = checkedCol;
                                viableLoci(nextViable,5) = pj; % number of steps into the channel the viable is at
                                nextViable = nextViable + 1; % increments row number for viableLoci
                            elseif requireElevAdvantage == 1 % if we require a pseudo gradient advantage...
                                if (lowelevation(checkedRow,checkedCol) - lowelevation(checkedRow,checkedCol + 1)) > ...
                                        (lowelevation(checkedRow,checkedCol) - lowelevation(pathHist(pj+1,1),pathHist(pj+1,2)))
                                    viableLoci(nextViable,1) = pj; % saves old pathHist row number
                                    viableLoci(nextViable,2) = checkedRow; % adds coords of loci to viableLoci list
                                    viableLoci(nextViable,3) = checkedCol;
                                    viableLoci(nextViable,5) = pj; % number of steps into the channel the viable is at
                                    nextViable = nextViable + 1; % increments row number for viableLoci
                                end
                            else
                                error('Invalid requireElevAdvantage selection (must be 0 or 1).')
                            end
                        end
                    end
                elseif (SEedD == 1  || SEedC == 1)...
                        && gridspace(checkedRow - 1, (checkedCol)) ~= (cGen)
                    if requireElevAdvantage == 0
                        viableLoci(nextViable,1) = pj; % saves old pathHist row number
                        viableLoci(nextViable,2) = checkedRow; % adds coords of loci to viableLoci list
                        viableLoci(nextViable,3) = checkedCol;
                        viableLoci(nextViable,5) = pj; % number of steps into the channel the viable is at
                        nextViable = nextViable + 1; % increments row number for viableLoci
                    elseif requireElevAdvantage == 1 % if we require a pseudo gradient advantage....
                        if (lowelevation(checkedRow,checkedCol) - lowelevation(checkedRow - 1,checkedCol)) > ...
                                (lowelevation(checkedRow,checkedCol) - lowelevation(pathHist(pj+1,1),pathHist(pj+1,2)))
                            viableLoci(nextViable,1) = pj; % saves old pathHist row number
                            viableLoci(nextViable,2) = checkedRow; % adds coords of loci to viableLoci list
                            viableLoci(nextViable,3) = checkedCol;
                            viableLoci(nextViable,5) = pj; % number of steps into the channel the viable is at
                            nextViable = nextViable + 1; % increments row number for viableLoci
                        end
                    else
                        error('Invalid requireElevAdvantage selection (must be 0 or 1).')
                    end
                elseif (SEedL == 1 || SEedC == 1) ...
                        && gridspace(checkedRow, (checkedCol - 1)) ~= (cGen)
                    % moving L is valid (elevation and generation)
                    if checkedCol > 3 % otherwise, goes out of bounds
                        if ((gridspace((checkedRow - 1), (checkedCol - 2)) ~= 0 || gridspace((checkedRow - 1), (checkedCol - 1)) ~= 0 || gridspace((checkedRow - 1), checkedCol) ~= 0) || ...
                                (gridspace(checkedRow, (checkedCol - 2)) ~= 0 && gridspace((checkedRow - 1), (checkedCol - 3)) ~= 0) || ...
                                (pathStep(checkedRow, (checkedCol - 1)) < pathStep(checkedRow,(checkedCol - 2)))) % moving R does not lead to a terminal cell
                            if requireElevAdvantage == 0
                                viableLoci(nextViable,1) = pj; % saves old pathHist row number
                                viableLoci(nextViable,2) = checkedRow; % adds coords of loci to viableLoci list
                                viableLoci(nextViable,3) = checkedCol;
                                viableLoci(nextViable,5) = pj; % number of steps into the channel the viable is at
                                nextViable = nextViable + 1; % increments row number for viableLoci
                            elseif requireElevAdvantage == 1 % if we require a pseudo gradient advantage...
                                if (lowelevation(checkedRow,checkedCol) - lowelevation(checkedRow,checkedCol - 1)) > ...
                                        (lowelevation(checkedRow,checkedCol) - lowelevation(pathHist(pj+1,1),pathHist(pj+1,2)))
                                    viableLoci(nextViable,1) = pj; % saves old pathHist row number
                                    viableLoci(nextViable,2) = checkedRow; % adds coords of loci to viableLoci list
                                    viableLoci(nextViable,3) = checkedCol;
                                    viableLoci(nextViable,5) = pj; % number of steps into the channel the viable is at
                                    nextViable = nextViable + 1; % increments row number for viableLoci
                                end
                            else
                                error('Invalid requireElevAdvantage selection (must be 0 or 1).')
                            end
                        end
                    end
                end
            end
        end
        
        % if you found at least one viable avulsion locus, or if this is the first runthrough
        if sum(viableLoci(:,2)) ~= 0 || j == 1
            
            % if this isn't the first runthrough (not pretty, but it helped with "when do we initialize these variables?" I forget why, specifically.
            if j ~= 1
                
                % tracks viableLoci in a 2D accumulative space for Fig. 14
                for vltr = 1:sum(viableLoci(:,2)~=0)
                    viableLociTracker(viableLoci(vltr,2),viableLoci(vltr,3)) = viableLociTracker(viableLoci(vltr,2),viableLoci(vltr,3)) + 1;
                end
                
                % increment the avulsion generation
                cGen = cGen + 1;
                
                % pick one random entry from the list of viable loci
                newRandEntry = randperm(sum(viableLoci(:,1)~=0));
                currentEntry = viableLoci(newRandEntry(1),1);
                currentRow = viableLoci(newRandEntry(1),2);
                currentCol = viableLoci(newRandEntry(1),3); % sets the new avulsion locus as the starting point
                pickedStep = viableLoci(newRandEntry(1),5); % number of steps along the channel the viable is at
                
                % adds the new avulsion locus to the avulsionLoci table
                avulsionLoci(cGen,1) = j;
                avulsionLoci(cGen,2) = currentRow;
                avulsionLoci(cGen,3) = currentCol;
                avulsionLoci(cGen,4) = currentEntry;
                avulsionLoci(cGen,5) = pickedStep;
                
                % resets progradation length measurements
                progradLength(cGen,1) = 0;
                progradLength(cGen,2) = j;
                progradLength(cGen,3) = currentRow;
                
                % saves old data, for restoration if the avulsion fails
                oldStepDirections = stepDirections;
                oldPathHist = pathHist;
                oldPathStep = pathStep;
                oldGridspace = gridspace;
                oldRepulsionTracker = repulsionTracker;
                
                % increments the value of the new starting point
                gridspace(currentRow,currentCol) = cGen+1;
                didAvulse(j,1) = 1;
                
                % save a "before" for your superelevation matrix
                oldSETmtx = SETmtx;
                % remove the accumulated tracked superelevation for the cell that avulsed
                SETmtx(currentRow,currentCol) = 0;
                
                % writes over every stepDirections entry past the new starting point with zeroes.
                for sdt = currentEntry:(sum(pathHist(:,1)~=0))
                    stepDirections(sdt,1) = 0;
                    stepDirections(sdt,2) = 0;
                end
                
                % writes over every pathHist entry past the new ending point with zeroes.
                % if you don't do this, and an older path had more steps than your current
                % path, it can avulse from a historic channel-belt anomalously
                for k = (currentEntry+1):(sum(pathHist(:,1)~=0))
                    pathHist(k,1) = 0;
                    pathHist(k,2) = 0;
                end
                
                for m = 1:currentEntry % updates channel upstream of avulsion node to the newest generation (i.e. active)
                    gridspace(pathHist(m,1),pathHist(m,2)) = cGen;
                end
            end
            
            % Re-Initialized Variables:
            % these three are used for pathHist, in order to see if a successful move
            % was made or not (and thus if a new row in PathHist should be filled)
            currentStep = 1;
            lastRow = 0;
            lastCol = 0;
            
            %resetting stepsThisGen
            if j == 1 % first runthrough
                stepsThisGen = 1;
            else
                stepsThisGen = avulsionLoci(cGen,4);
            end
            
            for i = 1:numSteps %% PATHFINDING
                if i == numSteps - 1
                    error('You spent far too long pathfinding in an i-loop. You''re probably stuck.')
                end
                
                if j==1 && i==1 % if it's the first run-through...
                    currentCol = centerCol; % variables to track where the pointer is
                    currentRow = (gridLength - 1);
                    gridspace((gridLength - 1),centerCol) = 2; % set the top-centre as the starting point
                end
                
                % checks for out-of-bounds; ends if it hits a wall
                if (currentCol > 1) && (currentCol < gridWidth) && (currentRow > 1)
                    % && (currentRow < gridLength) <- I think I can safely not check for hitting the upper bound. If I'm wrong, I'll get a chuckle at least.
                    
                    atLeftWall = 0; % resetting border flags
                    atRightWall = 0;
                    
                    if currentCol == 3 % if at the left-most allowable cell
                        atLeftWall = 1;
                    elseif currentCol == gridWidth - 2 % if at the right-most allowable cell
                        atRightWall = 1;
                    elseif currentCol < 3 || currentCol > (gridWidth - 2)
                        error('You managed to get out of bounds, buddy.')
                    end
                    
                    % updates the current space to 'previously occupied during generation j'
                    gridspace(currentRow,currentCol) = cGen;
                    
                    % check to make sure you aren't duplicating an entry in pathHist!
                    if lastRow ~= currentRow || lastCol ~= currentCol
                        % adds the current location to the path history matrix
                        pathHist((currentEntry-1)+currentStep,1) = currentRow;
                        pathHist((currentEntry-1)+currentStep,2) = currentCol;
                        % increments currentStep
                        currentStep = currentStep + 1;
                    end
                    
                    % saves the current position so that it can be compared during the next loop
                    lastRow = currentRow;
                    lastCol = currentCol;
                    
                    % picks a random number, 0 to 1, in order to pick a PATHFINDING direction
                    randDirection = rand;
                    
                    % if this is the first step of an avulsion, only move to a non-cGenminus1 cell.
                    if j ~= 1 && i == 1
                        
                        % if using slope-weighted random walk, recalculate firstStepOddsR et al
                        if randomWalkMode == 1
                            % measures the change in elevation, adds infinitesimal noise to tiebreak
                            FSOdropR = lowelevation(currentRow,currentCol+1)-lowelevation(currentRow,currentCol)+rand*10e-10;
                            FSOdropD = lowelevation(currentRow-1,currentCol)-lowelevation(currentRow,currentCol)+rand*10e-10;
                            FSOdropL = lowelevation(currentRow,currentCol-1)-lowelevation(currentRow,currentCol)+rand*10e-10;
                            
                            % rank the directions
                            [FSOranked, FSOranks] = sort([FSOdropR FSOdropD FSOdropL]);
                            
                            % assign the weights accordingly
                            firstStepOddsR = FSOwalkWeights(FSOranks==1);
                            firstStepOddsD = FSOwalkWeights(FSOranks==2) + firstStepOddsR;
                            firstStepOddsL = FSOwalkWeights(FSOranks==3) + firstStepOddsD;
                        end
                        
                        % actually pick the next step
                        while newPathFound == 0 % we know at least one step is viable; keep trying it until we find one
                            randDirection2 = rand;
                            if (randDirection2 < (firstStepOddsR))
                                if gridspace(currentRow, (currentCol + 1)) < (cGen - 1) && atRightWall == 0 % checks to ensure it's avulsing into a non-active, non-out-of-bounds cell
                                    if ((gridspace((currentRow - 1), (currentCol + 2)) ~= 0 || gridspace((currentRow - 1), (currentCol + 1)) ~= 0 || gridspace((currentRow - 1), currentCol) ~= 0) || ...
                                            (gridspace(currentRow, (currentCol + 2)) ~= 0 && gridspace((currentRow - 1), (currentCol + 3)) ~= 0) || ...
                                            (pathStep(currentRow, (currentCol + 1)) < pathStep(currentRow,(currentCol + 2))))
                                        % moving right does not lead to a terminal cell
                                        if gridspace(currentRow, (currentCol + 1)) ~= 0 % checks to see if decided mvmt right goes into old path
                                            captured = 1;
                                        end
                                        currentCol = currentCol + 1; % moves right
                                        firstStepOldNum = gridspace(currentRow,currentCol);
                                        newPathFound = 1;
                                        stepDirections(stepsThisGen,1) = 1;
                                        stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                    else
                                        continue % if movement would go to a terminal route, try again
                                    end
                                else
                                    continue % if backtracking or rolling out of bounds, try again
                                end
                            elseif (randDirection2 < (firstStepOddsD))
                                if gridspace((currentRow - 1), currentCol) < (cGen - 1) % checks to ensure it's avulsing into a non-active cell
                                    if gridspace((currentRow - 1), currentCol) ~= 0 % checks to see if decided mvmt down goes into a historic channel
                                        captured = 1;
                                    end
                                    currentRow = currentRow - 1; % moves D
                                    firstStepOldNum = gridspace(currentRow,currentCol);
                                    newPathFound = 1;
                                    stepDirections(stepsThisGen,1) = 1;
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                    progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                else
                                    continue % if backtracking, try again
                                end
                            elseif (randDirection2 < (firstStepOddsL))
                                if gridspace(currentRow, (currentCol - 1)) < (cGen - 1) && atLeftWall == 0 % checks to ensure it's avulsing into a non-active, non-out-of-bounds cell
                                    if ((gridspace((currentRow - 1), (currentCol - 2)) ~= 0 || gridspace((currentRow - 1), (currentCol - 1)) ~= 0 || gridspace((currentRow - 1), currentCol) ~= 0) || ...
                                            (gridspace(currentRow, (currentCol - 2)) ~= 0 && gridspace((currentRow - 1), (currentCol - 3)) ~= 0) || ...
                                            (pathStep(currentRow, (currentCol - 1)) < pathStep(currentRow,(currentCol - 2))))
                                        % moving left does not lead to a terminal cell
                                        if gridspace(currentRow, (currentCol - 1)) ~= 0 % checks to see if decided mvmt left goes into a historic channel
                                            captured = 1;
                                        end
                                        currentCol = currentCol - 1; % moves L
                                        firstStepOldNum = gridspace(currentRow,currentCol);
                                        newPathFound = 1;
                                        stepDirections(stepsThisGen,1) = 1;
                                        stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                    else
                                        continue % if movement would go to a terminal route, try again
                                    end
                                else
                                    continue % if backtracking or rolling out of bounds, try again
                                end
                            else % you shouldn't be able to get this error.
                                error('Something went wrong with the first step avulsionOddsD/L/R being less than rand.');
                            end
                        end
                        newPathFound = 0; % resets the while-looping flag variable
                    
                    else % if this isn't the first step of an avulsion...
                        
                        if captured == 0 % if this ISN'T a REOCCUPATION...
                            
                            % reset the flags we're about to use below
                            noGoR = 0;
                            noGoDR = 0;
                            noGoD = 0;
                            noGoDL = 0;
                            noGoL = 0;
                            
                            if j~= 1 % don't do this the first generation, because you have no previous channels to resist reoccupying
                                % okay, so first, gather info on which of the 5 possible moves would lead to a resistant reoccupation
                                if currentCol >= gridWidth - 2 % if at the right wall
                                    noGoR = 1;
                                elseif (isChannel(currentRow, (currentCol + 1)) ~= 0 && ... % or if the cell is prev.occ. AND resistant
                                        (highelevation(currentRow, (currentCol + 1)) - lowelevation(currentRow,currentCol)) > depth(stepsThisGen)*resistanceFactor)
                                    noGoR = 1; % you can't go R!
                                    repulsionTracker(currentRow,1) = repulsionTracker(currentRow,1) + 1;
                                end
                                
                                if currentCol >= gridWidth - 2 % if at the right wall
                                    noGoDR = 1;
                                elseif (isChannel((currentRow - 1), (currentCol + 1)) ~= 0 && ... % or if the cell is prev.occ. AND resistant
                                        (highelevation((currentRow - 1), (currentCol + 1)) - lowelevation(currentRow,currentCol)) > depth(stepsThisGen)*resistanceFactor)
                                    noGoDR = 1; % you can't go DR!
                                    repulsionTracker(currentRow+1,1) = repulsionTracker(currentRow+1,1) + 1;
                                end
                                
                                if isChannel((currentRow - 1), currentCol) ~= 0 && ... % if the cell is prev.occ. AND resistant
                                        (highelevation((currentRow-1),currentCol) - lowelevation(currentRow,currentCol)) > depth(stepsThisGen)*resistanceFactor
                                    noGoD = 1; % you can't go D!
                                    repulsionTracker(currentRow+1,1) = repulsionTracker(currentRow+1,1) + 1;
                                end
                                
                                if currentCol <= 3 % if at the left wall
                                    noGoDL = 1;
                                elseif (isChannel((currentRow - 1), (currentCol - 1)) ~= 0 && ... % or if the cell is prev.occ. AND resistant
                                        (highelevation((currentRow - 1), (currentCol - 1)) - lowelevation(currentRow,currentCol)) > depth(stepsThisGen)*resistanceFactor)
                                    noGoDL = 1; % you can't go DL!
                                    repulsionTracker(currentRow+1,1) = repulsionTracker(currentRow+1,1) + 1;
                                end
                                
                                if currentCol <= 3 % if at the left wall
                                    noGoL = 1;
                                elseif (isChannel(currentRow, (currentCol - 1)) ~= 0 && ... % or if the cell is prev.occ. AND resistant
                                        (highelevation(currentRow, (currentCol - 1)) - lowelevation(currentRow,currentCol)) > depth(stepsThisGen)*resistanceFactor)
                                    noGoL = 1; % you can't go L!
                                    repulsionTracker(currentRow,1) = repulsionTracker(currentRow,1) + 1;
                                end
                                
                                % additionally, see if moving DR or DL would cross a historic channel diagonally
                                if isChannel(currentRow, (currentCol + 1)) ~= 0 && ... % R is historic
                                        isChannel((currentRow - 1), (currentCol + 1)) == 0 && ... % DR is virgin
                                        isChannel((currentRow - 1), currentCol) ~= 0 && ... D is historic
                                        gridspace(currentRow, (currentCol + 1)) ~= cGen % the 'blocker' on R isn't itself
                                    noGoDR = 1;
                                end
                                
                                if isChannel(currentRow, (currentCol - 1)) ~= 0 && ... % L is historic
                                        isChannel((currentRow - 1), (currentCol - 1)) == 0 && ... % DL is virgin
                                        isChannel((currentRow - 1), currentCol) ~= 0 && ... D is historic
                                        gridspace(currentRow, (currentCol - 1)) ~= cGen % the 'blocker' on L isn't itself
                                    noGoDL = 1;
                                end
                                
                                % okay, additionally, check if moving R or L would mean going upstream
                                if gridspace(currentRow, (currentCol + 1)) == cGen % check R
                                    noGoR = 1;
                                end
                                
                                if gridspace(currentRow, (currentCol - 1)) == cGen % check L
                                    noGoL = 1;
                                end
                                
                                % check to see if at least one move is viable
                                if noGoR == 1 && noGoDR == 1 && noGoD == 1 && noGoDL == 1 && noGoL == 1
                                    
                                    % you are trapped, and the avulsion fails!
                                    % restores old data
                                    stepDirections = oldStepDirections;
                                    pathHist = oldPathHist;
                                    pathStep = oldPathStep;
                                    gridspace = oldGridspace;
                                    repulsionTracker = oldRepulsionTracker;
                                    aTryTracker(j,1) = 1; % flags the avulsion as failed
                                    cGen = cGen - 1;
                                    SETmtx = oldSETmtx;
                                    break
                                    
                                end
                            end
                            
                            % if the cell is currently beside a historic channel, SNAPTUIT, unless doing so would lead to a resistant cell
                            if (gridspace((currentRow - 1), currentCol) ~= 0 && i > 2 && noGoD == 0) ... % checking D
                                    || currentRow == 2 % OR IF YOU'RE ADJACENT TO THE BOTTOM
                                captured = 1;
                                currentRow = currentRow - 1; % moves D
                                stepDirections(stepsThisGen,1) = 1;
                                stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                progradLength(cGen,1) = progradLength(cGen,1) + 1;
                            elseif gridspace((currentRow - 1), (currentCol + 1)) ~= 0 && i > 2 && noGoDR == 0 % checking DR
                                captured = 1;
                                currentRow = currentRow - 1; % moves DR
                                currentCol = currentCol + 1;
                                stepDirections(stepsThisGen,2) = 1;
                                stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                progradLength(cGen,1) = progradLength(cGen,1) + 1;
                            elseif gridspace((currentRow - 1), (currentCol - 1)) ~= 0 && i > 2 && noGoDL == 0 % checking DL
                                captured = 1;
                                currentRow = currentRow - 1; % moves DL
                                currentCol = currentCol - 1;
                                stepDirections(stepsThisGen,2) = 1;
                                stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                progradLength(cGen,1) = progradLength(cGen,1) + 1;
                            elseif gridspace(currentRow, (currentCol + 1)) ~= cGen && ... % don't check upstream cell
                                    gridspace(currentRow, (currentCol + 1)) ~= 0 && i > 2 && noGoR == 0 % checking R is historic
                                captured = 1;
                                currentCol = currentCol + 1; % moves R
                                stepDirections(stepsThisGen,1) = 1;
                                stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                progradLength(cGen,1) = progradLength(cGen,1) + 1;
                            elseif gridspace(currentRow, (currentCol - 1)) ~= cGen && ... % don't check upstream cell
                                    gridspace(currentRow, (currentCol - 1)) ~= 0 && i > 2 && noGoL == 0 % checking L is historic
                                captured = 1;
                                currentCol = currentCol - 1; % moves L
                                stepDirections(stepsThisGen,1) = 1;
                                stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                progradLength(cGen,1) = progradLength(cGen,1) + 1;
                            else % if no non-resistant snaptuit cells are present
                                % figures out which direction to move, executes it
                                
                                % if not weighting random walk by slope
                                if randomWalkMode == 0
                                    
                                    if j == 2
                                        oddsR = primeOddsR;
                                        oddsD = primeOddsD;
                                        oddsL = primeOddsL;
                                    end
                                    
                                elseif randomWalkMode == 1 % if slope-weighted
                                    
                                    % calculate drop in elevation in each direction
                                    dropR = lowelevation(currentRow,currentCol+1)-lowelevation(currentRow,currentCol)+rand*10e-10;
                                    dropDR = lowelevation(currentRow-1,currentCol+1)-lowelevation(currentRow,currentCol)+rand*10e-10;
                                    dropD = lowelevation(currentRow-1,currentCol)-lowelevation(currentRow,currentCol)+rand*10e-10;
                                    dropDL = lowelevation(currentRow-1,currentCol-1)-lowelevation(currentRow,currentCol)+rand*10e-10;
                                    dropL = lowelevation(currentRow,currentCol-1)-lowelevation(currentRow,currentCol)+rand*10e-10;
                                    
                                    % rank the drops
                                    [oddsRanked, oddsRanks] = sort([dropR dropDR dropD dropDL dropL]);
                                    
                                    % give weights accordingly
                                    oddsR = walkWeights(oddsRanks==1);
                                    oddsDR = walkWeights(oddsRanks==2);
                                    oddsD = walkWeights(oddsRanks==3);
                                    oddsDL = walkWeights(oddsRanks==4);
                                    oddsL = walkWeights(oddsRanks==5);
                                    
                                end
                                
                                % let's roll until we get a viable direction, now that we know at least one exists and we have them ranked:
                                if (randDirection < (oddsR)) && atRightWall == 0 % if you rolled R
                                    if gridspace(currentRow, (currentCol + 1)) ~= cGen && noGoR == 0 % checks to see if it's backtracking
                                        if gridspace(currentRow, (currentCol + 1)) ~= 0 % checks to see if decided mvmt R goes into old path
                                            captured = 1;
                                            progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                        end
                                        currentCol = currentCol + 1; % moves R
                                        stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                        stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                    else
                                        continue % if backtracking or resistant, try again
                                    end
                                elseif (randDirection < (oddsR+oddsDR)) && atRightWall == 0 % if you rolled DR
                                    if noGoDR == 0
                                        %                     if (gridspace((currentRow - 1), currentCol) == 0 || gridspace(currentRow, (currentCol + 1)) == 0) && noGoDR == 0
                                        if gridspace((currentRow - 1), (currentCol + 1)) ~= 0 % checks to see if decided mvmt DR goes into old path
                                            captured = 1;
                                            progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                        end
                                        currentRow = currentRow - 1; % moves DR
                                        currentCol = currentCol + 1;
                                        stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                        stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                    else
                                        continue % if DR is resistant
                                    end
                                elseif (randDirection < (oddsR+oddsDR+oddsD))
                                    if noGoD == 0
                                        if gridspace((currentRow - 1), currentCol) ~= 0 % checks to see if decided mvmt D goes into old path
                                            captured = 1;
                                            progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                        end
                                        currentRow = currentRow - 1; % moves D
                                        stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                        stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                    else
                                        continue % if D is resistant
                                    end
                                elseif (randDirection < (oddsR+oddsDR+oddsD+oddsDL)) && atLeftWall == 0
                                    if noGoDL == 0
                                        if gridspace((currentRow - 1), (currentCol - 1)) ~= 0 % checks to see if decided mvmt DL goes into old path
                                            captured = 1;
                                            progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                        end
                                        currentRow = currentRow - 1; % moves DL
                                        currentCol = currentCol - 1;
                                        stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                        stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                    else
                                        continue % if DL is resistant
                                    end
                                elseif (randDirection < (oddsR+oddsDR+oddsD+oddsDL+oddsL)) && atLeftWall == 0
                                    if gridspace(currentRow, (currentCol - 1)) ~= cGen && noGoL == 0 % checks to see if it's backtracking
                                        if gridspace(currentRow, (currentCol - 1)) ~= 0 % checks to see if decided mvmt L goes into old path
                                            captured = 1;
                                            progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                        end
                                        currentCol = currentCol - 1; % moves L
                                        stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                        stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        progradLength(cGen,1) = progradLength(cGen,1) + 1;
                                    else
                                        continue % if backtracking, try again
                                    end
                                elseif atLeftWall == 1 || atRightWall == 1
                                    continue % you rolled a direction that would have taken you out of bounds, so re-roll!
                                else
                                    % you really shouldn't be able to get this error
                                    error('Something went wrong with the cumulative odds being less than rand during pathfinding.');
                                end
                            end
                        elseif captured == 1 % if this is captured pathfinding...
                            
                            % this section gathers floodplain elevations for all neighbouring cells, so that they can be compared and used
                             % for steepest-descent pathfinding during reoccupation
                            if gridspace(currentRow, (currentCol + 1)) ~= 0 && ... % right isn't floodplain...
                                    ((gridspace(currentRow, (currentCol + 1)) == cGen && pathStep(currentRow,currentCol) < pathStep(currentRow,(currentCol + 1))) || ...
                                    (gridspace(currentRow, (currentCol + 1)) ~= cGen && ...
                                    ((gridspace((currentRow - 1), (currentCol + 2)) ~= 0 || gridspace((currentRow - 1), (currentCol + 1)) ~= 0 || gridspace((currentRow - 1), currentCol) ~= 0) || ...
                                    (gridspace(currentRow, (currentCol + 2)) ~= 0 && gridspace((currentRow - 1), (currentCol + 3)) ~= 0) || ...
                                    (pathStep(currentRow, (currentCol + 1)) < pathStep(currentRow,(currentCol + 2))))))
                                elevR = lowelevation(currentRow, (currentCol + 1));% ...and isn't terminal, then use its elevation
                            else
                                elevR = 9999; % otherwise, use an insurmountable elevation
                            end
                            
                            if gridspace((currentRow - 1), (currentCol + 1)) ~= 0 % DR isn't floodplain
                                elevDR = lowelevation((currentRow - 1), (currentCol + 1)); % then use its elevation
                            else
                                elevDR = 9999; % otherwise, use an insurmountable elevation
                            end
                            
                            if gridspace((currentRow - 1), currentCol) ~= 0 % D isn't floodplain
                                elevD = lowelevation((currentRow - 1), currentCol); % then use its elevation
                            else
                                elevD = 9999; % otherwise, use an insurmountable elevation
                            end
                            
                            if gridspace((currentRow - 1), (currentCol - 1)) ~= 0 % DL isn't floodplain
                                elevDL = lowelevation((currentRow - 1), (currentCol - 1)); % then use its elevation
                            else
                                elevDL = 9999; % otherwise, use an insurmountable elevation
                            end
                            
                            if gridspace(currentRow, (currentCol - 1)) ~= 0 && ... % left isn't floodplain...
                                    ((gridspace(currentRow, (currentCol - 1)) == cGen && pathStep(currentRow,currentCol) < pathStep(currentRow,(currentCol - 1))) || ...
                                    (gridspace(currentRow, (currentCol - 1)) ~= cGen && ...
                                    ((gridspace((currentRow - 1), (currentCol - 2)) ~= 0 || gridspace((currentRow - 1), (currentCol - 1)) ~= 0 || gridspace((currentRow - 1), currentCol) ~= 0) || ...
                                    (gridspace(currentRow, (currentCol - 2)) ~= 0 && gridspace((currentRow - 1), (currentCol - 3)) ~= 0) || ...
                                    (pathStep(currentRow, (currentCol - 1)) < pathStep(currentRow,(currentCol - 2))))))
                                elevL = lowelevation(currentRow, (currentCol - 1)); % ...and isn't terminal, then use its elevation
                            else
                                elevL = 9999; % otherwise, use an insurmountable elevation
                            end
                            
                            % now that we have elevations for all of the neighbours, we'll look at which direction we came from, sort all resulting
                            % non-upstream options to find the greatest descent, and move there
                            if gridspace(currentRow, (currentCol + 1)) == cGen % if you just came from R
                                lowestMoveElev = min([elevDR,elevD,elevDL,elevL]); % sort the elevations of all possible moves
                                
                                if lowestMoveElev == 9999 % catches if the run is about to mess up
                                    captured = 0; % eject from the channelbelt, back to floodplain pathfinding
                                    continue % go back to determine your next (unconfined) step
                                end
                                
                                if lowestMoveElev == elevDR && gridspace((currentRow - 1), (currentCol + 1)) ~= 0
                                    currentRow = currentRow - 1;
                                    currentCol = currentCol + 1; % moves DR
                                    stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                elseif lowestMoveElev == elevD && gridspace((currentRow - 1), currentCol) ~= 0
                                    currentRow = currentRow - 1; % moves D
                                    stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                elseif lowestMoveElev == elevDL && gridspace((currentRow - 1), (currentCol - 1)) ~= 0
                                    currentRow = currentRow - 1;
                                    currentCol = currentCol - 1; % moves DL
                                    stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                elseif lowestMoveElev == elevL && gridspace(currentRow, (currentCol - 1)) ~= 0
                                    currentCol = currentCol - 1; % moves L
                                    stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                else
                                    error('Something went wrong after coming from the left during reoccupation pathfinding.')
                                end
                            elseif gridspace(currentRow, (currentCol - 1)) == cGen % if you just came from L
                                
                                lowestMoveElev = min([elevR,elevDR,elevD,elevDL]); % sort the elevations of all possible moves
                                
                                if lowestMoveElev == 9999 % catches if the run is about to mess up
                                    captured = 0; % eject from the channelbelt, back to floodplain pathfinding
                                    continue % go back to determine your next (unconfined) step
                                end
                                
                                if lowestMoveElev == elevR && gridspace(currentRow, (currentCol + 1)) ~= 0
                                    currentCol = currentCol + 1; % moves R
                                    stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                elseif lowestMoveElev == elevDR && gridspace((currentRow - 1), (currentCol + 1)) ~= 0
                                    currentRow = currentRow - 1;
                                    currentCol = currentCol + 1; % moves DR
                                    stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                elseif lowestMoveElev == elevD && gridspace((currentRow - 1), currentCol) ~= 0
                                    currentRow = currentRow - 1; % moves D
                                    stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                elseif lowestMoveElev == elevDL && gridspace((currentRow - 1), (currentCol - 1)) ~= 0
                                    currentRow = currentRow - 1;
                                    currentCol = currentCol - 1; % moves DL
                                    stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                else
                                    error('Something went wrong after coming from the left during reoccupation pathfinding.')
                                end
                            elseif gridspace((currentRow + 1), (currentCol - 1)) == cGen ... % if you just came from UL, U, or UR
                                    || gridspace((currentRow + 1), (currentCol)) == cGen ...
                                    || gridspace((currentRow + 1), (currentCol + 1)) == cGen
                                lowestMoveElev = min([elevR,elevDR,elevD,elevDL,elevL]); % sort the elevations of all possible moves
                                
                                if lowestMoveElev == 9999 % catches if the run is about to mess up
                                    captured = 0; % eject from the channelbelt, back to floodplain pathfinding...
                                    continue % go back to determine your next (unconfined) step
                                end
                                
                                if lowestMoveElev == elevD && gridspace((currentRow - 1), currentCol) ~= 0 % checks down first (in case of ties)
                                    currentRow = currentRow - 1; % moves D
                                    stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                    stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                elseif lowestMoveElev == elevDR || lowestMoveElev == elevDL % if down isn't steepest or viable, check DL and DR
                                    if elevDR ~= elevDL % if DR & DL are unequal elevation, pick the lower
                                        if lowestMoveElev == elevDR && gridspace((currentRow - 1), (currentCol + 1)) ~= 0
                                            currentRow = currentRow - 1;
                                            currentCol = currentCol + 1; % moves DR
                                            stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                            stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        elseif lowestMoveElev == elevDL && gridspace((currentRow - 1), (currentCol - 1)) ~= 0
                                            currentRow = currentRow - 1;
                                            currentCol = currentCol - 1; % moves DL
                                            stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                            stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        else
                                            error('Error during reoccupation pathfinding, after it moved down and elevations were unequal.')
                                        end
                                    elseif elevDR == elevDL % if DR and DL are equal elevation, roll to see which you pick
                                        if (randDirection < capturedOddsDL) % rolls DL, so picks DL
                                            if lowestMoveElev == elevDL && gridspace((currentRow - 1), (currentCol - 1)) ~= 0 % if DL is lower elev
                                                currentRow = currentRow - 1;
                                                currentCol = currentCol - 1; % moves DL
                                                stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                                stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                            else
                                                error('Reoccupation pathfinding messed up after rolling left.')
                                            end
                                        elseif (randDirection < capturedOddsDR) % rolls DR, so picks DR
                                            if lowestMoveElev == elevDR && gridspace((currentRow - 1), (currentCol + 1)) ~= 0 % if DR is lower elev
                                                currentRow = currentRow - 1;
                                                currentCol = currentCol + 1; % moves DR
                                                stepDirections(stepsThisGen,2) = 1; % we moved diagonally this timestep
                                                stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                            else
                                                error('Reoccupation pathfinding messed up after rolling right.')
                                            end
                                        else
                                            % do nothing
                                        end
                                    end
                                elseif lowestMoveElev == elevR || lowestMoveElev == elevL % finally, if R or L are lowest
                                    if elevR ~= elevL % if R & L are unequal elevation, pick the lower
                                        if lowestMoveElev == elevR && gridspace(currentRow, (currentCol + 1)) ~= 0
                                            currentCol = currentCol + 1; % moves R
                                            stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                            stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        elseif lowestMoveElev == elevL && gridspace(currentRow, (currentCol - 1)) ~= 0
                                            currentCol = currentCol - 1; % moves L
                                            stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                            stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                        else
                                            error('Error during reoccupation pathfinding, after it moved down and elevations were unequal.')
                                        end
                                    elseif elevR == elevL % if R and L are equal elevation, roll to see which you pick
                                        if (randDirection < capturedOddsL) % rolls L, so picks L
                                            if lowestMoveElev == elevL && gridspace(currentRow, (currentCol - 1)) ~= 0 % if L is lower elev
                                                currentCol = currentCol - 1; % moves L
                                                stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                                stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                            else
                                                error('Reoccupation pathfinding messed up after rolling left.')
                                            end
                                        elseif (randDirection < capturedOddsR) % rolls R, so checks R before L
                                            if lowestMoveElev == elevR && gridspace(currentRow, (currentCol + 1)) ~= 0 % if R is lower elev
                                                currentCol = currentCol + 1; % moves R
                                                stepDirections(stepsThisGen,1) = 1; % we moved orthogonally this timestep
                                                stepsThisGen = stepsThisGen + 1; % add one to the count for successful moves
                                            else
                                                error('Reoccupation pathfinding messed up after rolling right.')
                                            end
                                        else
                                            error('Error during reoccupation pathfinding, after it moved down and elevations were equal.')
                                        end
                                    else
                                        error('Error in captured pathfinding where R or L are the lowest moves.')
                                    end
                                else
                                    error('Did we run out of options during captured pathfinding after a downward move?')
                                end
                            else
                                error('Apparently, we didn''t come from L, UL, U, UR, or R during captured pathfinding.')
                            end
                        else
                            error('Captured should have been 1 or 0, but it wasn''t.')
                        end % end of 'If captured == 0 elseif captured == 1'
                    end
                    
                    % sets the pathStep value for the current cell (upstream vs downstream)
                    pathStep(currentRow,currentCol) = stepCounter;
                    stepCounter = stepCounter + 1;
                    
                else % this code is what happens if you hit a wall; you update and avulse
                    gridspace(currentRow,currentCol) = cGen; %updates the final space to j
                    
                    if captured == 0
                        progradLength(cGen,1) = progradLength(cGen,1) + 1;
                    end
                    
                    break
                end
                
                % updates current space to j+1
                gridspace(currentRow,currentCol) = cGen+1;
                
                % updates Fig1, if our settings suggest we would have to this timestep
                if blackboxMode == 1 && (animateFloodplainPathfinding == 1 || (j ~= 1 && i==1 && aTryTracker(j-1,1) ~= 1))
                    figure(1)
                    if animateCapturedPathfinding == 0
                        if captured == 0 % only animates when navigating through new cells, not for reoccupation
                            if colorSetting == 4
                                imagesc(superelevation); % plots the next frame
                                colormap bone % sets the colour
                                caxis([seCbarLow seCbarHigh]); % sets the colourbar limits (for elevation) at the start
                            elseif colorSetting == 5
                                imagesc(leveesuperelevation); % plots the next frame
                                colormap(customParula) % sets the colour
                                caxis([seCbarLow seCbarHigh]); % sets the colourbar limits (for elevation) at the start
                            else
                                error('Invalid colour bar setting.')
                            end
                            axis equal % ensures that cells plot as squares, not rectangles
                            axis xy
                            axis ([0 300 0 300]) % sets the bounds
                            colorbar % displays the colourbar legend
                            drawnow % updates the frame without making a new window
                        end
                    else
                        if colorSetting == 4
                            imagesc(superelevation); % plots the next frame
                            colormap bone % sets the colour
                            caxis([seCbarLow seCbarHigh]); % sets the colourbar limits
                        elseif colorSetting == 5
                            imagesc(leveesuperelevation); % plots the next frame
                            colormap lines % sets the colour
                            caxis([seCbarLow seCbarHigh]); % sets the colourbar limits (for elevation) at the start
                        else
                            error('Invalid colour bar setting.')
                        end
                        axis equal % ensures that cells plot as squares, not rectangles
                        axis xy
                        axis ([0 300 0 300]) % sets the bounds
                        colorbar % displays the colourbar legend
                        drawnow % updates the frame without making a new window
                    end
                end
            end %ends the whole river pathfinding loop, ready to avulse
        end

        % sets up our 1D diffusion problem:       
        % set the number of grid points and build a cell-center grid
        N = (sum(pathHist(:,1)~=0)); % number of cells, recalculated based on new river path length
        
        % establish the new unit length based on the new total length
        numOrthos = (sum(stepDirections(:,1)~=0));
        numDiags = (sum(stepDirections(:,2)~=0));
        
        % updates spatial step to be the new average step length, accounting for orthogonal and diagonal steps
        L = numOrthos * orthoSize + numDiags * diagSize;
        dx = L/N;
        Ltracker(j) = L;
        
        % saves the old x
        x0 = x;
        x = -.5*dx:dx:L+.5*dx; % the 1D set of spaced gridpoints in the 'x' direction
        x = x'; % turn x into a column vector.
        
        % redefine sigma (can't preserve between runs b/c changing length)
        sigma = zeros(sum(pathHist(:,1)~=0)+2,1);
        
        % sets sigma (looking up subsidence along betaSubMtx and mapping it by row)
        if subsidence > 0
            sigma(1) = ghostSigma;
            sigma((sum(pathHist(:,1)~=0))+2) = ghostSigma2;
            for sgc = 2:(sum(pathHist(:,1)~=0)+1)
                sigma(sgc,1) = betaSubMtx(pathHist(sgc-1,1),1);
            end
        else
            sigma = zeros(length(N)+2);
        end
        
        q = zeros(length(x),1)+qo; % specific discharge of water, m^2/yr, i.e. rate*basin length. For L=100,000 m, and rain = 1 m/yr, q = 100,000 m2/yr
        D =(8.*q.*Acoeff.*sqrt(cf))/(C0*(2.65-1)); % fluvial diffusivity, m^2/yr
        
        % Load Dm with average values D(j-1/2) and Dp with D(j+1/2)
        Dm = zeros(N+2,1);
        Dp = zeros(N+2,1); % Make the column vectors
        Dm(2:N+1) = .5*(D(2:N+1)+D(1:N)); % average j and j-1
        Dp(2:N+1) = .5*(D(2:N+1)+D(3:N+2)); % average j and j+1
        
        % for debugging, keep old eta (for comparison)
        eta0 = eta;
        ghostElev2 = eta(length(eta));
        
        % redefine eta (can't preserve between runs b/c changing length)
        ghostElev = eta(1);
        eta = zeros(sum(pathHist(:,1)~=0)+2,1);
        eta(1) = ghostElev;
        eta(N+2) = ghostElev2;
        
        % load the newest elevations along the active channel path
        for ee = 2:sum(pathHist(:,1)~=0)+1
            eta(ee) = lowelevation(pathHist(ee-1,1),pathHist(ee-1,2));
        end
        
        % Determine channel slope
        for slp = 1:(N-1)
            slope(slp) = abs((eta(slp+1) - eta(slp)) / dx);
        end
        
        % sets the slope for the final (unused) cell equal to the penultimate one
        slope(N) = slope(N-1);
        
        % QC channel slope by constraining outliers
        for slp2 = 1:(N)
            if slope(slp2) > maxslope
                slope(slp2) = maxslope;
            elseif slope(slp2) < minslope
                slope(slp2) = minslope;
            end
        end
        
        Q = zeros(N,1)+(qo*basinwidth/3.154e7);
        
        % do solve for river parameters along the active channel
        for dpt = 1:(sum(pathHist(:,1)~=0))
            
            % solving via shields
            To(dpt) = shields * ((rohS - rohW) * g * D50); % as of now,
            % this doesn't actually vary along-channel, so we could calculate it once outside of the loop for
             % computational efficiency (if this were the 90s)
            
            % calculate velocity via empirical bed shear stress eq'n
            velocity(dpt) =  sqrt(To(dpt)/(cf*rohW)); % same thing
            
            % calculate depth via open channel analytical shear stress eq'n
            depth(dpt) = To(dpt) / (rohW * g * slope(dpt));
            
            % calculate depth via definition of discharge
            chanwidth(dpt) = Q(dpt)/(depth(dpt)*velocity(dpt));
            
        end
        
        % writes over slope, depth, velocity, and width values past the new ending point with zeroes.
        for sow = (N+1):max([(sum(slope(:,1)~=0)),(sum(velocity(:,1)~=0)),(sum(depth(:,1)~=0)),(sum(chanwidth(:,1)~=0))])
            slope(sow) = 0;
            velocity(sow) = 0;
            depth(sow) = 0;
            chanwidth(sow) = 0;
        end
        
        % writes over later values of depth using the last value of depth
        lastDepth = depth(N);
        for unsow = (N+1):gridLength*2
            depth(unsow) = lastDepth;
        end
        
        if incisionRule == 1 % incises new cells down one channel depth from HIGHELEV
            if j == 1 % if you just picked your first path
                for dsir2 = 1:sum(pathHist(:,1)~=0)
                    lowelevation(pathHist(dsir2,1),pathHist(dsir2,2)) = highelevation(pathHist(dsir2,1),pathHist(dsir2,2)) - depth(dsir2);
                end
                % load the newest elevations along the active channel path
                for ee = 2:sum(pathHist(:,1)~=0)+1
                    eta(ee) = lowelevation(pathHist(ee-1,1),pathHist(ee-1,2));
                end
            elseif didAvulse(j) == 1 && aTryTracker(j,1) ~= 1 % if you just avulsed
                for dsir = (avulsionLoci(cGen,4)+1):sum(pathHist(:,1)~=0) % for all cells downstream of the avulsion node
                    % if incising one channel depth down from the highelev would lower
                     % the cell (i.e. cell is not already a confined valley)
                    if (highelevation(pathHist(dsir,1),pathHist(dsir,2)) - depth(dsir)) < lowelevation(pathHist(dsir,1),pathHist(dsir,2))
                        lowelevation(pathHist(dsir,1),pathHist(dsir,2)) = highelevation(pathHist(dsir,1),pathHist(dsir,2)) - depth(dsir);
                    end
                end
                % load the newest elevations along the active channel path
                for ee = 2:sum(pathHist(:,1)~=0)+1
                    eta(ee) = lowelevation(pathHist(ee-1,1),pathHist(ee-1,2));
                end
            elseif aTryTracker(j,1) == 1
                % error('You apparently just pathfound without avulsing?')
            end
        else
            error('newIncisionRule must be set to 0 or 1.')
        end
        
        % update Fig1 if appropriate
        if blackboxMode == 2 && rem(j,animateEvery) == 0 % only every X frames!
            
            figure(1)
            if colorSetting == 4
                imagesc(superelevation); % plots the next frame
                colormap bone % sets the colour
                caxis([seCbarLow seCbarHigh]); % sets the colourbar limits (for elevation) at the start
            elseif colorSetting == 5
                imagesc(leveesuperelevation); % plots the next frame
                colormap(customParula) % sets the colour
                caxis([seCbarLow seCbarHigh]); % sets the colourbar limits (for elevation) at the start
            else
                error('Invalid colour bar setting.')
            end
            axis equal % ensures that cells plot as squares, not rectangles
            axis xy
            axis([0 300 0 300])
            colorbar % displays the colourbar legend
            drawnow % updates the frame without making a new window
            title(['Elapsed Time: ',num2str(dt*j),' years'])
            xlabel('Distance (units of 500m)')
            ylabel('Distance from mountainfront (units of 500m)')
            hcbar = colorbar;
            set(get(hcbar,'label'),'string','Elevation above detrended floodplain (m)');
            
            % create a fancy animated GIF showing the model run
            if createGIF == 1
                
                currentFrame = getframe(fig1); % look at Fig1
                currentIm = frame2im(currentFrame); % convert it to Im
                [imIndexed,cm] = rgb2ind(currentIm,256); % convert it to RGB index w/ 256 colours
                
                if j == animateEvery % if this is the first frame being written
                    imwrite(imIndexed,cm,fullfile(exportDirectory,exportFolder,[GIFfilename,'.gif']),'gif','DelayTime',0.1,'Loopcount',inf);
                else % if this is any other frame being written
                    imwrite(imIndexed,cm,fullfile(exportDirectory,exportFolder,[GIFfilename,'.gif']),'gif','DelayTime',0.1,'WriteMode','append');
                end
                
            end
            
            % if you are saving pictures of Fig1 over time
            if exportOutput == 1 && rem(j,exportEvery) == 0
                exportFilename = [exportName,'_Year_',num2str(j*10),'.png'];
                saveas(fig1,fullfile(exportDirectory,exportFolder,exportFilename))
            end
            
        end
        
        % get back to setting up the 1D diffusion model
        % (sediment variables?)
        const = 2*dx^2 / dt;
        eq_slope = qsin/D(1);
        
        % Create the sparse matrices A and B
        aneg1 = [-Dm(2:N+1); 0;0]; %extra zeros are needed to pad row N+2 for boundary conditions
        a0 = [0;const + (Dm(2:N+1)+Dp(2:N+1));0]; %extra zeros are needed to pad row N+2 for boundary conditions
        a1 = [0;0; -Dp(2:N+1)];%extra zeros are needed to pad row N+2 for boundary conditions
        A = spdiags([aneg1 a0 a1],[-1 0 1],N+2,N+2);
        
        bneg1 = [Dm(2:N+1); 0;0]; %extra zeros are needed to pad row N+2 for boundary conditions
        b0 = [0;const - (Dm(2:N+1)+Dp(2:N+1));0]; %extra zeros are needed to pad row N+2 for boundary conditions
        b1 = [0;0; Dp(2:N+1)];%extra zeros are needed to pad row N+2 for boundary conditions
        B = spdiags([bneg1 b0 b1],[-1 0 1],N+2,N+2);
        
        % load the boundary conditions into A and B A(1,1)=1 and A(1,2)=1 makes first equation eta1-eta2=r(1), 
         % this means that r(1) is the elevation difference between the cells which we set as the equilibrium slope 
         % needed to transport sediment supplied.
        A(1,1) = 1; 
        A(1,2) = -1; 
        B(1,1) = 0; 
        % T(0) = 0
        A(N+2,N+1) = 0.5; 
        A(N+2,N+2) = 0.5; 
        B(N+2,N+2) = 0; 
        % T(L)=0
        
        %resets captured to 0
        captured = 0;
    end % end of the avulsion/timestep loop!!
    
    % start of what we do when we complete the final timestep
     % update figure 1
    if colorSetting == 4 % only spend the memory/time to make and calculate these matrices if you're gonna use 'em.
        for fpemr = 1:gridLength
            fpelevmtx(fpemr,:) = lowelevation(fpemr,1);
        end
        superelevation = lowelevation - fpelevmtx;
    elseif colorSetting == 5
        for fpemr = 1:gridLength
            fpelevmtx(fpemr,:) = lowelevation(fpemr,1);
        end
        leveesuperelevation = highelevation - fpelevmtx;
    end
    
    figure(1)
     % final animation step
    if colorSetting == 4
        imagesc(superelevation); % plots the next frame
        colormap(customParula) % sets the colour
        caxis([seCbarLow seCbarHigh]); % sets the colourbar limits (for elevation) at the start
    elseif colorSetting == 5
        imagesc(leveesuperelevation); % plots the next frame
        colormap(customParula) % sets the colour
        caxis([seCbarLow seCbarHigh]); % sets the colourbar limits (for elevation) at the start
    else
        error('Invalid colour bar setting.')
    end
    axis equal % ensures that cells plot as squares, not rectangles
    axis xy
    axis ([0 300 0 300]) % figure boundaries
    colorbar % displays the colourbar legend
    drawnow % updates the frame without making a new window
    
    % housekeeping: remove trailing zeros from some variables
    finalAvulsionLoci = zeros(sum(avulsionLoci(:,1)~=0)+1,3);
    finalprogradLength = zeros(sum(avulsionLoci(:,1)~=0)+1,3);
    
    % likewise, cleans up any unused rows in avulsionLoci before plotting
    for al = 1:sum(avulsionLoci(:,1)~=0)+1
        finalAvulsionLoci(al,1) = avulsionLoci(al,1);
        finalAvulsionLoci(al,2) = avulsionLoci(al,2);
        finalAvulsionLoci(al,3) = avulsionLoci(al,5);
    end
    
    % also cleans up any unused rows in progradLength before plotting
    for pl1 = 1:sum(progradLength(:,1)~=0)
        finalprogradLength(pl1,3) = progradLength(pl1,3); % currentRow
        finalprogradLength(pl1,2) = progradLength(pl1,2); % j
        finalprogradLength(pl1,1) = progradLength(pl1,1); % progradation length
    end
    
    % math for calculating width of the megafan
    % initialize variables
    Ncc = zeros(1,gridWidth);
    for ncf = 1:gridWidth
        Ncc(ncf) = nnz(gridspace(3:gridLength,ncf));
    end
    NccLeft = Ncc(1:(ceil(gridWidth/2)));
    NccRight = Ncc(gridWidth:-1:(ceil(gridWidth/2)));
    
    % calculate the left boundary of the fan
    totchLeft = sum(NccLeft);
    norchLeft = cumsum(NccLeft)/totchLeft;
    boundLeft = find(norchLeft>0.1,1); % 90th percentile of just the left half of the data (so 95th overall)
    
    % calculate the left boundary of the fan
    totchRight = sum(NccRight);
    norchRight = cumsum(NccRight)/totchRight;
    boundRight = gridWidth-find(norchRight>0.1,1); % 90th percentile of just the left half of the data (so 95th overall)
    
    % calculate megafan width
    megafanWidth = (boundRight-boundLeft)*cellSize/1000; % in km
    
    % math for calculating length of the megafan
     % initialize variables
    Ncr = zeros((gridLength-2),1);
    for nrf = 1:(gridLength-2)
        Ncr(nrf) = nnz(gridspace(nrf+2,:)); 
    end
    
    % calculate the downstream boundary of the fan)
    totchDown = sum(Ncr);
    norchDown = cumsum(Ncr)/totchDown;
    boundDown = (gridLength-2)-find(norchDown>0.2,1); % 80th percentile downstream for the whole data
    
    % calculate megafan length
    megafanLength = boundDown * cellSize/1000; % in km
    
    % we already created and drew Fig1 earlier. Now we just update the title with some final stats
    figure(1)
    title(['W: ',num2str(megafanWidth), 'km, L:', num2str(megafanLength), 'km, R:', num2str((megafanWidth/megafanLength)), ', A: ', num2str(nnz(gridspace)/4), ' sq.km'])
    
    % likewise, Fig2 is already drawn. Just update the labels
    fig2 = figure(2);
    title('Final river profile')
    xlabel('along-channel distance (m)')
    ylabel('elevation (m)')
    
    % create Figure 3, which shows four panels as subplots:
    figure(3)
    fig3 = figure(3);
    set(fig3,'Units','normalized');
    set(fig3,'Position',[0.5 0 1 1]); % configured for my desktop monitor setup
    
    if sum(avulsionLoci(:,1)~=0) ~= 0 % for debugging, lets runs work with no avulsions
        if plotAvulsionLoci == 1
            subplot(2,2,1) % plots avulsion distance from the mountain-front through time
            plot(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,1),finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,2),'b.')
            hold on
            title('y-values of avulsion loci vs avulsion #')
            xlabel('timestep #')
            ylabel('y-values of avulsion loci')
            axis([1 timeSteps 0 gridLength])
        end
    end
    
    subplot(2,2,2) % plots proximal/apex elevations over time
    plot(elevTracker(:,1))
    hold on
    title('Elevations of eta(1)[blue] and eta(2)[red]')
    xlabel('timestep #')
    ylabel('elevation (m)')
    plot(elevTracker(:,2),'r')
    axis([1 timeSteps 0 700])
    
    subplot(2,2,3) % plots distal/exit elevations over time
    plot(elevTracker(:,3),'g')
    hold on
    title('Elevations of eta(N+1)[green] and eta(N+2)[black]')
    xlabel('timestep #')
    ylabel('elevation (m)')
    plot(elevTracker(:,4),'k')
    axis([1 timeSteps -2 2])
    
    subplot(2,2,4) % plots river length over time
    plot(Ltracker)
    title('River length by timestep')
    xlabel('timestep #')
    ylabel('River length (m)')
    axis([1 timeSteps 0 400000])
    
    % prepare some variables re: avulsion locations for following figures
    mov95_5k = zeros(length(1:sum(finalAvulsionLoci(:,1)~=0)),1); % initialize
    
    for lpp = 5001:sum(finalAvulsionLoci(:,1)~=0)-5000 % 50kyr moving window average of the 95th percentile
        mov95_5k(lpp) = prctile(finalAvulsionLoci(lpp-5000:lpp+5000,2),90);
    end
    
    mov05_5k = zeros(length(1:sum(finalAvulsionLoci(:,1)~=0)),1);
    
    for lpp = 5001:sum(finalAvulsionLoci(:,1)~=0)-5000 % 50kyr moving window average of the 5th percentile
        mov05_5k(lpp) = prctile(finalAvulsionLoci(lpp-5000:lpp+5000,2),10);
    end
    
    mov95_500 = zeros(length(1:sum(finalAvulsionLoci(:,1)~=0)),1);
    
    for lpp = 501:sum(finalAvulsionLoci(:,1)~=0)-500 % 5kyr moving window average of the 95th percentile
        mov95_500(lpp) = prctile(finalAvulsionLoci(lpp-500:lpp+500,2),90);
    end
    
    mov05_500 = zeros(length(1:sum(finalAvulsionLoci(:,1)~=0)),1);
    
    for lpp = 501:sum(finalAvulsionLoci(:,1)~=0)-500 % 5kyr moving window average of the 5th percentile
        mov05_500(lpp) = prctile(finalAvulsionLoci(lpp-500:lpp+500,2),10);
    end
    
    fig4 = figure(4); % a bigger/detailed-er version of fig 2.1, the avulsion distances from mountain-front over time
    set(fig4,'Units','normalized');
    set(fig4,'Position',[0 0 1.0 0.5]); % configured for my desktop monitor setup
    subplot(1,2,1) % plots the datapoints
    plot(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,1),finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,2),'b.')
    hold on
    title('y-values of avulsion loci vs avulsion #')
    xlabel('timestep #')
    ylabel('y-values of avulsion loci')
    axis([1 timeSteps 0 gridLength])
    subplot(1,2,2) % just the trendlines
    plot(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,1),movmean(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,2),501),'b:')
    hold on
    plot(finalAvulsionLoci(501:sum(finalAvulsionLoci(:,1)~=0)-500,1),mov05_500(501:sum(finalAvulsionLoci(:,1)~=0)-500,1),'k:')
    plot(finalAvulsionLoci(501:sum(finalAvulsionLoci(:,1)~=0)-500,1),mov95_500(501:sum(finalAvulsionLoci(:,1)~=0)-500,1),'r:')
    plot(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,1),movmean(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,2),5001),'b')
    plot(finalAvulsionLoci(5001:sum(finalAvulsionLoci(:,1)~=0)-5000,1),mov05_5k(5001:sum(finalAvulsionLoci(:,1)~=0)-5000,1),'k')
    plot(finalAvulsionLoci(5001:sum(finalAvulsionLoci(:,1)~=0)-5000,1),mov95_5k(5001:sum(finalAvulsionLoci(:,1)~=0)-5000,1),'r')
    
    title('y-values of avulsion loci vs avulsion #')
    xlabel('timestep #')
    ylabel('y-values of avulsion loci')
    axis([1 timeSteps 0 gridLength])
    
    fig5 = figure(5); % shows a sample down-domain transect along the domain centreline of both levee and bed elevations
    set(fig5,'Units','normalized');
    set(fig5,'Position',[0 0 0.75 0.75]); % configured for my desktop monitor setup
    plot(highelevation(gridLength-1:-1:2,ceil(gridWidth/2)))
    hold on
    plot(lowelevation(gridLength-1:-1:2,ceil(gridWidth/2)))
    % hold off
    legend('leveetop elevation','channelbed elevation')
    title('Elevation profile along megafan dip')
    xlabel(['distance (',num2str(cellSize),'m cells)'])
    ylabel('elevation (m)')
    
    fig6 = figure(6); % histograms for where avulsions occured
    set(fig6,'Units','normalized');
    set(fig6,'Position',[0 0 1.0 0.5]); % configured for my desktop monitor setup
    subplot(1,2,1) % 2d histogram
    avulsion2DHist = histcounts2([avulsionLoci(2:sum(avulsionLoci(:,1)~=0),2);1;gridLength;gridLength;1],[avulsionLoci(2:sum(avulsionLoci(:,1)~=0),3);gridWidth;gridLength;1;1],[gridWidth gridLength]);
    imagesc(avulsion2DHist); % plots the next frame
    colormap(flipud(bone))
    axis equal % ensures that cells plot as squares, not rectangles
    axis xy
    caxis([0 100]) % need to find a way to make this variable
    colorbar % displays the colourbar legend
    title('2D avulsion histogram')
    xlabel(['distance (',num2str(cellSize),'m cells)'])
    xlim([0 300])
    ylim([0 300])
    
    fig6;
    subplot(1,2,2) % 1d histogram, going down-domain (by row, binned)
    binWidth = 10;
    numBins = ceil((gridLength-2)/binWidth);
    probhist = histogram(gridLength-avulsionLoci(2:sum(avulsionLoci(:,1)~=0),2)-1,numBins,'DisplayStyle','stairs');
    xlim([0 gridLength-2])
    ylim([0 1.001*max(probhist.Values())])
    hold on
    bar(flip(repulsionTracker),'FaceColor','none','EdgeColor','r','LineWidth',0.25)
    netSE = zeros(gridLength-2,1);
    periodicity = zeros(gridLength-2,1);
    predNumAvuls = zeros(gridLength-2,1); % predicts and plots the number of avulsions expected via subsidence alone
    for anar = 1:gridLength-2
        netSE(gridLength-1-anar) = tfinal*betaSubMtx(anar,1) - (tfinal*depositionMtx(anar,1) + tfinal*diffDepositionPreMtx(anar,1));
    end
    for anar3 = 1:gridLength-2
        predNumAvuls(anar3) = (netSE(anar3)/meandepth);
    end
    
    if sum(predNumAvuls) > (tfinal/checkViableEveryYears)
        predNumAvuls = predNumAvuls * ((tfinal/checkViableEveryYears)/sum(predNumAvuls));
    end
    
    movsumPredNumAvuls = movsum(predNumAvuls,gridLength/(numBins*2));
    plot(movsumPredNumAvuls,'g-')
    title('Downstream avulsion histogram')
    title(['#Avulsions: ',num2str(sum(avulsionLoci(:,1)~=0)), ', Theo. Max:', ...
        num2str(tfinal/checkViableEveryYears),', Ratio: ', num2str(sum(avulsionLoci(:,1)~=0)/(tfinal/checkViableEveryYears))])
    
    xlabel(['distance downdomain (',num2str(cellSize),'m cells)'])
    ylabel(['# Avulsions after ',num2str(tfinal),' years'])
    
    fig7 = figure(7); % this figure shows the generation in which each (active or unhealed/abandoned) cell was last visited
    set(fig7,'Units','normalized');
    set(fig7,'Position',[0 0 0.850 0.90]);
    imagesc(gridspace); % plots the next frame
    colormap jet % sets the colour
    caxis([0 numAvulsions]); % sets the colourbar limits (for #avulsions) at the start
    axis equal % ensures that cells plot as squares, not rectangles
    axis xy
    title('Generation in which each cell was last visited')
    xlabel(['distance (',num2str(cellSize),'m cells)'])
    colorbar % displays the colourbar legend
    drawnow % updates the frame without making a new window
    
    if trackLobeSwitching == 1 % shows the positional history along-strike of the channel over time at 8 distances from the mountain-front
        fig8 = figure(8);
        csvwrite(fullfile(exportDirectory,exportFolder,['lobeTracker',exportName,'.csv']),lobeTracker) % this file can take a hot minute to export
        plot(lobeTracker(:,2),lobeTracker(:,1),'.')
        hold on
        plot(lobeTracker(:,3),lobeTracker(:,1),'.')
        plot(lobeTracker(:,4),lobeTracker(:,1),'.')
        plot(lobeTracker(:,5),lobeTracker(:,1),'.')
        plot(lobeTracker(:,6),lobeTracker(:,1),'.')
        plot(lobeTracker(:,7),lobeTracker(:,1),'.')
        plot(lobeTracker(:,8),lobeTracker(:,1),'.')
        plot(lobeTracker(:,9),lobeTracker(:,1),'.')
        legend('25','50','100','150','200','225','250','275') % in units of cells
        xlim([0 301])
        xlabel('distance across basin (units of cells)')
        ylabel('timestep #')
        title('position of channel at each timestep, at distances (#cells) from the mountain-front')
    end
    
    fig9 = figure(9); % shows relative distance up- or down-stream between successive avulsions
    deltaLoci = zeros(length(finalAvulsionLoci)-1,1);
    for dll = 2:length(deltaLoci)
        deltaLoci(dll) = finalAvulsionLoci(dll+1,3)-finalAvulsionLoci(dll,3);
    end
    % shuffle finalAvulsionLoci's repeatedly to build a bigger list
    shuffledAvulsionLoci = zeros(length(deltaLoci)*10,1); % this is the same loci, repeated, but shuffled, to destroy any autocorrelation
    for sdl = 1:10
        shuffledAvulsionLoci((sdl-1)*length(finalAvulsionLoci(:,3))+1:sdl*length(finalAvulsionLoci(:,3))) = finalAvulsionLoci(randperm(length(finalAvulsionLoci(:,3))),3);
    end
    shuffledDeltaLoci = zeros(length(shuffledAvulsionLoci)-1,1);
    for dls = 2:length(shuffledDeltaLoci)
        shuffledDeltaLoci(dls) = shuffledAvulsionLoci(dls+1)-shuffledAvulsionLoci(dls);
    end
    
    dlh3 = histogram(deltaLoci,'DisplayStyle','stairs','LineStyle','-');
    hold on
    dlh4 = histogram(shuffledDeltaLoci,'DisplayStyle','stairs','LineStyle',':');
    
    dlh3.Normalization = 'probability';
    dlh3.BinWidth = 5;
    dlh4.Normalization = 'probability';
    dlh4.BinWidth = 5;
    title(['skewnesses = ' num2str(skewness(deltaLoci)) ', ' num2str(skewness(shuffledDeltaLoci))])
    legend('observed/modeled','from ~random')
    xlabel('distance downstream of previous avulsion (#cells)')
    ylabel('relative probability')
    
    fig10 = figure(10); % same as Fig9, but local to the avulsion instead of over the whole domain
    narrowDeltaLoci = deltaLoci(deltaLoci<11);
    narrowDeltaLoci = narrowDeltaLoci(narrowDeltaLoci>-11);
    dln3 = histogram(narrowDeltaLoci(ceil(0.00*length(narrowDeltaLoci)+0.1):ceil(0.10*length(narrowDeltaLoci))));
    dln3.Normalization = 'probability';
    hold on
    narrowShuffledDeltaLoci = shuffledDeltaLoci(shuffledDeltaLoci<11);
    narrowShuffledDeltaLoci = narrowShuffledDeltaLoci(narrowShuffledDeltaLoci>-11);
    dln4 = histogram(narrowShuffledDeltaLoci(ceil(0.00*length(narrowShuffledDeltaLoci)+0.1):ceil(0.10*length(narrowShuffledDeltaLoci))));
    dln4.Normalization = 'probability';
    title(['skewness = ' num2str(skewness(narrowDeltaLoci(ceil(0.00*length(narrowDeltaLoci)+0.1):ceil(0.10*length(narrowDeltaLoci)))))])
    xlabel('distance downstream of previous avulsion (#cells)')
    ylabel('relative probability')
    
    fig11 = figure(11); % time between each avulsion
    deltaTimes = zeros(length(finalAvulsionLoci)-1,1);
    for dll = 2:length(deltaLoci)
        deltaTimes(dll) = finalAvulsionLoci(dll+1,1)-finalAvulsionLoci(dll,1);
    end
    
    dlt1 = histogram(deltaTimes(2:ceil(length(deltaTimes)/10)));
    hold on
    dlt3 = histogram(deltaTimes(2:end));
    
    dlt1.BinWidth = 1;
    dlt3.BinWidth = 1;
    
    title(['skewnesses = ' num2str(skewness(deltaTimes(2:ceil(length(deltaTimes)/10)))) ', ' num2str(skewness(deltaTimes(2:end)))])
    legend('first 10%','all time')
    xlabel('time since last avulsion (#timesteps)')
    ylabel('# avulsions')
    
    fig12 = figure(12); % active or abandoned cells per row, with distance from the mountain-front
    cellsPerRow = zeros(1,gridLength);
    for f12 = 1:gridLength
        cellsPerRow(f12) = sum(isChannel(f12,:)~=0);
    end
    plot(cellsPerRow,1:length(cellsPerRow))
    hold on
    plot(movmean(cellsPerRow,10),1:length(cellsPerRow))
    ylim([0,306])
    title('active or abandoned cells/row, with distance from mountain-front')
    ylabel('distance from mountain-front (cells)')
    xlabel('# cells per row')
    
    fig13 = figure(13); % a repeat of Fig4.2, but just the lower panel; makes it easier to export/open in other software
    plot(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,1),movmean(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,2),501),'b:')
    hold on
    plot(finalAvulsionLoci(501:sum(finalAvulsionLoci(:,1)~=0)-500,1),mov05_500(501:sum(finalAvulsionLoci(:,1)~=0)-500,1),'k:')
    plot(finalAvulsionLoci(501:sum(finalAvulsionLoci(:,1)~=0)-500,1),mov95_500(501:sum(finalAvulsionLoci(:,1)~=0)-500,1),'r:')
    plot(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,1),movmean(finalAvulsionLoci(2:sum(finalAvulsionLoci(:,1)~=0)+1,2),5001),'b')
    plot(finalAvulsionLoci(5001:sum(finalAvulsionLoci(:,1)~=0)-5000,1),mov05_5k(5001:sum(finalAvulsionLoci(:,1)~=0)-5000,1),'k')
    plot(finalAvulsionLoci(5001:sum(finalAvulsionLoci(:,1)~=0)-5000,1),mov95_5k(5001:sum(finalAvulsionLoci(:,1)~=0)-5000,1),'r')
    
    title('y-values of avulsion loci vs avulsion #')
    xlabel('timestep #')
    ylabel('y-values of avulsion loci')
    axis([1 timeSteps 0 gridLength])
    
    fig14 = figure(14); % 2D avulsion viability histogram
    imagesc(viableLociTracker)
    colormap(flipud(bone)) % sets the colour
    caxis([0 50]); % sets the colourbar limits (for elevation) at the start
    axis equal % ensures that cells plot as squares, not rectangles
    axis xy
    axis ([0 300 0 300])
    colorbar % displays the colourbar legend
    drawnow % updates the frame without making a new window
    xlabel(['distance (',num2str(cellSize),'m cells)'])
    
    % saving and exporting figures and data
    saveas(fig1,fullfile(exportDirectory,exportFolder,['Fig1_png_',exportName,'.png']))
    saveas(fig1,fullfile(exportDirectory,exportFolder,['Fig1_fig_',exportName,'.fig']))
    saveas(fig1,fullfile(exportDirectory,exportFolder,['Fig1_pdf_',exportName,'.pdf']))
    
    saveas(fig2,fullfile(exportDirectory,exportFolder,['Fig2_png_',exportName,'.png']))
    saveas(fig2,fullfile(exportDirectory,exportFolder,['Fig2_fig_',exportName,'.fig']))
    
    saveas(fig3,fullfile(exportDirectory,exportFolder,['Fig3_png_',exportName,'.png']))
    saveas(fig3,fullfile(exportDirectory,exportFolder,['Fig3_fig_',exportName,'.fig']))
    
    saveas(fig4,fullfile(exportDirectory,exportFolder,['Fig4_png_',exportName,'.png']))
    saveas(fig4,fullfile(exportDirectory,exportFolder,['Fig4_fig_',exportName,'.fig']))
    
    saveas(fig5,fullfile(exportDirectory,exportFolder,['Fig5_png_',exportName,'.png']))
    saveas(fig5,fullfile(exportDirectory,exportFolder,['Fig5_fig_',exportName,'.fig']))
    
    saveas(fig6,fullfile(exportDirectory,exportFolder,['Fig6_png_',exportName,'.png']))
    saveas(fig6,fullfile(exportDirectory,exportFolder,['Fig6_fig_',exportName,'.fig']))
    saveas(fig6,fullfile(exportDirectory,exportFolder,['Fig6_pdf_',exportName,'.pdf']))
    saveas(fig6,fullfile(exportDirectory,exportFolder,['Fig6_svg_',exportName,'.svg']))
    
    saveas(fig7,fullfile(exportDirectory,exportFolder,['Fig7_png_',exportName,'.png']))
    saveas(fig7,fullfile(exportDirectory,exportFolder,['Fig7_fig_',exportName,'.fig']))
    
    csvwrite(fullfile(exportDirectory,exportFolder,['finalGridspace',exportName,'.csv']),gridspace)
    
    if trackLobeSwitching == 1
        saveas(fig8,fullfile(exportDirectory,exportFolder,['Fig8_png_',exportName,'.png']))
        saveas(fig8,fullfile(exportDirectory,exportFolder,['Fig8_fig_',exportName,'.fig']))
    end
    
    if saveAvulsionLoci == 1 % this file can take a hot minute to export
        csvwrite(fullfile(exportDirectory,exportFolder,['avulsionLoci',exportName,'.csv']),avulsionLoci)
    end
    
    saveas(fig9,fullfile(exportDirectory,exportFolder,['Fig9_png_',exportName,'.png']))
    saveas(fig9,fullfile(exportDirectory,exportFolder,['Fig9_fig_',exportName,'.fig']))
    saveas(fig9,fullfile(exportDirectory,exportFolder,['Fig9_svg_',exportName,'.svg']))
    saveas(fig9,fullfile(exportDirectory,exportFolder,['Fig9_pdf_',exportName,'.pdf']))
    
    saveas(fig10,fullfile(exportDirectory,exportFolder,['Fig10_png_',exportName,'.png']))
    saveas(fig10,fullfile(exportDirectory,exportFolder,['Fig10_fig_',exportName,'.fig']))
    saveas(fig10,fullfile(exportDirectory,exportFolder,['Fig10_svg_',exportName,'.svg']))
    saveas(fig10,fullfile(exportDirectory,exportFolder,['Fig10_pdf_',exportName,'.pdf']))
    
    saveas(fig11,fullfile(exportDirectory,exportFolder,['Fig11_png_',exportName,'.png']))
    saveas(fig11,fullfile(exportDirectory,exportFolder,['Fig11_fig_',exportName,'.fig']))
    saveas(fig11,fullfile(exportDirectory,exportFolder,['Fig11_svg_',exportName,'.svg']))
    saveas(fig11,fullfile(exportDirectory,exportFolder,['Fig11_pdf_',exportName,'.pdf']))
    
    saveas(fig12,fullfile(exportDirectory,exportFolder,['Fig12_png_',exportName,'.png']))
    saveas(fig12,fullfile(exportDirectory,exportFolder,['Fig12_fig_',exportName,'.fig']))
    
    saveas(fig13,fullfile(exportDirectory,exportFolder,['Fig13_png_',exportName,'.png']))
    saveas(fig13,fullfile(exportDirectory,exportFolder,['Fig13_fig_',exportName,'.fig']))
    saveas(fig13,fullfile(exportDirectory,exportFolder,['Fig13_pdf_',exportName,'.pdf']))
    
    saveas(fig14,fullfile(exportDirectory,exportFolder,['Fig14_png_',exportName,'.png']))
    saveas(fig14,fullfile(exportDirectory,exportFolder,['Fig14_fig_',exportName,'.fig']))
    saveas(fig14,fullfile(exportDirectory,exportFolder,['Fig14_pdf_',exportName,'.pdf']))
    
end