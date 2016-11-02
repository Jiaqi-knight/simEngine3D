function state = dynamicsAnalysis(sys,timeStart,timeEnd,timeStep) % perform dynamics analysis
% perform dynamics analysis on the system. System must be not be
% fully constrained (nDOF > 0).
% 
% dynamics analysis peformed using the 6 steps outline on slide 43 of 
% lecture 10/17/16, ME751_f2016
% 
% inputs: measured in seconds.
%   timeStart : the starting time of simulation, usually 0
%   timeEnd   : ending time of simulation
%   timeStep  : step size of time during stimulation (optional)
%
% output:
%   state = cell array of states of system

% if sys.nDOF == 0
%     error('For  dynamics analysis, the total number of degrees of freedom must be greater than zero.');
% end
if ~exist('timeStep','var') || isempty(timeStep)
    timeStep = 10^-2; % default time step
end

timeGrid = timeStart:timeStep:timeEnd; % establish time grid
state = cell(length(timeGrid),1); % preallocate for saved state


%%%%%%%%%%%
% check initial conditions (ME751_f2016, slide 41, lecture 10/17/16)
iT = 1; % at timeStart
sys.setSystemTime(timeGrid(iT)); % set system time 
sys.checkInitialConditions();
sys.findInitialAccelerations();

% store initial system state: (iT = 1)
state{1} = sys.getSystemState();

%%%%%%%%%%% GOOD WITH ABOVE
% using initial conditions, 
% find first time step using Quasi-Newton method and BDF of order 1
iT = 2; % first integration step
sys.setSystemTime(timeGrid(iT)); % set system time 
sys.QN_BDF1(timeStep,state{1}); % send initial state
state{2} = sys.getSystemState();

return

% iterate throughout the time grid:
for iT = 3:length(timeGrid)
    t = timeGrid(iT); % current time step
    sys.setSystemTime(t); % set system time 
    
    
    disp(['Dynamics analysis completed for t = ' num2str(t) ' sec.']);
end

end %dynamicsAnalysis