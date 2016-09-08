% Function that loads all the parameters. To explore variability of one of
% the parameters just remove it from here and change the program to convert
% it in input of "solspace" and "beatmodes".

function parameters(varargin)
%Declare global variables and parser
global Sp D0 chib kp k xi xitrn xirot
global bc motor
p = inputParser;

%INPUT PARAMETERS
% Default parameters
defLength            = 12;        % in um
defFrequency         = 30;        % in 1/s
defBasalStiffness    = 0.;        % in pN/um
defBasalFriction     = 0.;        % in pN*s/um
defHeadFrictionRot   = 0.;        % in pN*s*um
defHeadFrictionTrans = 0.;        % in pN*s/um
defSlidingStiffness  = 0.;        % in pN/um^2
defSlidingFricction  = 0.;        % in pN*s/um^2
defPivotStiffness    = 0.;        % in pN
defAsymmetry         = 0.0;       % in rad/um
defRadius            = 0.06;      % in um
defBendingRigidity   = 400;       % in pN*um^2
defNormalFriction    = 0.0034;    % in pN*s/um^2
defMotor             = 'sliding';
defBoundaries        = 'free-free';

% Available boundary conditions and motor models
availBoundaries = {'free-free','clamped-free','swim-free','pivot-free'};
availMotors     = {'sliding','curvature','dyn-curvature'};

%Parse inputs
addOptional(p,'Length',defLength,@isnumeric);
addOptional(p,'Frequency',defFrequency,@isnumeric);
addOptional(p,'BasalStiffness',defBasalStiffness,@isnumeric);
addOptional(p,'BasalFriction',defBasalFriction,@isnumeric);
addOptional(p,'HeadFrictionRot',defHeadFrictionRot,@isnumeric);
addOptional(p,'HeadFrictionTrans',defHeadFrictionTrans,@isnumeric);
addOptional(p,'SlidingStiffness',defSlidingStiffness,@isnumeric);
addOptional(p,'SlidingFriction',defSlidingFricction,@isnumeric);
addOptional(p,'PivotStiffness',defPivotStiffness,@isnumeric);
addOptional(p,'Asymmetry',defAsymmetry,@isnumeric);
addOptional(p,'Radius',defRadius,@isnumeric);
addOptional(p,'BendingRigidity',defBendingRigidity,@isnumeric);
addOptional(p,'NormalFriction',defNormalFriction,@isnumeric);
addOptional(p,'Motor',defMotor,@(x) any(validatestring(x,availMotors)));
addOptional(p,'Boundaries',defBoundaries,...
              @(x) any(validatestring(x,availBoundaries)));
parse(p,varargin{:});

%% GLOBAL PARAMETERS
% Calculate dimensionless and other global variables
bc    = p.Results.Boundaries;
motor = p.Results.Motor;
Sp    = 2*pi*p.Results.Frequency*p.Results.NormalFriction...
        *p.Results.Length^4/p.Results.BendingRigidity;
D0    = p.Results.Asymmetry*p.Results.Length;
chib  = (p.Results.BasalStiffness...
        +1i*2*pi*p.Results.Frequency*p.Results.BasalFriction)...
        *(p.Results.Length*p.Results.Radius^2/p.Results.BendingRigidity);
k     = p.Results.SlidingStiffness...
        *(p.Results.Length^2*p.Results.Radius^2/p.Results.BendingRigidity);
xi    = p.Results.SlidingFriction...
        *(p.Results.Length^2*p.Results.Radius^2/p.Results.BendingRigidity);
kp    = p.Results.PivotStiffness...
        *(p.Results.Length*p.Results.Radius^2/p.Results.BendingRigidity);
xitrn = 1i*2*pi*p.Results.Frequency*p.Results.HeadFrictionRot...
        *p.Results.Length^3/p.Results.BendingRigidity;
xirot = 1i*2*pi*p.Results.Frequency*p.Results.HeadFrictionTrans...
        *p.Results.Length/p.Results.BendingRigidity;
