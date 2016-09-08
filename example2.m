%% INFO
% This is an example file on how to calculate flagellar beats. This is done
% in three steps: (1) choose parameters, (2) explore phase space, (3)
% calculate solutions.
clear all;

%% CHOOSE PARAMETERS
% Default parameters are defined in "parameters.m". We can specify as many
% as desired using a structure. We give some Sliding Control (SC), one
% Curvature Control (CC) and one Dynamic Curvature Control (DCC) examples.

%% SC examples
% Example 1, standing waves: using the default symmetric bc (free-free) and
% the default basal compliance (null) the problem is symmetric and only
% standing waves arise.
% inputStanding = struct('Frequency',10,...
%                        'Length',60,...
%                        'BendingRigidity',1700,...
%                        'Motor','sliding');
                   
% Example 2, freely swimming flagella(Fig 5.2C1 of my thesis): by simply
% adding a basal compliance polar symmetry is broken, and traveling modes
% appear (the first travels forward, the second backwards).
% inputFreeSliding = struct('Frequency',10,...
%                           'Length',60,...
%                           'BendingRigidity',1700,...
%                           'Motor','sliding',...
%                           'BasalStiffness',50000.,...
%                           'BasalFriction',400,...
%                           'Boundaries','free-free'); 

% Example 3, clamped head forward traveling (HFSP parameters): for a
% clamped base only by carefully tuning the basal compliance we obtain a
% forward traveling mode.
inputSperm = struct('Frequency',21,...
                    'Length',58,...
                    'BendingRigidity',10.,...
                    'Motor','sliding',...
                    'BasalStiffness',10.,...
                    'BasalFriction',1274.,...
                    'Boundaries','clamped-free');

%% CC examples
% This example corresponds to typical curvature control beats
% inputCurvature = struct('Motor','curvature',...
%                         'SlidingStiffness',15000);
                    
%% DCC example
% This example corresponds to the Chlamy parameters from eLife paper
% inputChlamyDynamicCurvature = struct('Frequency',50,...
%                              'Length',12,...
%                              'Asymmetry',-0.25,...
%                              'BasalStiffness',50000.,...
%                              'BendingRigidity',400,...
%                              'Motor','dyn-curvature');
%% EXPLORE PHASE-SPACE
% Load parameters. For defaults, use "parameters" with no arguments
parameters(inputSperm);

% For sliding control, we explore the [chi1,chi2] space
chi1=-200:1:0;
chi2=-150:.25:0;
space = solspace(chi1,chi2);

% For curvature control, 
% beta1=-2:.1:2;
% beta2=-16:.1:16;
% space = solspace(beta1,beta2);

% For dynamic curvature control, which fails to properly find minima
% chi1=0:.1:40;     % refinement chi1=36:.005:40;
% beta2=-15:.01:0;  % refinement beta2=-14:.005:-13.5;
% space = solspace(chi1,beta2);

% Since the data is stored, we can plot it as wished
figure;
h=surfc(space.xyz{1},space.xyz{2},space.xyz{3}');
set(h,'LineStyle','none');

%% CALCULATE SOLUTIONS
% Taking the seeds from the phase space we can calculate the modes
solutions = beatmodes(space.seeds);

% If we have a good guess of the seed, we can just pass it to "beatmodes"
% for example for dynamic curvature control we can choose as seed 
% mysolution = beatmodes([20;-10]);

%% STUDY A SOLUTIONS
% We can now pick a solution and check the error, which should be small
% s1 = solutions(1);
% disp(s1.err);
% 
% % We can also plot the solution over arc-length "s"
% s = 0:0.01:1;
% psi1 = s1.A(1)*exp(s1.k(1)*s)+s1.A(2)*exp(s1.k(2)*s)...
%      + s1.A(3)*exp(s1.k(3)*s)+s1.A(4)*exp(s1.k(4)*s);
% figure; plot(s,abs(psi1));
% figure; plot(s,unwrap(angle(psi1)));