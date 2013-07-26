function [tY, op, d1MCCt] = SIMMpp(p_laser, p_cutX, k_tens, n)
% This function sovles a system on ODEs describing the molecular regulatory network
% that underlies the spindle assembly checkpoint (SAC). 
% The model is based on a model published by He et al. 2011 in PNAS.
% This implementation of the model is designed to mimick a set of
% experiments, which look at the dynamics of checkpoint re-engagement in
% response to a stimulus (chromosome detachment from mitotic spindle).
% The stimulus is modelled as a sudden reduction in chromosome attachment
% (Xtens) 

%Give the user some feedback!
fprintf('------------------------------------------------------------------\n');
fprintf('\tRUNNING CHECKPOINT RE-ENGAGEMENT SIMULATION #%d \n', n);
fprintf('\t\t p_plaser = %d\n', p_laser);
fprintf('\t\t p_cutX = %d\n', p_cutX);
fprintf('\t\t k_tens = %d\n', k_tens);
fprintf('------------------------------------------------------------------\n\n');

verbo=1;

%Define the vector holding the initial conditions
% CycB = Y(1);
% kTa = Y(2);
% MCCt = Y(3);
% APCMCC = Y(4);
% MCCu = Y(5);
% Sec = Y(6);
% Xtens = Y(7);
% CAPP = Y(8);
initialConditions=[0.47, 0.19, 1.99, 0.98, 0.09, 1.8, 1, 0.2];

%Definition of parameters
%p_laser = 10;
%p_cutX = 0.7;
p_Mad2t = 2;
p_CAPPt = 1;
p_APCt = 1;
p_Xt = 1;
%k_tens = 0.02;
k_aNt = 5;
k_u = 10;
k_du = 5;
k_cat = 0.5;
k_as = 100;
k_di = 0.5
k_scyc = 0.01;
k_dcyc = 0.01
k_dcycapc = 0.5;
k_an = 1;
k_in = 2
k_ssec = 0.1;
k_dsec = 0.05;
k_dsecapc = 0.5;
k_app = 0.05;
k_ipp = 5;



%Pack parameters into a vector
parameters=[p_laser, p_cutX, p_Mad2t, p_CAPPt, p_APCt, p_Xt, k_tens, k_aNt, k_u, k_du, k_cat, k_as, k_di, k_scyc, k_dcyc, k_dcycapc, k_an, k_in, k_ssec, k_dsec, k_dsecapc, k_app, k_ipp];

%Call the integration routine
[tY, op, d1MCCt] = integration(verbo,parameters, initialConditions, n);

%=========================================================
% Integration
%=========================================================

function [tY, op, d1MCCt] = integration(v,parameters,initialConditions, n)

if v==1
    fprintf('Setting up...\n');
end
tic
%trans defines the thransient integration time, allow system to reach ss
%before actually starting the integration
trans=0;

%total time in minutes
tend=500;
%time step
tstep=0.1;

fprintf('\t\t');
toc

%call the nested integration routine run
[t1,t2,Y1,Y2] = run(initialConditions, parameters ,trans, tend, tstep);


fprintf('\nSummarising results...\n');
tic
%Merge the arays holding time and concentration values from pre-cutting
%and after cutting into a single array.
t = cat(1,t1,t2);
Y = cat(1,Y1,Y2);
tY=[t,Y];             

%Find the rate of change of MCCt
d1MCCt = diff(Y(:,3));
%Add a zero at the end to give it the same dimension as the data array
d1MCCt = cat(1,d1MCCt,0);

%Define a time window around the laser cutting, for the local min and max
%to be compared
stimulusID = round(parameters(1) * 10);

if parameters(1) < 10
    lowerID = 1;
else
    lowerID = stimulusID - 99;
end

upperID = stimulusID + 99;
lD_d1MCCt = d1MCCt(lowerID:upperID);
aD_d1MCCt = d1MCCt(stimulusID :upperID);

fprintf('max(d1MCCt) in stimulus domain (+- 99) = %d\n', abs(max(lD_d1MCCt)));
fprintf('mean(d1MCCt) after stimulus (:stimulus:+99dt) = %d\n\t\t', mean(aD_d1MCCt));
%%Use the comparison on local min & max around laser cutting to evaluate,
%wheter the checkpoin is responsive to re-engagement.
if (abs(max(lD_d1MCCt)) + abs(min(lD_d1MCCt))) >  abs(min(lD_d1MCCt)) && max(Y((stimulusID:upperID),3)) >= 0.6 && mean(aD_d1MCCt) > - 0.0001
    %Opinion is a string printed to the screen
    opinion = sprintf('Checkpoint is RESPONSIVE to re-engagement');
    %op is a logical variable used by the higher level sweeping script
    op = 1;

else
    opinion = sprintf('Checkpoint is NOT RESPONSIVE to re-engagement');
    op = 0;
end
toc

fprintf('\nOPINION:\n');
fprintf(opinion);
fprintf('\n------------------------------------------------------------------\n');

% % ============================================================================================
% % Visualisation
% % ============================================================================================
% set(figure(1),'Position', [50 100 600 400]);
% 
% clf;
% 
% % CycB = Y(1);
% % kTa = Y(2);
% % MCCt = Y(3);
% % APCMCC = Y(4);
% % MCCu = Y(5);
% % Sec = Y(6);
% % Xtens = Y(7);
% 
% plotname = sprintf('130627 He Reengagement Responsiveness #%d', n);
% figure(1); 
% plot(t,Y(:,1),'k',t,Y(:,3),'g',t,Y(:,4),'b',t,Y(:,6),'c',t,Y(:,7),'r'); 
% xlabel('Time (min)');
% ylabel('Concentrations');
% xlim([0 100]);
% ylim([0 2]);
% title(plotname,... 
%   'FontWeight','bold')
% set(gca,'xtick',[0:10:tend]);
% h = legend('CycB','MCCt','APCMCC', 'Sec', 'Xtens', 2);
% 
% figure(2);
% plot(t,d1MCCt);
% xlim([0 100]);
% ylim([-0.05 0.05])





% ============================================================================================
% Run
% ============================================================================================

function [t1,t2,Y1,Y2]=run(initialConditions,parameters,trans,tend,tstep)

fprintf('\nInitialising calculations...\n');
tic
ttrans = [0:tstep:trans];
tspan1 = [0:tstep:parameters(1)];
tspan2 = [0:tstep:(tend-tspan1(end))];

option = odeset('RelTol', 1e-5);


if trans > 0 
    [t1, Y1] = ode23s(@dYdt,ttrans,initialConditions,option,parameters);
    x0=Y(end,:);
end
fprintf('\t\t');
toc

fprintf('\nCheckpoint disengagement...\n');
tic
[t1, Y1] = ode23s(@dYdt,tspan1,initialConditions,option,parameters);

initialConditions2 = Y1(end,:);
initialConditions2(7) = parameters(2);
fprintf('\t\t');
toc

fprintf('\nLaser cutting after %d s...\n', t1(end));
fprintf('New initial conditions:\n CycB:\t\t %d \n kTa:\t\t %d \n MCCt:\t\t %d \n APCMCC:\t %d \n CAPP:\t\t %d \n Sec:\t\t %d \n Xtens:\t\t %d \n', initialConditions2);
tic
[t2, Y2] = ode15s(@dYdt,tspan2,initialConditions2,option,parameters);
fprintf('\t\t');
toc
t2 = t2 + t1(end) + tstep;


% ============================================================================================
% dxdt
% ============================================================================================

function dYdt = dYdt(t,Y,parameters)
CycB = Y(1);
kTa = Y(2);
MCCt = Y(3);
APCMCC = Y(4);
MCCu = Y(5);
Sec = Y(6);
Xtens = Y(7);
CAPP = Y(8);

%Map the model parameters to the elements of the parameter vector
p_laser = parameters(1);
p_cutX = parameters(2);
p_Mad2t = parameters(3);
p_CAPPt = parameters(4);
p_APCt = parameters(5);
p_Xt = parameters(6);
k_tens = parameters(7);
k_aNt = parameters(8);
k_u = parameters(9);
k_du = parameters(10);
k_cat = parameters(11);
k_as = parameters(12);
k_di = parameters(13);
k_scyc = parameters(14);
k_dcyc = parameters(15);
k_dcycapc = parameters(16);
k_an = parameters(17);
k_in = parameters(18);
k_ssec = parameters(19);
k_dsec = parameters(20);
k_dsecapc = parameters(21);
k_app = parameters(22);
k_ipp = parameters(23);

APC = p_APCt - APCMCC;

%Function for APC

%Define the model ODEs
dMCCt = k_aNt * kTa * ( p_Mad2t - MCCt)/(0.1 + p_Mad2t - MCCt) - k_u * APC * MCCu;
dMCCu = k_cat * APCMCC - ( k_du + k_u * APC) * MCCu;
dSec = k_ssec - (k_dsec + k_dsecapc * APC) * Sec;
dAPCMCC = k_as * ( p_APCt - APCMCC ) * ( MCCt - APCMCC - MCCu ) - ( k_di + k_cat ) * APCMCC;
dCycB = k_scyc - ( k_dcyc + k_dcycapc * APC ) * CycB;
dkTa = k_an * CycB * ( 1 - Xtens - kTa ) - k_in * p_CAPPt * kTa;
dXtens = k_tens * heaviside(1-Xtens);
dCAPP = k_app * ( p_CAPPt - CAPP) - k_ipp * CycB * CAPP;


%Map 
dYdt = [dCycB; dkTa; dMCCt; dAPCMCC; dMCCu; dSec; dXtens; dCAPP];

