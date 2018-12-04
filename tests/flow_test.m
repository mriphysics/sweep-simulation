function [dat, tissue, RF, motion] = flow_test()
%% Laurence Jackson, BME, KCL, 2018
% 
% script to simulate blood flow contrast when using the SWEEP method
% 
% output:
%   dat(:,1) = simulation results for tissue component
%   dat(:,2) = simulation results for blood flow component
% 

saveloc = 'simresults/flowsim.mat';

%% set up parameter objects
load('tests/settings_flow.mat'); % load bulk of simulation parameters

if(~exist('simresults','dir')); mkdir simresults; end

%% sim
Rf.npulses = 90*5;
Rs = linspace(0,1,11);
fs = linspace(-40e-3,40e-3,9); % TODO: run with 5 flowrates

for ii = 1:length(Rs)
    for ff = 1:length(fs)
        RF.swp = Rs(ii);
        RF.flip = 70;
        
        motion.flow = 0e-3;
        tissue.T1 = 1820; tissue.T2 = 99;
        [dat{ii,ff,1}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % tissue
        
        motion.flow = fs(ff);
        tissue.T1 = 1550; tissue.T2 = 275;
        [dat{ii,ff,2}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % flow
        close all; % suppress image outputs
    end
end

disp(['Saving simulation results at ' saveloc])
save(saveloc,'dat','-v7.3')


end