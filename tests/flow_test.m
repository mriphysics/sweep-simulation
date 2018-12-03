function [dat, tissue, RF, motion] = flow_test()
%% Laurence Jackson, BME, KCL, 2018
% 
% script to simulate blood flow contrast when using the SWEEP method
% 
% output:
%   dat(:,1) = simulation results for tissue component
%   dat(:,2) = simulation results for blood flow component
% 

%% set up tissue parameter object
load('tests/settings_flow.mat'); % load bulk of simulation parameters

if(~exist('simresults','dir')); mkdir simresults; end
    
%% sim
Rs = linspace(0,1,11);
for ii = 1:length(Rs)
    RF.swp = Rs(ii);
    RF.flip = 70;
    
    motion.flow = 0e-3;
    tissue.T1 = 1820; tissue.T2 = 99;    
    [dat{ii,1}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % tissue
    
    motion.flow = -40e-3;
    tissue.T1 = 1550; tissue.T2 = 275;    
    [dat{ii,2}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % flow

end

save('simresults/flowsim_negative.mat','dat')

end