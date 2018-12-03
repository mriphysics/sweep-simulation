function [dat, tissue, RF, motion] = motion_test()
%% Laurence Jackson, BME, KCL, 2018
% 
% script to simulate the effects of respiratory motion of the pulse
% excitation profile with SWEEP
% 
% output:
%   dat = simulation results for tissue component
% 

%% load paramters
load('tests/settings_motion.mat'); % load bulk of simulation parameters

%% sim
Rs = linspace(0,1,11);

for ii = 1:length(Rs)
    RF.swp = Rs(ii);
    
    RF.seq = 'bssfp';
    [dat{ii}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % bffe
    save('simresults/motionsim_bffe.mat','dat')
    
    RF.seq = 150;
    [dat{ii,2}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % SPGR
    save('simresults/motionsim_spgr.mat','dat')
end

end