function [dat, tissue, RF, motion] = motion_test()
%% Laurence Jackson, BME, KCL, 2018
% 
% script to simulate the effects of respiratory motion of the pulse
% excitation profile with SWEEP
% 
% output:
%   dat = simulation results for tissue component
% 

clear

%% load paramters
load('tests/settings_motion.mat'); % load bulk of simulation parameters

%% sim
Rs = linspace(0,1,11);
zero_resp_phase = 500;
resp = [zeros([1,zero_resp_phase]), sin(linspace(0,(2*pi*(RF.npulses-zero_resp_phase).*RF.TR*0.001*(motion.respfreq)), RF.npulses - zero_resp_phase)).*(motion.respmag./2)]; 
motion.custom = resp;
plot(motion.custom)

% bssfp
RF.seq = 'bssfp';
for ii = 1:length(Rs)
    
    RF.swp = Rs(ii);
    [dat{ii}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % bffe

end
save('simresults/motionsim_bffe.mat','dat')

% spgr
RF.seq = 150;
RF.TR = 15;
RF.flip = 10;
for ii = 1:length(Rs)
    
    RF.swp = Rs(ii);
    [dat{ii}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % bffe

end
save('simresults/motionsim_spgr.mat','dat')

end