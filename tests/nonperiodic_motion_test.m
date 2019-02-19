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

% define custom non-periodic motion
motion.custom = zeros(1,RF.npulses);
nskips = 8;
skiploc = floor(linspace(400, 1400,nskips));
skipamp = linspace(1e-3,100e-3,nskips);
skipdur = 50;
for ii = 1:nskips
    x = skiploc(ii)-skipdur:skiploc(ii)+skipdur;
    y = normpdf(x,skiploc(ii));
    y = ((y- min(y)) / ((max(y) - min(y)))) .* skipamp(ii);
    motion.custom(x) = y;
end

plot(motion.custom)

for ii = 1:length(Rs)
    RF.swp = Rs(ii);
    
    RF.seq = 'bssfp';
    [dat{ii}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % bffe
    save('simresults/motion_nonperiodic_sim_bffe.mat','dat')
    
%     RF.seq = 150;
%     [dat{ii,2}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % SPGR
%     save('simresults/motionsim_spgr.mat','dat')
end

end