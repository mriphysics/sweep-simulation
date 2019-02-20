function [dat, tissue, RF, motion] = nonperiodic_motion_test()
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
motion.custom = zeros(1, RF.npulses);
nskips = 8;
skiploc = floor(linspace(400, 1400,nskips));

maxamp = 80e-3;
skipamp = exp(0.9*(1:nskips));
skipamp = ((skipamp- min(skipamp)) / ((max(skipamp) - min(skipamp)))) .* maxamp;

skipdur = 25;

for ii = 1:nskips
    x = (skiploc(ii)-skipdur):(skiploc(ii)+skipdur);
    y = normpdf(x,skiploc(ii), skipdur/2);
    y = ((y- min(y)) / ((max(y) - min(y)))) .* skipamp(ii);
    motion.custom(x) = y;
end
plot(motion.custom)

 RF.seq = 'bssfp';
for ii = 1:length(Rs)
    
    RF.swp = Rs(ii);
    [dat{ii}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % bffe

end

save('simresults/motionsim_nonperiodic_bffe.mat','dat','-v7.3')

end