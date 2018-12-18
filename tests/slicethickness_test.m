function [dat, tissue, RF, motion] = slicethickness_test()
%% Laurence Jackson, BME, KCL, 2018
% 
% Script to investigate the slice thickness of sweep versus sweep rate

%% simulation
if(~exist('simresults','dir')); mkdir simresults; end

Rs = linspace(0,1,11);
load('tests/settings_thickenss_bSSFP.mat');
RF.npe = 1e12;
for ii = 1:length(Rs)
    RF.swp = Rs(ii);
    [dat{ii}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion);
    close all
end
save('simresults/bssfp_Rs01_nomotion.mat','dat')

clearvars -except Rs

load('tests/settings_thickenss_SPGR.mat')
RF.npe = 1e12;
for ii = 1:length(Rs)
    RF.swp = Rs(ii);
    [dat{ii}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion);
    close all
end
save('simresults/SPGR_Rs01_nomotion.mat','dat')

%% figures
