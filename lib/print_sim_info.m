function print_sim_info(tissue, RF, motion)
%% Laurence Jackson, BME, KCL, 2018
% 
% Script to print info about simulation objects to console
% Looks tidiers hidden in different function
%% print commands
fprintf('\n\t::SIMULATION SUMMARY::\t\n')
fprintf('Tissue paramters:\n')
fprintf(' \tT1 \t\t\t\t%d\n',tissue.T1)
fprintf(' \tT2 \t\t\t\t%d\n',tissue.T2)

fprintf('RF paramters:\n')
fprintf(' \tSweep rate: \t%.2f%%\n',RF.swp)
fprintf(' \tNumber RFs: \t%d\n',RF.npulses)
fprintf(' \tFlip angle: \t%.2f\n',RF.flip)
fprintf(' \tSlice thk: \t\t%.2fmm\n',RF.thk*1e3)

fprintf('Motion paramters:\n')
fprintf(' \tResp Freq: \t\t%.2fHz\n',motion.respfreq);
fprintf(' \tResp Mag: \t\t%.2fmm\n',motion.respmag*1e3);
fprintf(' \tFlow:  \t\t\t%.2fmm/s\n',motion.flow*1e3)

fprintf('Sequence definition:\n')
if RF.seqspec == 1
    fprintf(' \tseq order: \t\t%s\n',RF.dynorder);
    fprintf(' \tNum slices: \t%d\n',RF.nslice);
    fprintf(' \tNum dynamics \t%d\n',RF.ndyn);
end


end