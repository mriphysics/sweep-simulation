%% Sweep_sim_main
% Laurence Jackson, BME, KCL, 2018
%
% Script for running sweep simulation
% 
% Usage:
%   
% 
% Laurence.Jackson@kcl.ac.uk 

%% 
clear 

addpath(genpath('lib'))
addpath(genpath('bin'))
addpath(genpath('tests'))
addpath(genpath('remove'))


%% set up tissue parameter object
tissue.T1 = 1412; % simulated T1
tissue.T2 = 50; % simulated T2

%% set up RF parameter object
RF.npulses = 600 ; % number of pulses in simulation
RF.npe = 90; % number of phase encodes per slice

% Basic sequence parameters
RF.seq = 150; % phase cycling method string::<'bssfp'> or a double::<angle_in_degrees> phase cycling angle for SPGR
RF.swp = 0.0; % sweep rate as percentage of slice thickness moved per TR
RF.thk = 4.0e-3; % nominal slice thickness
RF.slicegap = 0.0*1e-3; % slice gap (used if RF.swp == 0)

RF.block=0;

if ischar(RF.seq) % balanced sequence
    RF.TR = 6;
    RF.flip = 22;
else % non-balanced sequence
    RF.TR = 15;
    RF.flip = 10;
end

% Sequence order
RF.seqspec = 0; % specify sequence - enabaling this option will override RF.npulses and use RF.npe to simulate a given number of slices and dynamics
RF.dynorder = 'slices'; % Interleave order of npe encoding ['slices' ; 'dynamics']
RF.sliceorder = 'odd-even'; % Interleave order of slices ['ascending' ; 'descending' ; 'odd-even' ; 'random']
RF.ndyn = 1;
RF.nslice = 4;

RF.catalysation = []; % vector of catalysation pulses to apply 


%% set up motion paramter object
motion.flow = 0e-3; % [m/s] blood flow (must be +ve at the moment)
motion.respfreq = 0.3; % resp frequency in Hz
motion.respmag = 2e-3; % resp magnitude in mm (note: this is magnitude of motion in the through-plane direction)

%% Tests
% loads a mat file with preset RF, tissue and motion settings uncomment the test to run.

% load('static_test');
% load('respiration_test');

%% Other - add loop/intercept individual parameters here to explore any specific setting
Rs = linspace(0,1,11);
for ii = 1:length(Rs)
    RF.swp = Rs(ii);
    RF.flip = 70;
    
    motion.flow = 0e-3;
    tissue.T1 = 1820; tissue.T2 = 99;    
    [dat{ii,1}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % tissue
    
    motion.flow = 40e-3;
    tissue.T1 = 1550; tissue.T2 = 275;    
    [dat{ii,2}, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion); % flow

end
save('simresults/flowsim.mat','dat')

%% default sim
[dat, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion);

%% === Results =============================================================================
FontSize1 = 12;
FontSize2 = 16;

% Display flipmat
figure();imagesc(dat.tissue.vec.*1000-dat.RF.range,1:RF.npulses,rad2deg(dat.flipmat)); hold on;
colormap 'parula'; caxis([0 RF.flip+10]); c = colorbar; c.Label.String = 'Flip angle';set(gca,'fontsize', FontSize1);
ylabel('RF pulse #', 'FontSize', FontSize2);xlabel('z location (mm)', 'FontSize', FontSize2);hold off;

% Display s0
figure();imagesc(dat.tissue.vec.*1000-RF.range,1:RF.npulses,abs(dat.s0)); hold on;
colormap 'jet'; c = colorbar; c.Label.String = 'S/S_0';set(gca,'fontsize', FontSize1);
ylabel('RF pulse #', 'FontSize', FontSize2);xlabel('z location (mm)', 'FontSize', FontSize2);hold off;

% Plot total signal
for pls = 1:(RF.npulses)
    dat.sig(pls) = trapz(abs(dat.s0(pls,:)));
end
figure();box on; set(gca,'fontsize', FontSize1); ylabel('Total signal(A.U)', 'FontSize', FontSize2);xlabel('RF pulse number', 'FontSize', FontSize2); hold on;
plot(1:RF.npulses,dat.sig,'linewidth', 1.5); 
xlim([0 RF.npulses]);
kcs = (linspace((RF.npe/2),(RF.nslice*RF.ndyn*RF.npe)-(RF.npe/2),RF.nslice*RF.ndyn)); % centres of kspace
plot(kcs,dat.sig(kcs),'o','MarkerFaceColor', 'r','MarkerEdgeColor', 'k'); 
hold off;

% Display s0 profiles
sz = size(dat.s0);
figure();box on; set(gca,'fontsize', FontSize1); ylabel('S/S0', 'FontSize', FontSize2);xlabel('z location (mm)', 'FontSize', FontSize2);xlim([-RF.range (max(tissue.vec)*1e3)-RF.range]);hold on;
for tt = 1:1:sz(1)
    cla;
    plot( dat.tissue.vec.*1000-RF.range,abs(dat.s0(tt,:)));ylim([0 0.3]);hold on
%     plot( dat.tissue.vec.*1000-RF.range,dat.flipmat(tt,:));ylim([0 0.3]);hold on
    drawnow;
    set(gcf,'Color',[1 1 1]); % make background white
    f = getframe(gcf);
    imgif{tt} = frame2im(f);
end

% Display s0 profiles
sz = size(dat1.s0);
figure();box on; set(gca,'fontsize', FontSize1); ylabel('S/S0', 'FontSize', FontSize2);xlabel('z location (mm)', 'FontSize', FontSize2);xlim([-RF1.range (max(tissue1.vec)*1e3)-RF1.range]);hold on;
for tt = 1:1:sz(1)
    cla;
    plot( dat1.tissue.vec.*1000-RF1.range,abs(dat1.s0(tt,:)));ylim([0 0.3]);hold on
    plot( dat2.tissue.vec.*1000-RF1.range,abs(dat2.s0(tt,:)));ylim([0 0.3]);hold on
    drawnow;
    set(gcf,'Color',[1 1 1]); % make background white
    f = getframe(gcf);
    imgif{tt} = frame2im(f);
end
