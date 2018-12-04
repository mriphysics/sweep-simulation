function [dat, tissue, RF, motion] = sweep_sim_EPG_2(tissue, RF, motion)
%% Laurence Jackson, BME, KCL, 2018
%
% Simulates pulse profile moving across tissue with motion and flow
% consideration
%
% INPUTS::
%   tissue  - tissue struct with fields
%       T1      - T1 of tissue
%       T2      - T2 of tissue
%       length  - length of tissue to simulate
%
%   RF      - sequence struct with fields
%       profile - RF profile (flip radians vs time)
%       TR      - TR
%       thk     - nominal slice thickness
%       swp     - sweep rate
%
% OUTPUTS::
%       dat     - return data structure with fields
%       s0      - signal
%       flipmat - applied flips (time vs zloc)
%

%% Hidden options - shouldnt need to be changed in most cases
% tissue_multiplier = 5; % resolution of tissue vector
RF.range = 10; % +/-mm to simulate rf pulse over (i.e. extend of sidebands to include)

offload = 1; % offload to remote machine if set up

%% calculated values

[RF.profile, zz,~] = pulse_profile(RF); % simulate pulse profile

if RF.block ==1
    RF.profile(RF.profile>(0.5*(max(RF.profile)))) = max(RF.profile);
    RF.profile(RF.profile<(0.5*(max(RF.profile)))) = 0;
end

if RF.seqspec == 1
    RF.npulses = RF.npe * RF.ndyn *  RF.nslice;
    warning('RF.npulses is being overridden by RF.seqspec and RF.npe -- npulses is now %d',RF.npulses)
else
    RF.nslice = ceil(RF.npulses./RF.npe)+1; % estimated number of slices to make sure enough tissue is simulated
end

RF.pulseshift = 0.01.*RF.swp.*RF.thk; % pulse shift in m (RF.swp = % slice moved per pulse)

RF.flip = rad2deg(max(RF.profile)); % check approximate flip angle
RF.sweepdur = RF.npulses.*RF.TR*0.001; % sweep duration in seconds

%% calculate motion vectorss
motion_resp = sin(linspace(0,2*pi*RF.sweepdur*(motion.respfreq),RF.npulses)).*(motion.respmag./2);

%% Introduce flow component
motion.flow_per_pulse = 0;
if ~motion.flow==0   
    motion.flow_per_pulse = (motion.flow * RF.sweepdur) / RF.npulses; % displacement per pulse from flow
    motion.flow_dist = motion.flow_per_pulse .* RF.npulses;
end

%% define tissue
tissue.min = min([0,-max(motion_resp),(motion.flow_per_pulse.*RF.npulses)]);
tissue.max = max([max(motion_resp),(motion.flow_per_pulse.*RF.npulses)]);
tissue.length = (RF.range*1e-3) + tissue.min + ((RF.thk + RF.slicegap) * RF.nslice) + (RF.pulseshift.*RF.npulses) + (abs(motion.flow_per_pulse).*RF.npulses);
tissue_resolution = (tissue.length*1e3) * 100; % 100 elements per mm
tissue.vec = linspace(tissue.min,tissue.length,tissue_resolution);

%% Introduce flow component
% motion.flow_per_pulse = 0;
% if ~motion.flow==0

    % extend tissue.vec to include flowing spins
%     tissue.vec_short = tissue.vec;
%     tissue.vec_short = tissue.vec;
%     tissue.vec = linspace(0,tissue.length+motion.flow_dist,RF.npulses.*tissue_multiplier.*(motion.flow_dist/tissue.length));
%     tissue.vec_long = tissue.vec; % store this for later
% end

%% print final simulation paramters
print_sim_info(tissue, RF, motion)

%% Produce flipmat
zzabs = (zz) + abs(min(zz));
flipmat = zeros(RF.npulses,length(tissue.vec));
sliceshift = 0;
firstpulse = [1];
sliceidx = 1;
dynidx = 1;

switch RF.sliceorder
    case 'ascending'
        slicev = 0:(RF.nslice-1);
    case 'descending'
        slicev = (RF.nslice-1):-1:0;
    case 'odd-even'
        v = 0:(RF.nslice-1);
        v_odd = v(rem(v,2)~=0);
        v_even = v(rem(v,2)==0);
        slicev = [v_even v_odd];
    case 'random'
        vv = 0:(RF.nslice-1);
        slicev = vv(randperm(length(vv)));
    otherwise
        error('Check RF.sliceorder definitition')
end

for puls = 1:RF.npulses
    if RF.seqspec == 1 % npulses defined by nslices and ndyns
        switch RF.dynorder
            case 'slices'
                
                if RF.pulseshift == 0 && mod(puls-1,RF.npe) == 0 % not sweep and first pulse in 2D k-space
                    
                    sliceshift = slicev(sliceidx).*(RF.thk + RF.slicegap);
                    
                    if sliceidx < RF.nslice
                        sliceidx = sliceidx + 1;
                    else
                        dynidx = dynidx + 1;
                        sliceidx = 1;
                    end
                    
                    firstpulse = [firstpulse; puls];
                    
                end
                
            case 'dynamics'
                if RF.pulseshift == 0 && mod(puls-1,RF.npe) == 0
                    sliceshift = slicev(sliceidx).*(RF.thk + RF.slicegap);
                    
                    if dynidx < RF.ndyn
                        dynidx = dynidx + 1;
                    else
                        sliceidx = sliceidx + 1;
                        dynidx = 1;
                    end
                    
                    firstpulse = [firstpulse; puls];
                end
        end
    else % normal behaviour
        if RF.pulseshift == 0 && mod(puls,RF.npe) == 0 % not sweep and first pulse in 2D k-space
            sliceshift = sliceshift + RF.thk + RF.slicegap;
            firstpulse = [firstpulse; puls];
        end
    end
    
    xx =  zzabs + (puls - 1).*RF.pulseshift + sliceshift + motion_resp(puls) + (puls - 1).*motion.flow_per_pulse; % where the pulse IS
    
    zq = find((tissue.vec >= xx(1)) & (tissue.vec <  (zzabs(end) + xx(end)))); % index of these locations in tissue vector
    
    flipvec = interp1(xx,RF.profile,tissue.vec(zq),'linear');
    flipmat(puls,zq) = flipvec;
    dat.offset(puls) = xx(1);
    
end
figure();imagesc(tissue.vec.*1000-RF.range,1:RF.npulses,rad2deg(flipmat));

flipmat(isnan(flipmat)==1) = 0; % remove nans

% Include catalysation pulses
if ~isempty(RF.catalysation)
    for rr = 1:length(RF.catalysation)
        for ff = 1:length(firstpulse)
            fliploc = firstpulse(ff)+(rr-1);
            if fliploc > size(flipmat,1)
                continue;
            end
            v = flipmat(fliploc,:);
            v = (v - min(v(:))) / (max(v(:)) - min(v(:)));
            flipmat(fliploc,:) = v*deg2rad(RF.catalysation(rr));
        end
    end
end

%% EPG
phi = RF_phase_cycle(length(flipmat(:,1)), RF.seq); % phase cycling scheme

if offload == 1
    SS.flipmat = flipmat;
    SS.phi = phi;
    SS.RF = RF;
    SS.tissue = tissue;
    s0_RF = send2remote('EPG_sim_offload',SS);
    delete('temp_struct.mat')
else
    parfor zz = 1:size(flipmat,2)
        s0_RF(:,zz) = EPG_GRE(flipmat(:,zz),phi,RF.TR,tissue.T1,tissue.T2);
    end
end

%% convert to scanner co-ordinates- space in which signals are measured
% tissue.vec = linspace(tissue.min,tissue.length,RF.npulses.*tissue_multiplier);  % redeclare tissue.vec to remove flow extension if it exists
sliceshift = 0;
s0 = zeros([RF.npulses,length(tissue.vec)]);
for puls = 1:RF.npulses
    
%     if ~motion.flow==0
%         xx =  tissue.vec_long + sliceshift - motion_resp(puls) - (puls - 1).*motion.flow_per_pulse;
%     else
        xx =  tissue.vec + sliceshift - motion_resp(puls) - ((puls - 1).*motion.flow_per_pulse);
%     end
    
    zq = find((tissue.vec >= xx(1)) & (tissue.vec <  xx(end)));% index of these locations in tissue vector
    
    flipvec = interp1(xx,s0_RF(puls,:),tissue.vec(zq),'linear');
    s0(puls,zq) = flipvec;
    
%     if ~motion.flow==0
%         xx_pr =  tissue.vec_long - (puls - 1).*RF.pulseshift - sliceshift - motion_resp(puls) - (puls - 1).*motion.flow_per_pulse; % where the pulse IS
%         zq_pr = find((tissue.vec_long >= xx_pr(1)) & (tissue.vec_long <  xx_pr(end))); % index of these locations in tissue vector
%         qq = linspace(tissue.vec_long(zq_pr(1)),tissue.vec_long(zq_pr(end)),1000);
%     else
        xx_pr =  tissue.vec - (puls - 1).*RF.pulseshift - sliceshift - motion_resp(puls) - (puls - 1).*motion.flow_per_pulse; % where the pulse IS
        zq_pr = find((tissue.vec >= xx_pr(1)) & (tissue.vec <  xx_pr(end))); % index of these locations in tissue vector
        qq = linspace(tissue.vec(zq_pr(1)),tissue.vec(zq_pr(end)),1000);
%     end

    dat.profile(1,:,puls) = qq;
    dat.profile(2,:,puls) = interp1(xx_pr,s0_RF(puls,:),qq,'linear');
    
end

figure();imagesc(tissue.vec.*1000-RF.range,1:RF.npulses,abs(s0));

s0(isnan(s0)==1) = 0; % remove nans

%% Bring results into dat
dat.s0 = s0; % signal in scanner space
% dat.s0_RF = s0_RF; % signal in RF space - useful for debugging
dat.flipmat = flipmat;
dat.RF = RF;
dat.mnotion = motion;
dat.tissue = tissue;
