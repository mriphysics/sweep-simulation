function [pulse_prof,zz,mxy] = pulse_profile(RF)
% script to find slice profile for given RF shape
%


%% parse
if strcmpi(RF.seq,'bssfp')==1 % balanced pulse
    disp('Simulating RF pulse - SG_100_100_0')
    load('sg_100_100_0.mat'); % pulse shape
    TBP = 2.1362; % time-bandwidth product SG_150_100_167
    b1max =5.1834e-3; % b1max for 90deg flip in mt - use to control flip angle
    dur = 2.1056e-3;
elseif ~ischar(RF.seq)
    disp('Simulating RF pulse - SG_150_100_167')
    load('sg_150_100_167.mat'); % pulse shape
    TBP = 3.0469; % time-bandwidth product SG_150_100_167
    b1max = 13.4333e-3; % b1max for 90deg flip in mt - use to control flip angle
    dur = 1.2608e-3;
else
    error('check RF pulse definitition')
end

b1max_fact = b1max/90;
b1max =  b1max_fact*RF.flip;
 
pulse = pulse(:);
pulse = b1max*pulse/max(pulse);
M1=length(pulse);

dt = 6.4e-6;

%    interpolate M*dt to match pulse duration
M2 = floor(dur ./ dt);
pulse = interp1(1:length(pulse),pulse,linspace(1,length(pulse),M2),'nearest');
pulse = pulse(:);
M=length(pulse);

%c=consts;
c.gamma_uT = 267.5221;
gamma = c.gamma_uT*1e3;

ns=512;
z = linspace(-RF.range*1e-3,+RF.range*1e-3,ns); % update this as appropriate for pulse profile
pos = [z(:)*0 z(:)*0 z(:)];
G = zeros([M 3]);
Gsel = 2*pi*TBP/((M-1)*dt*RF.thk*gamma);
G(:,3) = Gsel;
[m,mz] = (blochsim_CK(pulse,G,pos,ones([ns 1]),zeros([ns 1]),'dt',dt));
%     plot(pos(:,3)*1e3,abs(m));



zz = pos(:,3);
mxy = abs(m);
pulse_prof = acos(mz);


end
