function y=DNLT(ftd,c,dt)

% ftd: values of the waveform for 0;
% dt, 2dt, ... , (N-1)dt;
% c: damping factor;
% dt: time step;

order = size(ftd);
N = order(2);
n = [0:N-1];
exp1 = exp(-c*dt-1i*pi/N).^n;
fw = ftd.*exp1;
y = dt*fft(fw);

