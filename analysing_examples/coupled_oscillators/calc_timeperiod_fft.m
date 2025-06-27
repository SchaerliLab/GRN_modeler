function T = calc_timeperiod_fft(t,c)
% calculate the time period with fft
% t: time it has to be uniformly distributed !!!!!!
% c: the variable

Y = fft(c);
L = length(t);
Fs = 1/(t(2)-t(1));

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

% plot(f(2:end),P1(2:end)) 
% title("Single-Sided Amplitude Spectrum of X(t)")
% xlabel("f (Hz)")
% ylabel("|P1(f)|")

pos = find(P1(2:end)==max(P1(2:end)),1,'first');
T = 1./f(1+pos);

end

