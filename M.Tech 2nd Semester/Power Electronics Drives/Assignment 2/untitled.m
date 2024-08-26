% Numerator
num = [(1/0.065)];
% Denominator
den = [(8*0.0047^3)  (8*0.0047^2)  (8*0.0047)  1];
% Transfer Function
G = tf(num, den)
% Plot Frequency Response
bode(G), grid