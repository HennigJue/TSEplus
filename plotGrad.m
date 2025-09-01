function Grad = plotGrad(Grad)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
t=[0 Grad.riseTime Grad.riseTime+Grad.flatTime Grad.riseTime+Grad.flatTime+Grad.fallTime];
s=[0 Grad.amplitude Grad.amplitude 0];
%figure
plot(t,s)
end