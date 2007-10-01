clear all
%A = -1.0/log(0.25);
A = 0.5;
fplot(@(x)[exp(-x/(A*100)), exp(-x/(A*1000))],[0 5000]);