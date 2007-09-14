clear all

StartPointX = -0.00000016;
EndPointX = 1.0E-14;
StartPointY = -1.15E-14;
EndPointY = 0.24;

result = dblquad("TestFunction2D",StartPointX,EndPointX,StartPointY,EndPointY);