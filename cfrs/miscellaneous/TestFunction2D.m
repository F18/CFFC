function f = TestFunction2D(x,y)
f = (y-0.24)*(y-1.5);
f = f*(y+0.25);
f = f*(x-0.25)*(x-0.5)*(x+0.25);
