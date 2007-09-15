import re, sys, os
from audioop import avg
from math import log10

# Paterns for identification
Order_pat = "\\s*Reconstruction_Order\\s*=\\s*"
Funct_pat = "\\s*Function_Definition\\s*=\\s*"
Norm_pat = "\\s*Norm[\w =]*(\\d+\.?\\d+(e-+)?\\d+)"

path = os.environ["Rsrc"]
os.chdir(path)
name = os.path.join(path,"reconstruct.in")
if os.path.isfile(name):
    input = file(name,"r")
    contents = input.readlines()
    input.close()

# set the tests
testFunction = ["Default", "Example1", "Example2", "Example3", "Example4"]
testFunction = ["Example3"]
Order = [1,2,3,4,5,6]

outContents = contents[0:len(contents)]
tol = 1.0e-15
for test in testFunction:
    Plot_maximum = []
    Plot_average = []
    for order in Order:
        Residual = []
        output = file("DataOut.dat","w")
        # parse the contents
        index = 0
        for line in contents:
            if re.search(Funct_pat,line):
                outContents[index+1] = test+"\n"
            if re.search(Order_pat,line):
                outContents[index+1] = str(order)+"\n"
            index += 1
        output.writelines(outContents)
        output = os.popen("reconstruct -f DataOut.dat","r")
        result = output.readlines()
        for line in result:
          #  print line,
            mo = re.search(Norm_pat, line)
            if mo:
           #     print "String is:", mo.group(1)
           #     print "Float is:", float(mo.group(1))
                Residual.append(float(mo.group(1)))
#        print Residual
        maximum = max(Residual)
        average = 0.0
        index = 0
        for r in Residual:
            if r > tol:
                index += 1
                average += r
        #average = sum(Residual)
        if index > 0:
            average /= index
        Plot_maximum.append([log10(order),log10(maximum)])
        Plot_average.append([log10(order),log10(average)])
    # write the tecplot file
    name = "ResVersusOrder_" + test + ".dat"
    output = file(name, "w")
    print >> output, "TITLE = ", test
    print >> output, "VARIABLES = order  residual"
    print >> output, "ZONE"
    for (first,second) in Plot_maximum:
        print >> output, first, second
    print >> output, "ZONE"
    for (first,second) in Plot_average:
        print >> output, first, second
    print Plot_maximum

output.close()
#print Plot_maximum
#print Plot_average

