import os, sys

fileName = sys.argv[1]
outputFileName = sys.argv[2]
count = 0
with open(fileName) as textFile:
    for line in textFile:
        #os.system('tar -xzvf '+line)
        count = count + 1
        myfile = open('RunRoot.sh', 'w')
        for shLines in open('RunRootMASTER.sh'):
            if 'XXXX' in shLines:
                shLines = shLines.replace('XXXX',line.rstrip())
            if 'YYYY' in shLines:
                shLines = shLines.replace('YYYY',outputFileName+'_'+str(count)+'.root')
            myfile.write(shLines)
        myfile.close()
        os.system('bash RunRoot.sh')
        print "Files done: ", count
    
    os.system('hadd '+outputFileName+'_All.root '+outputFileName+'_*.root')
    
