import os, sys, glob
#### run like: python OpenLog.py <output root file name>
outputFileName = sys.argv[1]
count = 0
os.system('bash MakeLib.sh')
os.system('ls *.tgz > file75V.txt')
fileName = 'file75V.txt'
with open(fileName) as textFile:
    for line in textFile.readlines():
        ### untar the main file
        os.system('tar -xzvf '+str(line.rstrip())+' > tarFile.txt')
        tarFile = open('tarFile.txt')
        tarAllLine = tarFile.readlines()
        tarLine = tarAllLine[1].split('/')
        print 'tarLine[0] : ', tarLine[0]
        filePath = tarLine[0].rstrip()+'/log.RAWtoALL'
        
        count = count + 1
        print 'filePath: ', filePath
        os.system('grep -i "Arka" '+filePath+' > '+outputFileName+str(count)+'.txt')
        
        myfile = open('RunRoot.sh', 'w')
        for shLines in open('RunRootMASTER.sh'):
            if 'XXXX' in shLines:
                shLines = shLines.replace('XXXX', outputFileName+str(count)+'.txt')
            if 'YYYY' in shLines:
                shLines = shLines.replace('YYYY',outputFileName+'_'+str(count)+'.root')
            myfile.write(shLines)
        myfile.close()
        os.system('bash RunRoot.sh')
        print "Files done: ", count
        if(count > 120):
            break
        
os.system('hadd '+outputFileName+'_All.root '+outputFileName+'_*.root')
os.system('mv '+outputFileName+'_All.root /eos/user/a/asantra/ForTaka/')