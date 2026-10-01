from numpy import *
from scipy.integrate import *
from random import randint
#from pylab import *
import os, sys
from time import time
import glob
#from scipy.integrate import simpson

class analysis:

    #   Set dataPath and resultPath
    def __init__(self, dp, rp):
        self.dataPath = dp
        self.resultPath = rp

    #   Get list of directories in the path starting with 's_'
    def getDirs(self, path, sstr='s_'):
        dirs = [d for d in os.listdir(path) if os.path.isdir(path + '/' + d) and sstr in d]
        dirs.sort()
        return(dirs)

    #   Get data from dataFile in array
    def getData(self, dataFile):
       	data = []
       	f = open(dataFile, "r")
       	for l in f:
            if not l.startswith('#'): data.append([float(x) for x in l.strip('\n').strip(" ").split(" ")])
        return (data)

    #   Average Over all Seeds
    def avg_dat(self, inFile, outFile):
       	print("\nCalculating Average of", inFile)
        # Get list of directories in data path
       	dirs = self.getDirs(self.dataPath)
       	self.seeds = len(dirs)
       	print(self.seeds)

       	j=1
       	avg = genfromtxt(self.dataPath+'/'+dirs[0]+inFile)
        l = len(avg)
        for i in range(1,self.seeds):
            temp = genfromtxt(self.dataPath+'/'+dirs[i]+inFile, invalid_raise=False)
            if len(temp) != l:
                l = len(temp)
                print(i, dirs[i], l, temp[-1])
            else:
                avg += temp #genfromtxt(self.dataPath+'/'+dirs[i]+inFile)
                j += 1
        avg = avg/j #self.seeds
        print('j = ',j)

        of = self.resultPath + outFile
        print("Writing average to: " + of)
        savetxt(of, avg, fmt='%.6f')
        
    # --------------------------------------------------------
    # Compute calcium concentration from arrays
    # --------------------------------------------------------
    def computeCaConc(self, time, ca, step=5):

        c_tc = multiply(time, ca)

        dt = step * (time[1] - time[0])

        t_out = []
        c_out = []

        for i in range(0, len(time)-step-1, step):

            t_out.append(time[i])
            c_out.append((c_tc[i+step] - c_tc[i]) / dt)

        return array(t_out), array(c_out)

    #   Ca Concentration Calculation
    def conc_calc(self, step=5, inFile="/ca.dat", outFile="/CaConc"):

        time, ca = genfromtxt(
            self.resultPath + inFile,
            usecols=(0,1),
            unpack=True
        )

        print("Calculating Calcium Concentration...")

        t, conc = self.computeCaConc(time, ca, step)

        c_out = column_stack((t, conc))

        print("Writing Ca Conc. to file:" + outFile)

        savetxt(
            self.resultPath + outFile,
            c_out,
            fmt='%.6f'
        )

    # --------------------------------------------------------
    # Presynaptic (cytosolic) Ca2+ concentration from raw counts
    # --------------------------------------------------------
    #   Ca.Pre (ca.dat col 2) is an instantaneous molecule count, converted
    #   directly (no finite-difference, unlike Ca.Conc.az/conc_calc above).
    #   presynaptic_vol = 0.961 um^3 = 9.6e-16 L
    #   presynaptic_vol = 4*0.5*0.5 - er_vol in misc_*.mdl. Output is nM.
    def pre_ca_conc(self, inFile="/ca.dat", outFile="/pre_ca_conc.dat"):

        time, ca_pre = genfromtxt(
            self.resultPath + inFile,
            usecols=(0, 2),
            unpack=True
        )

        print("Calculating presynaptic Ca concentration...")

        pre_conc = ((ca_pre / 6.022e23) / 9.6e-16) * 1e9  # nM

        c_out = column_stack((time, pre_conc))

        print("Writing presynaptic Ca conc. to file:" + outFile)

        savetxt(
            self.resultPath + outFile,
            c_out,
            fmt='%.6f'
        )

    # --------------------------------------------------------
    # ER lumenal Ca2+ concentration from raw counts
    # --------------------------------------------------------
    #   Ca.ER (ca.dat col 3), same direct-conversion treatment as pre_ca_conc.
    #   er_vol = 3.9*0.1*0.1 um^3 = 3.9e-17 L, matches er_vol in misc_*.mdl.
    #   Output is uM 
    def er_ca_conc(self, inFile="/ca.dat", outFile="/er_ca_conc.dat"):

        time, ca_er = genfromtxt(
            self.resultPath + inFile,
            usecols=(0, 3),
            unpack=True
        )

        print("Calculating ER Ca concentration...")

        er_conc = ((ca_er / 6.022e23) / 3.9e-17) * 1e6  # uM

        c_out = column_stack((time, er_conc))

        print("Writing ER Ca conc. to file:" + outFile)

        savetxt(
            self.resultPath + outFile,
            c_out,
            fmt='%.6f'
        )

    #   Get Vesicle Release Statistics for PPF
    def relppf(self, isi, vdcc, resample=1000, tc=0.02, outFile="/result"): # isi in ms
       	n = 2
       	self.ts = [(i*isi+2.0)/1000.0 for i in range(n)]
       	alldirs = self.getDirs(self.dataPath)
       	ndirs = len(alldirs)
       	print('seeds: ', ndirs)

       	for d in [(self.dataPath + '/' + dir + '/dat/') for dir in alldirs]:
            os.system("cd " + d + "; cat vdcc.* > rel.dat")

        prs = []
        for r in range(resample):
            if (r+1)%100==0: print('resampling:', r+1)
            x=[randint(0,ndirs-1) for p in range(0,ndirs)]
            dirs = [alldirs[i] for i in x]

            nRel = [0]*n # [Rel1, Rel2,... Reln]
            pr = []
            cp = [0]*4 # [P00, P01, P10, P11]
            for d in [(self.dataPath + '/' + dir + '/dat/') for dir in dirs]:
                fpath = (d + "/rel.dat")
                f = open(fpath, 'r')

                #   Get no. of vesicles released after AP specified by self.ts
                temp = [0]*n
                p = [0, 0]
                time = 0
                for line in f:
                    time = float(line.strip("\n").split(" ")[0])
                    for i in range(len(self.ts)):
                        if (time>self.ts[i] and time<self.ts[i]+tc):
                            temp[i] = 1

                    '''
                    for i in range(len(self.ts)):
                        if (time>self.ts[i] and time<self.ts[i]+tc):
                            p[i] = 1
                    '''

                for i in range(n):
                    if temp[i] == 1: nRel[i] += 1

                '''
                # Calculate Conditional Release Probabilities
                if(p[0]==0 and p[1]==0): cp[0] += 1
                if(p[0]==0 and p[1]==1): cp[1] += 1
                if(p[0]==1 and p[1]==0): cp[2] += 1
                if(p[0]==1 and p[1]==1): cp[3] += 1
                '''
            '''
            print cp
            cp = [float(i)/ndirs for i in cp]
            print cp
            pp = [cp[0]/(cp[0]+cp[1]), cp[1]/(cp[0]+cp[1]), cp[2]/(cp[2]+cp[3]), cp[3]/(cp[2]+cp[3])]
            pp1 = concatenate((cp, pp), axis=0)
            '''
            for i in range(n):
                pr.append(nRel[i]/float(len(dirs)))
                print(pr)
            for i in range(1,n):
                pr.append(pr[i]/pr[0])

            '''
            for i in range(8):
                pr.append(pp1[i])
            '''

            prs.append(pr)

        m = mean(prs, axis=0)
        s = std(prs, axis=0)
        print(m,s)

        result = vstack([m,s])
        print(result)
        savetxt(self.resultPath + outFile, vstack([m,s]), fmt='%.6f')

        os.system("cat " + self.dataPath + "/*/dat/rel.dat > " + self.resultPath + "/vesRel")

        os.system("cat " + self.dataPath + "/*/dat/vdcc.async_*.dat > " + self.resultPath + "/asyncRel")
        os.system("cat " + self.dataPath + "/*/dat/vdcc.sync_*.dat > " + self.resultPath + "/syncRel")

    #   Get Vesicle Release Statistics for PTP
    def relptp(self, n, isi, resample=1000, tc=0.02, outFile="/result"): # isi in ms

        self.ts = [(i*isi)/1000.0 for i in range(n)]
        alldirs = self.getDirs(self.dataPath)
        ndirs = len(alldirs)
        print('seeds: ', ndirs)

        for d in [(self.dataPath + '/' + dir + '/dat/') for dir in alldirs]:
            if not os.path.exists(d+'/rel.dat'):
                os.system("cd " + d + "; cat vdcc.* > rel.dat")

        prs = []
        for r in range(resample):
            if (r+1)%100==0: print('resampling:', r+1)
            x=[randint(0,ndirs-1) for p in range(0,ndirs)]
            dirs = [alldirs[i] for i in x]

            nRel = [0]*len(self.ts) # [Rel1, Rel2,... Reln]
            pr = []
            ptp = []
            for d in [(self.dataPath + '/' + dir + '/dat/') for dir in dirs]:
                fpath = (d + "/rel.dat")
                f = open(fpath,'r')

                #   Get no. of vesicles released after AP specified by self.ts
                temp = [0]*n
                for line in f:
                    time = float(line.strip("\n").split(" ")[0])
                    for i in range(len(self.ts)):
                        if (time>self.ts[i] and time<self.ts[i]+tc):
                            #nRel[i] += 1
                            temp[i] = 1
                for i in range(n):
                    if temp[i] == 1: nRel[i] += 1

            for i in range(n):
                pr.append(nRel[i]/float(len(dirs)))
            for i in range(n):
                ptp.append(pr[i]/pr[0])
                pr.append(pr[i]/pr[0])

            prs.append(pr)

        m = mean(prs, axis=0)
        s = std(prs, axis=0)

        os.system("cat " + self.dataPath + "/*/dat/rel.dat > " + self.resultPath + "/vesRel")
        os.system("cat " + self.dataPath + "/*/dat/vdcc.async_*.dat > " + self.resultPath + "/asyncRel")
        os.system("cat " + self.dataPath + "/*/dat/vdcc.sync_*.dat > " + self.resultPath + "/syncRel")

        result = vstack([m,s])
        print(result)

        savetxt(self.resultPath + outFile, vstack([m,s]), fmt='%.6f')
        savetxt(self.resultPath + '/ptp', ptp, fmt='%.6f')
        
    
    def combineReleaseFiles(self):

        """
        Combine all vdcc.sync_*.dat and vdcc.async_*.dat files from every seed
        into syncRel and asyncRel in the result directory.
        """

        print("\nCombining release files...")

        # -------- Synchronus release --------
        sync_files = sorted(
            glob.glob(os.path.join(self.dataPath,
                                   "*/dat/vdcc.sync_*.dat"))
        )

        sync_out = os.path.join(self.resultPath, "syncRel")

        with open(sync_out, "w") as outfile:
            for fname in sync_files:
                with open(fname, "r") as infile:
                    outfile.write(infile.read())

        print(f"Wrote {len(sync_files)} files to {sync_out}")

        # -------- Asynchronous release --------
        async_files = sorted(
            glob.glob(os.path.join(self.dataPath,
                                   "*/dat/vdcc.async_*.dat"))
        )

        async_out = os.path.join(self.resultPath, "asyncRel")

        with open(async_out, "w") as outfile:
            for fname in async_files:
                with open(fname, "r") as infile:
                    outfile.write(infile.read())

        print(f"Wrote {len(async_files)} files to {async_out}")

        # -------- Vesicle release --------
        rel_files = sorted(
            glob.glob(os.path.join(self.dataPath,
                                   "*/dat/rel.dat"))
        )

        rel_out = os.path.join(self.resultPath, "vesRel")

        with open(rel_out, "w") as outfile:
            for fname in rel_files:
                with open(fname, "r") as infile:
                    outfile.write(infile.read())

        print(f"Wrote {len(rel_files)} files to {rel_out}")

     #   Calculate Current (pA)
    def fluxCurrent(self, inFile, outFile, step=10, ncharge=1, line=1):
        dataFile = self.resultPath + inFile
        data = self.getData(dataFile)

        charge = 1#1.602e-7 #pico Coulomb
        c_out = []
        dt=step*(data[1][0]-data[0][0])
        for i in range(0,len(data)-step-1,step):
            temp = [data[i+step][0]-dt/2]
            for l in range(1,line+1):
                temp.append(ncharge*charge*(data[i+step][l]-data[i][l])/dt)
            c_out.append(temp)

        of = self.resultPath + outFile
        print("\nWriting average to file: " + of)
        outfile = open(of,'w')
        for l in c_out:
            s = ""
            for d in l: s += str(d) + " "
            outfile.write(s + '\n')
        outfile.close()
