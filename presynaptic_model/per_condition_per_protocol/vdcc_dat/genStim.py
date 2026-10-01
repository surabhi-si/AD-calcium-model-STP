from pylab import *

# Note: stimuli are 5 ms wide.
isi = 100e-3
for n in [1]:  # n = number of pules in the train-stim
    d = 100e-3    # d = delay in first pulse

    infile=[0 for _ in range(9)]
    outfile=[0 for _ in range(9)]

    infile[0]="VDCC_PQ_C01.dat"
    infile[1]="VDCC_PQ_C10.dat"
    infile[2]="VDCC_PQ_C12.dat"
    infile[3]="VDCC_PQ_C21.dat"
    infile[4]="VDCC_PQ_C23.dat"
    infile[5]="VDCC_PQ_C32.dat"
    infile[6]="VDCC_PQ_C34.dat"
    infile[7]="VDCC_PQ_C43.dat"
    infile[8]="VDCC_PQ_Ca_1.2mM.dat"

    # Generate outfile
    for (i,f) in enumerate(infile):
        f = f.split(".dat")[0]
        #outfile[i] = f + f'_{int(isi*1000):d}ms_ppf.dat'
        outfile[i] = f + '_test.dat'
        print (outfile[i])

    # Generate Train-Stim
    for i in range(len(infile)):

        # Read The Stim
        idata = genfromtxt(infile[i])

                # Make The Stim
        ofile = open(outfile[i], 'w')
        duration = idata[-1][0]

        # Train-stim with number of spikes n and interspike interval isi
        for i in range(n):
            tdelay = d + i*isi
            for j in range(len(idata)):
                ofile.write(str(idata[j][0]+tdelay)+"\t"+str(idata[j][1])+"\n")
