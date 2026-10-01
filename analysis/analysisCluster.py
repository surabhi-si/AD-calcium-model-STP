#!/home/surabhi/analysisscripts/anaconda3/bin/python

#PBS -j oe
#PBS -t 1-11
#PBS -o localhost:$PBS_O_WORKDIR/
#PBS -q batch
#PBS -l walltime=72:00:00

import os, sys

key = 'PBS_ARRAYID'
fid = os.getenv(key)
print(fid)

from time import time
sys.path.insert(0, '/home/surabhi/analysisscripts/')
from modFunc import *

i = int(fid) - 1

d =  os.environ['d']
dataDirName = d # d comes from command line input
print(dataDirName)

dataType = 'train/control/'
path = "/home/surabhi/Output/model_testing/" + dataType + dataDirName

dataPath = path
resultPath = "/home/surabhi/results/" + dataType + dataDirName

if not os.path.exists(resultPath):
    os.makedirs(resultPath)

print("\n\ndataPath : ", dataPath)
print("resultPath : ", resultPath)

ns = analysis(dataPath, resultPath)

def azAvg():
    ns.avg_dat(inFile="/dat/az.dat", outFile="/az.dat")

def caAvg():

    ns.avg_dat(
        inFile="/dat/ca.dat",
        outFile="/ca.dat"
    )

    ns.conc_calc(
        inFile="/ca.dat",
        outFile="/CaConc"
    )

    ns.pre_ca_conc(
        inFile="/ca.dat",
        outFile="/pre_ca_conc.dat"
    )

    ns.er_ca_conc(
        inFile="/ca.dat",
        outFile="/er_ca_conc.dat"
    )

def rrpAvg():
    ns.avg_dat(inFile="/dat/rrp.dat", outFile="/rrp.dat")

def vdccFlux():
    ns.avg_dat(inFile="/dat/vdcc_pq_ca_flux.dat", outFile="/vdccCaFlux.dat")
    ns.fluxCurrent(inFile="/vdccCaFlux.dat", outFile="/vdccCaFluxRate.dat")

def pmcaAvg():
    ns.avg_dat(inFile="/dat/pmca&leak_ca_flux.dat", outFile="/pmca_leak.dat")

def calbAvg():
    ns.avg_dat(inFile="/dat/calbindin_mol.dat", outFile="/calB.dat")

def sercaFluxAvg():
    ns.avg_dat(inFile="/dat/serca_ca_flux.dat", outFile="/sercaCaFlux.dat")
    ns.fluxCurrent(inFile="/sercaCaFlux.dat", outFile="/sercaCaFluxRate.dat")

def sercaMolAvg():
    ns.avg_dat(inFile="/dat/serca_mol.dat", outFile="/sercaMol.dat")

def ryrFluxAvg():
    ns.avg_dat(inFile="/dat/ryr_ca_flux.dat", outFile="/ryrCaFlux.dat")
    ns.fluxCurrent(inFile="/ryrCaFlux.dat", outFile="/ryrCaFluxRate.dat", line=2)

def ryrMolAvg():
    ns.avg_dat(inFile="/dat/ryr_mol.dat", outFile="/ryrMol.dat")

def ppf():
    isi = int(dataDirName.split("I")[1].split("V")[0])
    vdcc = int(dataDirName.split("V")[1].split('C')[0])
    print('isi: ', isi, 'vdcc: ', vdcc)
    ns.relppf(isi, vdcc)
    #ns.relppf()
    
def rel():
    ns.combineReleaseFiles()

def ptp():
    n = 20
    isi = 1000/int(dataDirName.split('p')[1].split('hz')[0])
    ns.relptp(n, isi)

# RS sims
fl = [ptp, azAvg, caAvg, vdccFlux, pmcaAvg, calbAvg, sercaFluxAvg, sercaMolAvg, ryrFluxAvg, ryrMolAvg, rrpAvg] #use ppf instead of ptp for paired-pulse analysis

fl[i]()
