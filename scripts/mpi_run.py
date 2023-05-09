import SepVector
import Hypercube
import genericIO
import WEM
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# par = ["","par=geom.par"]
# io=genericIO.pyGenericIO.ioModes(par)
# parObj=io.getDefaultIO().getParamObj()
# slow = genericIO.defaultIO.getVector("Vel/gauss1.H")
wave = genericIO.defaultIO.getVector("Wav/wave.H")

#     # modeling
# dshot = parObj.getFloat("dsx")
# oshot = parObj.getFloat("osx")
# nshot = parObj.getInt("ns")
# orec = parObj.getFloat("orx")
# drec = parObj.getFloat("drx")
# nrec = parObj.getInt("nr")
# nt = wave.getHyper().getAxis(1).n
# dt = wave.getHyper().getAxis(1).d
# ot = 0.
# dataLocal = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,nrec,nshot],
#     													os=[ot,orec,oshot],
#     													ds=[dt,drec,dshot]))

# create a WEM propagator
# print ("I am process %s"%rank)
# propWEM = WEM.WEM(slow,dataLocal,parObj,wave)
# propWEM.forward(False,slow,dataLocal)

# if rank == 0:
#     dshot = parObj.getFloat("dsx")
#     oshot = parObj.getFloat("osx")
#     nshot = parObj.getInt("ns")
#     orec = parObj.getFloat("orx")
#     drec = parObj.getFloat("drx")
#     nrec = parObj.getInt("nr")
#     nt = wave.getHyper().getAxis(1).n
#     dt = wave.getHyper().getAxis(1).d
#     ot = 0.
    # dataFull = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,nrec,nshot],
    #     													os=[ot,orec,oshot],
    #     													ds=[dt,drec,dshot]))
    # comm.Gather(dataLocal.getNdArray(),dataFull.getNdArray(),root=0)
    # genericIO.defaultIO.writeVector(D+"/dataGauss1.H",dataFull)

#
# # dataWin = genericIO.regFile.writeWindow(data,n1=int(nt/2))
