# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:40:44 2015

@author: dcarbone
"""

import pytrnsys_spf.pdata.processIceEx as load
import os
import numpy as num
import math
import string
import matplotlib.pyplot as plt
import sys

#Minutes	Hour	A001	A002	FlowDUM_Lh T_kSP_z1 T_kSP_z2 T_kSP_z3 T_kSP_z4

class ProcessIceExNum(load.ProcessIceEx):
    
    def __init__(self,_name):

        load.ProcessIceEx.__init__(self,_name)

        self.cpBrine  = 3816.0
        self.rhoBrine = 1016.0
        self.dt  = 10.
        self.nHx = 1
        self.areaHx = 2 # m2
        self.calculateMFlow = False


    def getData(self):
        

        print self.namesVariables

        for i in range(self.numberOfVariables):
            
            name = string.replace(self.namesVariables[i]," ","")
            self.namesVariables[i] = name
            
        print self.namesVariables
            
        self.tIn  = self.get("t1i")     
        
#        print self.tIn
#        print self.numberOfDataPoints
      
        self.hour = self.get("TIME")
        self.tOut = self.get("tOutHx")

        self.mFlow   = self.get("MfrCsHx1i")*self.nHx #kg/s for one Hx (parallem mass flow not included)              
        self.mFlowSI = self.mFlow/3600. #kg/s for one Hx (parallem mass flow not included)      
        
        try:
            self.qPower = self.get("sumQHx")/1000. #kW   
        except:
            self.qPower = self.get("qHxTotal")/1000. #kW  
            
        self.qPowerCalc =  self.mFlowSI*self.cpBrine*(self.tIn-self.tOut)/1000.
        
        self.qDiffPower = self.qPowerCalc-self.qPower
        
        if(sum(self.qDiffPower)>1.):
#            print Warning("Error in Q calculation qSumError (kW):",sum(self.qDiffPower))                                
#            myDiffName = self.fileName+"-diffPower.dat"
#            
#            outfile=open(myDiffName,'w')
#            lines = "!qPowerCalc qPowerRead qdiff mFlowSI\n"
#            
#            for i in range(len(self.qPowerCalc)):
#                line = "%f %f %f %f\n"%(self.qPowerCalc[i],self.qPower[i],self.qDiffPower[i],self.mFlowSI[i])
#                lines = lines+line
#            outfile.writelines(lines)
#            outfile.close()    
            
#            raise ValueError("Error in Q calculation qSumError (kW):",sum(self.qDiffPower))                    
            print Warning("Error in Q calculation qSumError (kW):",sum(self.qDiffPower))                    
        
        
        self.sumQIce    = self.get("sumQIce")/1000.
        self.sumQFused  = self.get("sumQFused")/1000.
        self.sumQLoss   = self.get("sumQLoss")/1000.
        self.sumQAcum   = self.get("sumQAcum")/1000.
        self.sumQhx     = self.get("sumQhx")/1000.
                
        imb =  sum(self.sumQhx)-sum(self.sumQAcum)- sum(self.sumQLoss) - sum(self.sumQFused) + sum(self.sumQIce)
        
        print "BALANCE imb%f"%imb
        
        self.tAmb = self.get("tAmb")
        self.UA = self.get("UA1")*self.nHx/1000. #kW/K
        self.U  = self.UA*1000./(self.nHx*self.areaHx) #W/Km2
        

        self.tStoreAv = self.get("TAvgPCM")
        
        self.ratioIce = self.get("VIceRatio")
        self.MassIce = self.get("MassOfIce")
        self.radioIce = self.get("dOutIce")/2.0
        
        self.ratioIceCalc = self.ratioIce
        
        self.tStore1 = self.get("TsensorPcm1")
        self.tStore2 = self.get("TsensorPcm2")
        self.tStore3 = self.get("TsensorPcm3")
        self.tStore4 = self.get("TsensorPcm4")
        self.tStore5 = self.get("TsensorPcm5")

        self.fConstrained = self.get("fConstrained")
        
        self.numberOfData = len(self.tIn)
        
        self.UaLog = num.zeros(self.numberOfData)
        self.ULog = num.zeros(self.numberOfData)
        self.eff = num.zeros(self.numberOfData)
        self.lmtd = num.zeros(self.numberOfData)
        
        self.HeatFlow1 = num.zeros(self.numberOfData)
        self.HeatFlow2 = num.zeros(self.numberOfData)
        self.HeatFlow3 = num.zeros(self.numberOfData)
        
        self.Tsurface1 = num.zeros(self.numberOfData)
        self.Tsurface2 = num.zeros(self.numberOfData)
        self.Tsurface3 = num.zeros(self.numberOfData)
            
        self.TstoreTS1 = num.zeros(self.numberOfData)
        self.TstoreTS2 = num.zeros(self.numberOfData)
        self.TstoreTS3 = num.zeros(self.numberOfData)
            
        self.ThxTS1 = num.zeros(self.numberOfData)
        self.ThxTS2 = num.zeros(self.numberOfData)
        self.ThxTS3 = num.zeros(self.numberOfData)
        
    def calculateEnergy(self,correctMFlow=False):


        self.eAcumQice   =  num.cumsum(self.sumQIce)*(self.dt/3600.) #kWh
        self.eAcumQFused =  num.cumsum(self.sumQFused)*(self.dt/3600.) #kWh
        self.eAcumQLoss  =  num.cumsum(self.sumQLoss)*(self.dt/3600.) #kWh
        self.eAcumQAcum  =  num.cumsum(self.sumQAcum)*(self.dt/3600.) #kWh
        self.eAcumQhx    =  num.cumsum(self.sumQhx)*(self.dt/3600.) #kWh   
        
        self.eHxPerCyclePos = num.zeros(self.numberOfData)
        self.eHxPerCycleNeg = num.zeros(self.numberOfData)
        self.eIcing = num.zeros(self.numberOfData)      
        self.qHeatPcm = num.zeros(self.numberOfData)      
        self.qCoolPcm = num.zeros(self.numberOfData)      
        
        self.eAcum = self.eAcumQhx    
                   
        time = self.hour[0]                
        
        for i in range(self.numberOfData):
        
            time = time + self.dt/3600.
            self.hour[i] = time
        
                        
            if(self.qPower[i]>=0.):
                self.eHxPerCyclePos[i] = self.eHxPerCyclePos[i-1] + self.qPower[i]*(self.dt/3600.)
                self.eHxPerCycleNeg[i] = 0.0
                self.qHeatPcm[i] = self.qPower[i]
                self.qCoolPcm[i] = 0.0
                      
            else:
                self.eHxPerCycleNeg[i] = self.eHxPerCycleNeg[i-1] + self.qPower[i]*(self.dt/3600.)
                self.eHxPerCyclePos[i] = 0.0
                self.qCoolPcm[i] = self.qPower[i]
                self.qHeatPcm[i] = 0.0
                   
            if(self.ratioIce[i]>0.):
                self.eIcing[i] = self.eIcing[i-1] + self.qPower[i]*(self.dt/3600.)
            else:
                self.eIcing[i] = self.eIcing[i-1]
#            if(self.hour[i]>=20. and self.hour[i]<=30):
#                print "q:%f mFlow:%f to=:%f ti=%f"%(self.qPower[i],self.mFlowSI[i],self.tOut[i],self.tIn[i])
            
        #CP FIX AS IN TRNSYS, RHO NO BECASUE WE USEIT TO CALC Kg/h as INPUT FOR TRNSYS
        
#        self.qPower = -self.mFlowSI*self.cpBrine*(self.tOut-self.tIn)/1000. #kW
         
        for i in range(self.numberOfData):
            #        I limited to 0 becasue if the sensor is frozen then we get a wrong water T

            self.eff[i]  = min(1,abs(self.tIn[i]-self.tOut[i])/abs(self.tIn[i]-self.tStoreAv[i]))
            self.lmtd[i] = abs(self.getTLog(self.tIn[i],self.tOut[i],self.tStoreAv[i]))      
            
            self.UaLog[i] = abs(self.qPower[i])/self.lmtd[i] #kW/K     
            self.ULog[i] = self.UaLog[i]*1000./(self.nHx*self.areaHx) #W/m2K
            
#        self.UaLog = self.UA    
#        self.ULog  = self.U
            
        
    def plotYearlyPcm(self,path):
            
        N = 2
        width = 0.65        # the width of the bars

        ind = num.arange(N)  # the x locations for the groups

        fig = plt.figure(1,figsize=(12,8))
        
        plot = fig.add_subplot(111)
        
        qHeatYear = num.zeros(N)
        qGainGroundYear = num.zeros(N)
        qReleaseYear = num.zeros(N)              
        qIceFormYear = num.zeros(N)

        qCoolYear = num.zeros(N)
        qLossGroundYear = num.zeros(N)
        qAcumYear = num.zeros(N)
        qMeltHxYear = num.zeros(N)              
        
        imbPos = num.zeros(N)
        imbNeg = num.zeros(N)
                                
#        imb = sum(self.sumQhx)-sum(self.sumQAcum)- sum(self.sumQLoss) - sum(self.sumQFused) + sum(self.sumQIce)
        
        qReleaseAcumNeg = max(sum(self.sumQAcum),0.0)/1000.
        qReleaseAcumPos = -min(sum(self.sumQAcum),0.0)/1000.
       
        qHeatYear[0]    = max(sum(self.sumQhx),0.0)/1000.
        
        qLossTotal = max(sum(self.sumQLoss),0.0)/1000.
        qGainTotal = -min(sum(self.sumQLoss),0.0)/1000.
        
        qGainGroundYear[0]   = qGainTotal
        qReleaseYear[0]      = qReleaseAcumPos        
        qIceFormYear[0]      = sum(self.sumQIce)/1000.
        
        qCoolYear[1]         = -min(sum(self.sumQhx),0.0)/1000.
        qLossGroundYear[1]   = qLossTotal
        qAcumYear[1]         = qReleaseAcumNeg
        qMeltHxYear[1]       = sum(self.sumQFused)/1000.   
        
        if(1):
            print " + qHeat:%f qGainWall:%f Qrelease:%f qIceForm:%f" % (qHeatYear[0],qGainGroundYear[0],qReleaseYear[0],qIceFormYear[0])
            print " - qCool:%f qLossWall:%f Qacum:%f    qMeltHx:%f" % (qCoolYear[1],qLossGroundYear[1],qAcumYear[1],qMeltHxYear[1])
        
        positive = qHeatYear[0] + qGainGroundYear[0] + qReleaseYear[0] + qIceFormYear[0]
        negative = qCoolYear[1] + qLossGroundYear[1] + qAcumYear[1] + qMeltHxYear[1]
        
        myImb = positive - negative         
        
        if(myImb>0.):
            imbNeg[1] = positive - negative 
        else:
            imbPos[0] = abs(positive - negative)            
                        
        barQHeat = plot.bar(ind-0.5*width, qHeatYear, width, color='r')
        barQGain = plot.bar(ind-0.5*width, qGainGroundYear, width, color='g',bottom=qHeatYear)
        barQRelease = plot.bar(ind-0.5*width, qReleaseYear, width, color='c',bottom=(qHeatYear+qGainGroundYear))
        barQIceForm = plot.bar(ind-0.5*width, qIceFormYear, width, color='0.85',bottom=(qHeatYear+qGainGroundYear+qReleaseYear))
        
        barImbPos = plot.bar(ind-0.5*width, imbPos, width, color='k',bottom=(qHeatYear+qGainGroundYear+qReleaseYear+qIceFormYear))
        
        barQCool   = plot.bar(ind-0.5*width, qCoolYear, width, color='b')
        barQLoss   = plot.bar(ind-0.5*width, qLossGroundYear, width, color='g',bottom=qCoolYear)
        barQAcum   = plot.bar(ind-0.5*width, qAcumYear, width, color='c',bottom=(qCoolYear+qLossGroundYear))        
        barQMeltHx = plot.bar(ind-0.5*width, qMeltHxYear, width, color='#CC00CC',bottom=(qCoolYear+qLossGroundYear+qAcumYear))        
        barImbNeg  = plot.bar(ind-0.5*width, imbNeg, width, color='k',bottom=(qCoolYear+qLossGroundYear+qAcumYear+qMeltHxYear))
                
        # add some            

        plot.set_ylabel('$Q$ $[MWh]$',size=15)
#        plot.set_xlabel('',size=10)
        
        box = plot.get_position()        
        plot.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        
        plot.set_title('Ice storage yearly energy flows')
        plot.set_xticks(ind)
        plot.set_xticklabels(('In', 'Out'))                       

        plot.legend( (barQHeat[0], barQGain[0],barQRelease[0], barQIceForm[0], barQCool[1], barQMeltHx[1],barImbPos[1]),\
        ('$Q_{+}$', '$Q_{gain(+)/loss(-)}$','$Q_{release(+)/acum(-)}$','$Q_{ice,form}$', '$Q_{-}$','$Q_{melt,hx}$','$Q_{melt,floating}$','$Q_{imb}$'),\
        bbox_to_anchor=(1.05,1),loc=2, borderaxespad=0.,fontsize=15)
        
        self.namePcmYearlyPlotPdf = 'IceTESYear.pdf'
        self.namePcmYearlyPlotPdfWithPath = '%s\%s' % (path,self.namePcmYearlyPlotPdf)

        plt.xlim([-0.5,1.5])
        plt.savefig(self.namePcmYearlyPlotPdfWithPath)
        plt.close()
         
if __name__ == '__main__':

    iceModel = "type861"

    name = "PCMOut.plt" 
    cmds = []

    path = ".\capillary_mats"
    case = "ClinaG-16Hx_T5_type861"

    path = r".\flat_plate"
    case = "ESSA-8Hx_T5_type861"

    path = r".\coils"
    case = "Coil_T5_type861"

    myFile = "%s\%s" % (path,name)

    nameLatex = myFile.split('.')[0]

    pilot= ProcessIceExNum(myFile)

    pilot.loadFile(skypChar="\n",verbose=True,splitArgument="\t",skypedLines=0)

    nameOut = "numPower"
    nameOutDownSized = "numPower-downSized"

    pilot.getData()
    pilot.calculateEnergy()
    pilot.printOutPower(path,name,nameOut)
    pilot.printDataEvery = 50
    pilot.printOutPower(path,name,nameOutDownSized)
    pilot.plotYearlyPcm(path)



        