import numpy as np
import matplotlib.pyplot as plt

from modules.tracker import *
from modules.naff import *
from modules.grid import *
from modules.tune_resonances import *
from modules.FMA import *
from modules.NAFF import *
from modules.fft_naff import *

############ Tune#####################


for i in range (2):
    if (i==0):
        data=[]
        data_init =[]
        data_f1   =[]
        data_final=[]
        #print "sin signal"
        t=np.arange(0.0,0.1,0.0001)
        for i in t:
              #data.append(8*np.sin(2*3.141592*1700.33*i)+2.0*np.sin(2*3.141592*2033.333*i))#+13*np.sin(2*3.141592*3000*1*i))
              #data.append((10*np.cos(2*3.141592*100*i+6.28/6)+13*np.cos(2*3.141592*300*i+6.28*0.2))*(1/(i+1.)**10))
              ##data.append((13*np.cos(2*3.141592*300*i+6.28*0.2*0.0)))#*(1/(i+1.)**10))
              data.append((1000*np.cos(2*3.141592*300.0*i))+ ((500*np.cos(2*3.141592*1000*i))))#*(1/(i+1.)**10))
              #data.append(8*np.sin(2*3.141592*100*i+1.*6.28/6))
              #data.append(8*np.sin(2*3.141592*1000*i+0.*6.28/6)+13*np.sin(2*3.141592*3000*i+0*6.28/5))
         #data.append(0.5*np.cos(2*3.1415*1000*i))
        coord=Vec_cpp()
        zero=Vec_cpp()
        coord.extend(i for i in data)
        zero.extend(0.0*i for i in data)
        naff=NAFF()
        tune_all=naff.get_f1(coord,zero)
        #print tune_all 
        #tune_all=get_f1(coord,zero)
        #for i in tune_all:
        #    print i
        ###print tune_all
        #plt.plot(data,linestyle='-')
        #plt.show()
        quit()
    elif (i==2):
        quit()
    #    data = []
    #    data_init =[]
    #    data_f1   =[]
    #    data_final=[]
    #    f=open('LHC_Data','r')
    #    for line in f:
    #        parts=line.split()
    #        data.append(float(parts[0]))
    #    coord=Vec_cpp()
    #    zero=Vec_cpp()
    #    coord.extend(i for i in data)
    #    zero.extend(0*i for i in data)
    #    tune_all=NAFF_f(coord,zero)
    #    for i in tune_all:
    #        print i
############################################## LHC DATA ##############################3
    f = open('text.dat','r')
    for line in f:
        parts = line.split()
        data_init.append( complex (float(parts[0]), float(parts[1])) )
        data_final.append( complex (float(parts[4]), float(parts[5])) )
        data_f1.append( complex (float(parts[2]), float(parts[3])) )
    t_final  = np.arange( len(data_final) )
    t_f1 = np.arange(len(data_f1))
    t_init = np.arange(len(data_init))
    sp_final = np.fft.fft(data_final)
    sp_f1= np.fft.fft(data_f1)
    sp_init= np.fft.fft(data_init)
    freq_final  = np.fft.fftfreq(t_final.shape[-1])
    freq_f1 = np.fft.fftfreq(t_f1.shape[-1])
    freq_init = np.fft.fftfreq(t_init.shape[-1])
    plt.plot(freq_init,(sp_init))
    plt.xlabel(r'Frequency', fontsize=20)
    plt.ylabel(r'Amplitude', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title('FFT for initial signal')
    plt.tight_layout()
    plt.grid(True)
    plt.show()
    plt.plot(freq_f1,(sp_f1))
    plt.xlabel(r'Frequency', fontsize=20)
    plt.ylabel(r'Amplitude', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title('FFT for data that will be subtracted C++')
    plt.tight_layout()
    plt.grid(True)
    plt.show()
    
    plt.plot(freq_final,(sp_final))
    plt.xlabel(r'Frequency', fontsize=20)
    plt.ylabel(r'Amplitude', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.title('FFT for final signal')
    plt.tight_layout()
    plt.grid(True)
    plt.show()
    
    plt.plot(data_init)
    #plt.show()
    plt.plot(data_f1)
    plt.plot(data_final,ms=4,linewidth=4,color='r')
    plt.show()
    quit()
    
    
    



