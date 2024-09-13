
import numpy as np
import matplotlib.pyplot as plt 
import scipy.signal as signal
import Filtres as filt
import math

def signal_resolution(ECG, PPG):
    resolution=0.0032
    ECG_gain=972
    for i in range(len(ECG)):
        if ECG[i]>512:
            ECG[i]=(ECG[i]-512)*resolution/ECG_gain*math.pow(10,3)
        elif ECG[i]==512:
            ECG[i]=0
        elif ECG[i]<512:
            ECG[i]=-1*(512-ECG[i])*resolution/ECG_gain*math.pow(10,3)
        PPG[i]=(PPG[i]-512)*resolution
    return (ECG, PPG)


def max_derivative(PPG,det_th_der,atime):
    i=0
    long=len(PPG)
    PPG_der=np.zeros(long)
    
    for i in range(long):
        if (i!=0 and i!=long-1 and i!=1 and i!=long-2):
            PPG_der[i]=(2*(PPG[i]-PPG[i-1])+PPG[i+2]-PPG[i-2])/8*10
        if i==2:
            PPG_der[0]=PPG_der[i]
            PPG_der[1]=PPG_der[i]
        if i==long-3:
            PPG_der[i+1]=PPG_der[i]
            PPG_der[i+2]=PPG_der[i]
    
    time_max_der=[]
    points_max_der=[]
    time_plot_max_der=[]
    x=0
    y=0
    th_der=0.05
    for m in range(len(PPG_der)):
        if (m!=0 and m!=len(PPG_der)-1):
            x=PPG_der[m]-PPG_der[m-1]
            y=PPG_der[m]-PPG_der[m+1]
        if (x>0 and y>0 and PPG_der[m]>det_th_der):
            if PPG[m]>0.2:
                points_max_der.append(PPG[m])
                time_max_der.append(atime[m])
                time_plot_max_der.append(atime[m])
    return time_max_der,points_max_der,time_plot_max_der

def time_scale(PPG):
    initial_time=0
    freq=1000
    end_time = (1/freq)*len(PPG)
    time=np.linspace(initial_time, end_time,len(PPG))
    return time

def ECG_peaks(det_th,filtered_ECG,time):
    pics_ECG=[]
    time_peaks=[]
    time_pics_ECG=[]
    
    min_p=0
    state=0
    x=0
    y=0

    for m in range(len(filtered_ECG)):
        if state==0:
            if (m!=0 and m!=len(filtered_ECG)-1):
                x=filtered_ECG[m]-filtered_ECG[m-1]
                y=filtered_ECG[m]-filtered_ECG[m+1]
            if (x>0 and y>0 and filtered_ECG[m]>det_th):
                pics_ECG.append(filtered_ECG[m])
                time_pics_ECG.append(time[m])
                time_peaks.append(time[m])
                pics_ECG.append(filtered_ECG[min_p])
                time_peaks.append(time[min_p])
                state=1
            if (x<0 and y<0):
                min_p=m
        else:
            x=filtered_ECG[m]-filtered_ECG[m-1]
            y=filtered_ECG[m]-filtered_ECG[m+1]
            if (x<0 and y<=0):
                pics_ECG.append(filtered_ECG[m])
                time_peaks.append(time[m])
                state=0

    return (pics_ECG,time_peaks, time_pics_ECG)

def heartrate_ECG(time_pics_ECG):
    heart_b_t_b_Ecg=0
    heart_Ecg=[]
    for z in range(len(time_pics_ECG)):
        if z!=len(time_pics_ECG)-1 :
            heart_Ecg.append(time_pics_ECG[z+1]-time_pics_ECG[z])
            heart_b_t_b_Ecg+=(time_pics_ECG[z+1]-time_pics_ECG[z])/len(time_pics_ECG)
    
    time_between_pulses_ECG = heart_b_t_b_Ecg #time between two pulses
    heart_rate_ECG = 60/heart_b_t_b_Ecg #heart rate 

    return time_between_pulses_ECG, heart_rate_ECG

def PPG_peaks(det_th, filtered_PPG, time):
    
    pics_pulsoximetre = []
    time_pics_pulso=[]

    x=0
    y=0

    for m in range(len(filtered_PPG)):
        if (m!=0 and m!=len(filtered_PPG)-1):
            x=filtered_PPG[m]-filtered_PPG[m-1]
            y=filtered_PPG[m]-filtered_PPG[m+1]
        if (x>0 and y>0 and filtered_PPG[m]>det_th):
            pics_pulsoximetre.append(filtered_PPG[m])
            time_pics_pulso.append(time[m])

    return (time_pics_pulso, pics_pulsoximetre)

def heartrate_PPG(time_pics_pulso):
    
    heart_b_t_b_pulso=0
    heart_pulso=[]
    for z in range(len(time_pics_pulso)):
        if z!=len(time_pics_pulso)-1 :
            heart_pulso.append(time_pics_pulso[z+1]-time_pics_pulso[z])
            heart_b_t_b_pulso+=(time_pics_pulso[z+1]-time_pics_pulso[z])/len(time_pics_pulso)
    
    time_between_pulses_PPG = heart_b_t_b_pulso #time between two pulses
    heart_rate_PPG = 60/heart_b_t_b_pulso #heart rate from ppg

    return (time_between_pulses_PPG, heart_rate_PPG)

def arrival_time(time_pics_ECG, time_pics_pulso,time_max_der,g):
    
    Arrival_time=[]
    av_atime=0
    lon=0
    lecg=len(time_pics_ECG)
    lppg=len(time_pics_pulso)
    if time_pics_ECG[0]>time_max_der[0]:
        time_max_der.pop(0)
        g=1
                
    if time_pics_ECG[lecg-1]>time_max_der[lppg-1-g]:
        time_pics_ECG.pop(lecg-1)
        g=0
    
    if (lecg>=lppg):
        lon=lppg
    else:
        lon=lecg

    for h in range(lon):
        Arrival_time.append(time_pics_pulso[h]-time_pics_ECG[h])
        av_atime+=(time_pics_pulso[h]-time_pics_ECG[h])/lon

    average_arrival_time = av_atime
    
    return(Arrival_time, average_arrival_time,g)


def meanBloodPressure (calibratedMBP, gammaConstant, calibratedPTT, calculatedPTT):
    MBP = calibratedMBP + (2/gammaConstant)*math.log(calibratedPTT/calculatedPTT)
    return MBP


def bloodPressure (MBP, calibratedPP, calibratedPTT, calculatedPTT):
    DBP = MBP - (1/3)*calibratedPP*math.pow((calibratedPTT/calculatedPTT),2)
    SBP = DBP + calibratedPP*math.pow((calibratedPTT/calculatedPTT),2)
    return (DBP, SBP)
