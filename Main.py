import numpy as np
import matplotlib.pyplot as plt
import time
import bitalino
import Filtres as filt
import ipywidgets
import functions as fun

mode=input("0 to use bitalino adquisition, 1 to use txt:")

if mode=='0':
    macAddress = "98:D3:71:FD:61:F2"
    
    device = bitalino.BITalino(macAddress)
    time.sleep(1)

    srate = 100
    nframes = 1000
    threshold = 5
    device.start(srate, [1,5])
    win_time=10
    len=win_time*srate

print ("START")
# Function to generate a sample segment
def generate_segment(length=10000):
    return np.random.rand(length)

# Function for continuous plotting
def continuous_plot(segment_length=10000, interval=10, display_length=20000):
    i=0
    x=0
    g=0
    t=0
    order = 8
    fs = 1000
    cutoff_ecg = 30
    cutoff_ppg = 10
    data=[]
    det_th=0.4
    det_th_der=0.04
    Ecg1=np.zeros(10000)
    Ecg2=np.zeros(10000)
    Ppg1=np.zeros(10000)
    Ppg2=np.zeros(10000)
    Ecg1_f=np.zeros(10000)
    Ecg2_f=np.zeros(10000)
    Ppg1_f=np.zeros(10000)
    Ppg2_f=np.zeros(10000)
    lecg=0
    lppg=0
    
    calibratedPTT=float(input("Escriu el valor del PTT de calibració:"))
    Cal_SBP=float(input("Escriu el valor del SBP de calibració:"))
    Cal_DBP=float(input("Escriu el valor del DBP de calibració:"))
    calibratedPP=Cal_SBP-Cal_DBP
    calibratedMBP=(Cal_SBP+Cal_DBP)/2
    gammaConstant = 0.017       # Constant gamma for MBP calculation

    heart_rate_ECG=0
    heart_rate_PPG=0
    average_arrival_time=0

    plt.ion()  # Turn on interactive mode for continuous updating

    try:
        signal_ecg = np.zeros(display_length)  # Initialize an array for the displayed signal
        signal_ppg = np.zeros(display_length)
        atime = np.zeros(display_length)
        while True:
            
            print(i)
            # Generate a new segment
            if mode=='0':
                data = device.read(nframes)
            else:
                with open("ECG_PPG_movement.txt", "r") as f:
                    for linia in f.readlines():
                        liniaStrip = linia.strip() #Eliminem caràcters als extrems
                        data.append(liniaStrip.split("\t"))  #Separem la cadena anterior
                data=data[t*10000:(t+1)*10000]
            if i==0:
                for j in range(10000):

                    Ecg1[j]=float(data[j][5])
                    Ppg1[j]=float(data[j][6])
                Ecg1,Ppg1=fun.signal_resolution(Ecg1,Ppg1)
                
                Ecg1_f = filt.butter_lowpass_filter(Ecg1, cutoff_ecg, fs, order)
                Ppg1_f = filt.butter_lowpass_filter(Ppg1, cutoff_ppg, fs, order)
                
                i=1
                
            else:
                for j in range(10000):
                    
                    Ecg2[j]=float(data[j][5])
                    Ppg2[j]=float(data[j][6])
                Ecg2,Ppg2=fun.signal_resolution(Ecg2,Ppg2)
                
                Ecg2_f = filt.butter_lowpass_filter(Ecg2, cutoff_ecg, fs, order)
                Ppg2_f = filt.butter_lowpass_filter(Ppg2, cutoff_ppg, fs, order)
                
                i=0
                
            
            # Update the signal array with the new segment
            if x==0:
                signal_ecg=np.concatenate((Ecg2_f,Ecg1_f))
                signal_ppg=np.concatenate((Ppg2_f,Ppg1_f))
            
                x=1
            else:  
                signal_ecg=np.concatenate((Ecg1_f,Ecg2_f))
                signal_ppg=np.concatenate((Ppg1_f,Ppg2_f))
                
                x=0
            
            print(max(signal_ecg))
            atime=np.linspace((t-1)*10,(t+1)*10,20000)
            pics_ECG,time_peaks, time_pics_ECG=fun.ECG_peaks(det_th,signal_ecg,atime)
            time_between_pulses_ECG, heart_rate_ECG=fun.heartrate_ECG(time_pics_ECG)
            time_pics_pulso, pics_pulsoximetre=fun.PPG_peaks(det_th,signal_ppg,atime)
            time_between_pulses_PPG, heart_rate_PPG=fun.heartrate_PPG(time_pics_pulso)
            time_max_der,points_max_der,time_plot_max_der=fun.max_derivative(signal_ppg,det_th_der,atime)
            print(type(time_pics_ECG))
            print(type(time_pics_ECG[0]))
            # lecg=len(time_pics_ECG)
            # lppg=len(time_max_der)
            # if time_pics_ECG[0]>time_max_der[0]:
            #     time_max_der.pop(0)
            #     h=1
                
            # if time_pics_ECG[lecg-1]>time_max_der[lppg-1-h]:
            #     time_pics_ECG.pop(lecg-1)
            #     h=0
            print(time_pics_ECG)
            print(time_max_der)
            
            Arrival_time, average_arrival_time,g=fun.arrival_time(time_pics_ECG, time_max_der,time_max_der,g)
            
            MBP = fun.meanBloodPressure (calibratedMBP, gammaConstant, calibratedPTT, average_arrival_time)
            (DBP, SBP) = fun.bloodPressure (MBP, calibratedPP, calibratedPTT, average_arrival_time)

            
            t+=1
            
            # Plot the signal
            plt.clf()  # Clear the previous plot
            plt.subplot(2,1,1)
            plt.plot(atime,signal_ecg)
            plt.plot(time_peaks,pics_ECG,"*",color='y')
            plt.title("Heart rate by ECG signal:"+str(round(heart_rate_ECG,0))+
                      "\t Heart rate by PPG signal:"+str(round(heart_rate_PPG,0))+
                      "\t The average transit time is:"+str(round(average_arrival_time,3))+
                      "\n The DBP is:"+str(round(DBP,3))+
                      "\n The SBP is:"+str(round(SBP,3)))
            plt.ylabel("ECG")
            plt.xlim(atime[0], atime[display_length-1])  # Set x-axis limits for the last 2000 samples
            plt.draw()

            plt.subplot(2,1,2)
            plt.plot(atime,signal_ppg)
            plt.plot(time_pics_pulso,pics_pulsoximetre,"*",color='y')
            plt.plot(time_plot_max_der,points_max_der,"*",color='g')
            plt.xlabel("Time (seg)")
            plt.ylabel("PPG")
            plt.xlim(atime[0], atime[display_length-1])  # Set x-axis limits for the last 2000 samples
            plt.draw()
            plt.pause(0.1)  # Pause for a short duration to update the plot

            time.sleep(interval)  # Wait for the next segment


            
            


    except KeyboardInterrupt:
        pass  # Stop the loop when KeyboardInterrupt (Ctrl+C) is detected

    finally:
        plt.ioff()  # Turn off interactive mode when done
        plt.show()

# Call the function
continuous_plot(segment_length=10000, interval=10, display_length=20000)


