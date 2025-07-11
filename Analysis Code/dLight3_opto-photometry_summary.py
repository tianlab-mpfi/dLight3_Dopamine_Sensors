# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 14:38:31 2025

@author: kenta.hagihara
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['pdf.use14corefonts'] = False
mpl.rcParams['font.family'] = 'Arial'

datadir='...\Analysis Data\Figure4_opto-photometry\'

#3.6
data0 = np.load(datadir + "765713_o_15000_40_trials.npy")
data1 = np.load(datadir + r"765713_p_15000_40_trials.npy")
data2 = np.load(datadir + r"765713_q_15000_40_trials.npy")
data3 = np.load(datadir + r"775991_o_15000_40_trials.npy")
data4 = np.load(datadir + r"775991_p_15000_40_trials.npy")
data5 = np.load(datadir + r"775991_q_15000_40_trials.npy")
data6 = np.load(datadir + r"771374_o_15000_40_trials.npy")
data7 = np.load(datadir + r"771374_p_15000_40_trials.npy")
data8 = np.load(datadir + r"771374_q_15000_40_trials.npy")
data9 = np.load(datadir + r"771376_o_15000_40_trials.npy")
data10 = np.load(datadir + r"771376_p_15000_40_trials.npy")
data11 = np.load(datadir + r"771376_q_15000_40_trials.npy")
data12 = np.load(datadir + r"765716_o_15000_40_trials.npy")
data13 = np.load(datadir + r"765716_p_15000_40_trials.npy")
data14 = np.load(datadir + r"765716_q_15000_40_trials.npy")
data15 = np.load(datadir + r"765718_o_15000_40_trials.npy")
data16 = np.load(datadir + r"765718_p_15000_40_trials.npy")
data17 = np.load(datadir + r"765718_q_15000_40_trials.npy")
data18 = np.load(datadir + r"775993_o_15000_40_trials.npy")
data19 = np.load(datadir + r"775993_p_15000_40_trials.npy")
data20 = np.load(datadir + r"775993_q_15000_40_trials.npy")

#DA3m
data21 = np.load(datadir + r"765648_o_15000_40_trials.npy")
data22 = np.load(datadir + r"765648_p_15000_40_trials.npy")
data23 = np.load(datadir + r"765648_q_15000_40_trials.npy")
data24 = np.load(datadir + r"765651_o_15000_40_trials.npy")
data25 = np.load(datadir + r"765651_p_15000_40_trials.npy")
data26 = np.load(datadir + r"765651_q_15000_40_trials.npy")
data27 = np.load(datadir + r"765649_o_15000_40_trials.npy")
data28 = np.load(datadir + r"765649_p_15000_40_trials.npy")
data29 = np.load(datadir + r"765649_q_15000_40_trials.npy")
data30 = np.load(datadir + r"765650_o_15000_40_trials.npy")
data31 = np.load(datadir + r"765650_p_15000_40_trials.npy")
data32 = np.load(datadir + r"765650_q_15000_40_trials.npy")
data33 = np.load(datadir + r"777479_o_15000_40_trials.npy")
data34 = np.load(datadir + r"777479_p_15000_40_trials.npy")
data35 = np.load(datadir + r"777479_q_15000_40_trials.npy")
data36 = np.load(datadir + r"777832_o_15000_40_trials.npy")
data37 = np.load(datadir + r"777832_p_15000_40_trials.npy")
data38 = np.load(datadir + r"777832_q_15000_40_trials.npy")
data39 = np.load(datadir + r"777833_o_15000_40_trials.npy")
data40 = np.load(datadir + r"777833_p_15000_40_trials.npy")
data41 = np.load(datadir + r"777833_q_15000_40_trials.npy")

#dLight3.8
data42 = np.load(datadir + r"F:\Analysis\dLight3_paper2\765720_o_15000_40_trials.npy")
data43 = np.load(datadir + r"F:\Analysis\dLight3_paper2\765720_p_15000_40_trials.npy")
data44 = np.load(datadir + r"F:\Analysis\dLight3_paper2\765720_q_15000_40_trials.npy")
data45 = np.load(datadir + r"F:\Analysis\dLight3_paper2\765721_o_15000_40_trials.npy")
data46 = np.load(datadir + r"F:\Analysis\dLight3_paper2\765721_p_15000_40_trials.npy")
data47 = np.load(datadir + r"F:\Analysis\dLight3_paper2\765721_q_15000_40_trials.npy")
data48 = np.load(datadir + r"F:\Analysis\dLight3_paper2\768978_o_15000_40_trials.npy")
data49 = np.load(datadir + r"F:\Analysis\dLight3_paper2\768978_p_15000_40_trials.npy")
data50 = np.load(datadir + r"F:\Analysis\dLight3_paper2\768978_q_15000_40_trials.npy")

preW=100
sampling_rate=20

def PSTHplot(PSTH, MainColor, SubColor, LabelStr):
    plt.plot(np.arange(np.shape(PSTH)[1])/20 - preW/sampling_rate, np.mean(PSTH.T,axis=1),label=LabelStr,color = MainColor)
    y11 =  np.mean(PSTH.T,axis=1) + np.std(PSTH.T,axis=1)/np.sqrt(np.shape(PSTH)[0])
    y22 =  np.mean(PSTH.T,axis=1) - np.std(PSTH.T,axis=1)/np.sqrt(np.shape(PSTH)[0])
    plt.fill_between(np.arange(np.shape(PSTH)[1])/20 - preW/sampling_rate, y11, y22, facecolor=SubColor, alpha=0.5)

#%% dLight3.6
Ave_36_10Hz=np.empty((7,400))
Ave_36_20Hz=np.empty((7,400))
Ave_36_5Hz=np.empty((7,400))

plt.figure(figsize=(8, 10))

for ii in range(7):
    LabelName = "3.6_animal_" + str(ii+1)
    plt.subplot(4,2,ii+1)
    PSTHplot(eval('data' + str(ii*3+0)), "orange", "darkorange", "10Hz")
    PSTHplot(eval('data' + str(ii*3+1)), "g", "darkgreen", "20Hz")
    PSTHplot(eval('data' + str(ii*3+2)), "magenta", "darkmagenta", "5Hz") 
    plt.xlabel('Time - Optostim (s)')
    plt.title(LabelName)
    plt.axvspan(0, 2, color = [0, 0, 0, 0.1])
    plt.ylim([-10, 250])

    Ave_36_10Hz[ii,:] = np.mean(eval('data' + str(ii*3+0)),axis=0)
    Ave_36_20Hz[ii,:] = np.mean(eval('data' + str(ii*3+1)),axis=0)
    Ave_36_5Hz[ii,:] = np.mean(eval('data' + str(ii*3+2)),axis=0)
    
plt.subplot(4,2,8)
PSTHplot(Ave_36_10Hz, "orange", "darkorange", "10Hz")
PSTHplot(Ave_36_20Hz, "g", "darkgreen", "20Hz")
PSTHplot(Ave_36_5Hz, "magenta", "darkmagenta", "5Hz") 
plt.ylim([-10, 250])
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.xlabel('Time - Optostim (s)')
plt.title("dLight3.6 summary")
plt.axvspan(0, 2, color = [0, 0, 0, 0.1])
plt.legend()

#%% DA3m
Ave_DA3m_10Hz=np.empty((7,400))
Ave_DA3m_20Hz=np.empty((7,400))
Ave_DA3m_5Hz=np.empty((7,400))

plt.figure(figsize=(8, 10))

for ii in range(7):
    LabelName = "GRAB-DA3m_animal_" + str(ii+1)
    plt.subplot(4,2,ii+1)
    PSTHplot(eval('data' + str(ii*3+21)), "orange", "darkorange", "10Hz")
    PSTHplot(eval('data' + str(ii*3+22)), "g", "darkgreen", "20Hz")
    PSTHplot(eval('data' + str(ii*3+23)), "magenta", "darkmagenta", "5Hz") 
    plt.xlabel('Time - Optostim (s)')
    plt.title(LabelName)
    plt.axvspan(0, 2, color = [0, 0, 0, 0.1])
    plt.ylim([-10, 250])
    plt.subplots_adjust(wspace=0.5, hspace=0.5)

    Ave_DA3m_10Hz[ii,:] = np.mean(eval('data' + str(ii*3+21)),axis=0)
    Ave_DA3m_20Hz[ii,:] = np.mean(eval('data' + str(ii*3+22)),axis=0)
    Ave_DA3m_5Hz[ii,:] = np.mean(eval('data' + str(ii*3+23)),axis=0)

plt.subplot(4,2,8)
PSTHplot(Ave_DA3m_10Hz, "orange", "darkorange", "10Hz")
PSTHplot(Ave_DA3m_20Hz, "g", "darkgreen", "20Hz")
PSTHplot(Ave_DA3m_5Hz, "magenta", "darkmagenta", "5Hz") 
plt.ylim([-10, 250])
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.xlabel('Time - Optostim (s)')
plt.title("DA3m summary")
plt.axvspan(0, 2, color = [0, 0, 0, 0.1])
plt.legend()

#%% dLight3.8 
Ave_38_10Hz=np.empty((3,400))
Ave_38_20Hz=np.empty((3,400))
Ave_38_5Hz=np.empty((3,400))

plt.figure(figsize=(8, 10))

for ii in range(3):
    LabelName = "dLight3.8_animal_" + str(ii+1)
    plt.subplot(4,2,ii+1)
    PSTHplot(eval('data' + str(ii*3+42)), "orange", "darkorange", "10Hz")
    PSTHplot(eval('data' + str(ii*3+43)), "g", "darkgreen", "20Hz")
    PSTHplot(eval('data' + str(ii*3+44)), "magenta", "darkmagenta", "5Hz") 
    plt.xlabel('Time - Optostim (s)')
    plt.title(LabelName)
    plt.axvspan(0, 2, color = [0, 0, 0, 0.1])
    plt.ylim([-10, 250])
    plt.subplots_adjust(wspace=0.5, hspace=0.5)
    
    Ave_38_10Hz[ii,:] = np.mean(eval('data' + str(ii*3+42)),axis=0)
    Ave_38_20Hz[ii,:] = np.mean(eval('data' + str(ii*3+43)),axis=0)
    Ave_38_5Hz[ii,:] = np.mean(eval('data' + str(ii*3+44)),axis=0)

plt.subplot(4,2,8)
PSTHplot(Ave_38_10Hz, "orange", "darkorange", "10Hz")
PSTHplot(Ave_38_20Hz, "g", "darkgreen", "20Hz")
PSTHplot(Ave_38_5Hz, "magenta", "darkmagenta", "50Hz") 
plt.ylim([-10, 250])
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.xlabel('Time - Optostim (s)')
plt.title("dLight3.8 summary")
plt.axvspan(0, 2, color = [0, 0, 0, 0.1])
plt.legend()

#%% Summary
plt.figure(figsize=(24,6)),

plt.subplot(1,6,1)
PSTHplot(Ave_DA3m_20Hz, "g", "darkgreen", "20Hz")
PSTHplot(Ave_DA3m_10Hz, "orange", "darkorange", "10Hz")
PSTHplot(Ave_DA3m_5Hz, "magenta", "darkmagenta", "5Hz") 
plt.ylim([-10, 300])
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.xlabel('Time - Optostim (s)')
plt.title("DA3m summary")
plt.axvspan(0, 2, color = [0, 0, 0, 0.1])
plt.legend()

plt.subplot(1,6,2)
PSTHplot(Ave_36_20Hz, "g", "darkgreen", "20Hz")
PSTHplot(Ave_36_10Hz, "orange", "darkorange", "10Hz")
PSTHplot(Ave_36_5Hz, "magenta", "darkmagenta", "5Hz") 
plt.ylim([-10, 300])
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.xlabel('Time - Optostim (s)')
plt.title("dLight3.6 summary")
plt.axvspan(0, 2, color = [0, 0, 0, 0.1])
plt.legend()

plt.subplot(1,6,3)
PSTHplot(Ave_38_20Hz, "g", "darkgreen", "20Hz")
PSTHplot(Ave_38_10Hz, "orange", "darkorange", "10Hz")
PSTHplot(Ave_38_5Hz, "magenta", "darkmagenta", "50Hz") 
plt.ylim([-10, 300])
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.xlabel('Time - Optostim (s)')
plt.title("dLight3.8 summary")
plt.axvspan(0, 2, color = [0, 0, 0, 0.1])
plt.legend()

#%% Peak amp
plt.subplot(1,6,4)

max36_5=np.max(Ave_36_5Hz,axis=1)
max36_10=np.max(Ave_36_10Hz,axis=1)
max36_20=np.max(Ave_36_20Hz,axis=1)
maxDA3m_5=np.max(Ave_DA3m_5Hz,axis=1)
maxDA3m_10=np.max(Ave_DA3m_10Hz,axis=1)
maxDA3m_20=np.max(Ave_DA3m_20Hz,axis=1)
max38_5=np.max(Ave_38_5Hz,axis=1)
max38_10=np.max(Ave_38_10Hz,axis=1)
max38_20=np.max(Ave_38_20Hz,axis=1)


plt.plot(np.tile([0.4, 0.8, 1.6],(7,1)).T, [maxDA3m_5, maxDA3m_10, maxDA3m_20],color=[0.5, 0.5, 0.5, 0.5])
plt.errorbar(0.4, np.mean(maxDA3m_5),np.std(maxDA3m_5)/np.sqrt(len(maxDA3m_5)), fmt='o', color="magenta")
plt.errorbar(0.8, np.mean(maxDA3m_10),np.std(maxDA3m_10)/np.sqrt(len(maxDA3m_10)), fmt='o', color="orange")
plt.errorbar(1.6, np.mean(maxDA3m_20),np.std(maxDA3m_20)/np.sqrt(len(maxDA3m_20)), fmt='o', color="g")

plt.plot(np.tile([2.4, 2.8, 3.6],(7,1)).T, [max36_5, max36_10, max36_20],color=[0.5, 0.5, 0.5, 0.5])
plt.errorbar(2.4, np.mean(max36_5),np.std(max36_5)/np.sqrt(len(max36_5)), fmt='o', color="magenta")
plt.errorbar(2.8, np.mean(max36_10),np.std(max36_10)/np.sqrt(len(max36_10)), fmt='o', color="orange")
plt.errorbar(3.6, np.mean(max36_20),np.std(max36_20)/np.sqrt(len(max36_20)), fmt='o', color="g")

plt.plot(np.tile([4.4, 4.8, 5.6],(3,1)).T, [max38_5, max38_10, max38_20],color=[0.5, 0.5, 0.5, 0.5])
plt.errorbar(4.4, np.mean(max38_5),np.std(max38_5)/np.sqrt(len(max38_5)), fmt='o', color="magenta")
plt.errorbar(4.8, np.mean(max38_10),np.std(max38_10)/np.sqrt(len(max38_10)), fmt='o', color="orange")
plt.errorbar(5.6, np.mean(max38_20),np.std(max38_20)/np.sqrt(len(max38_20)), fmt='o', color="g")
plt.xticks([1,3,5],["DA3m", "3.6", "3.8"])
plt.ylabel('Peak dF/F (%)')


#%% tau-off
from scipy.optimize import curve_fit
def exponential_func(x, a, b, c):
    return -a * np.exp(-b * x) + c

# Generate some data
tau_36=[]
tau_DA3m=[]
tau_38=[]

tauon_36=[]
tauon_DA3m=[]
tauon_38=[]

x_data_d = np.arange(0, 7, 1/20)
x_data_d_on = np.arange(0, 2, 1/20)

for ii in range(7):
    y_data_d = Ave_36_20Hz[ii, 140:280]
    y_data_d_on = Ave_36_20Hz[ii, 100:140]
    popt_d, pcov_d = curve_fit(exponential_func, x_data_d, y_data_d)
    popt_d_on, pcov_d = curve_fit(exponential_func, x_data_d_on, y_data_d_on)
    tau_36.append(np.log(2)/popt_d[1])
    tauon_36.append(np.log(2)/popt_d_on[1])

for ii in range(7):
    y_data_d = Ave_DA3m_20Hz[ii, 140:280]
    y_data_d_on = Ave_DA3m_20Hz[ii, 100:140]
    popt_d, pcov_d = curve_fit(exponential_func, x_data_d, y_data_d)
    popt_d_on, pcov_d = curve_fit(exponential_func, x_data_d_on, y_data_d_on)
    tau_DA3m.append(np.log(2)/popt_d[1])
    tauon_DA3m.append(np.log(2)/popt_d_on[1])

for ii in range(3):
    y_data_d = Ave_38_20Hz[ii, 140:280]
    y_data_d_on = Ave_38_20Hz[ii, 100:140]
    popt_d, pcov_d = curve_fit(exponential_func, x_data_d, y_data_d)
    popt_d_on, pcov_d = curve_fit(exponential_func, x_data_d_on, y_data_d_on)
    tau_38.append(np.log(2)/popt_d[1])
    tauon_38.append(np.log(2)/popt_d_on[1])

plt.subplot(1,6,6)
plt.plot([1] * 7, tau_DA3m, "o", color=[0.5, 0.5, 0.5])
plt.errorbar(1, np.mean(tau_DA3m),np.std(tau_DA3m)/np.sqrt(len(tau_DA3m)), fmt='o', color="green")

plt.plot([2] * 7, tau_36, "o", color=[0.5, 0.5, 0.5])
plt.errorbar(2, np.mean(tau_36),np.std(tau_36)/np.sqrt(len(tau_36)), fmt='o', color="green")

plt.plot([3] * 3, tau_38, "o", color=[0.5, 0.5, 0.5])
plt.errorbar(3, np.mean(tau_38),np.std(tau_38)/np.sqrt(len(tau_38)), fmt='o', color="green")

plt.ylim([0, 3])
plt.xlim([0, 4])
plt.xticks([1,2,3],["DA3m", "3.6", "3.8"])
plt.ylabel('Tau-off (s)')

plt.subplot(1,6,5)
plt.plot([1] * 7, tauon_DA3m, "o", color=[0.5, 0.5, 0.5])
plt.errorbar(1, np.mean(tauon_DA3m),np.std(tauon_DA3m)/np.sqrt(len(tauon_DA3m)), fmt='o', color="green")

plt.plot([2] * 7, tauon_36, "o", color=[0.5, 0.5, 0.5])
plt.errorbar(2, np.mean(tauon_36),np.std(tauon_36)/np.sqrt(len(tauon_36)), fmt='o', color="green")

plt.plot([3] * 3, tauon_38, "o", color=[0.5, 0.5, 0.5])
plt.errorbar(3, np.mean(tauon_38),np.std(tauon_38)/np.sqrt(len(tauon_38)), fmt='o', color="green")

plt.ylim([0, 1])
plt.xlim([0, 4])
plt.xticks([1,2,3],["DA3m", "3.6", "3.8"])
plt.ylabel('Tau-on (s)')

#%% stat test
from scipy.stats import ttest_ind

tt,pp = ttest_ind(tauon_DA3m, tauon_36)
tt,pp_off = ttest_ind(tau_DA3m, tau_36)

print("Tau-on (s), DA3m:" + str(np.mean(tauon_DA3m)) + " +- " + str(np.std(tauon_DA3m)/np.sqrt(len(tauon_DA3m))))
print("Tau-on (s), 36:" + str(np.mean(tauon_36)) + " +- " + str(np.std(tauon_36)/np.sqrt(len(tauon_36))))
print("Tau-on (s), 38:" + str(np.mean(tauon_38)) + " +- " + str(np.std(tauon_38)/np.sqrt(len(tauon_38))))
print("Tau-on (s), DA3m vs 3.6, p = " + str(pp))
print("Tau-off (s), DA3m vs 3.6, p = " + str(pp_off))

