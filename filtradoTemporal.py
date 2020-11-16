from IPython.display import display
import matplotlib.pyplot as plt
import numpy as np
import os
import shutil
import posixpath

import wfdb

def moving_average(x, w):
    a = np.convolve(x, np.ones(w), 'valid') / w
    b=[a, x[len(a):len(x)-1]]
    return b
    

# Lectura de un registro de WFDB empleado la función 'rdrecord' del objeto wfdb.Record.
record = wfdb.rdrecord('ecg-id-database-1.0.0/Person_01/rec_2', physical=False)
record.record_name = 'persona1Rec2'
# wfdb.plot_wfdb(record=record, title='ecg-id') 
# display(record.__dict__)
# display(record.d_signal)
rawSignal=np.array(record.d_signal[:,0])
# display(rawSignal)
signalMean=moving_average(rawSignal,30)
# signalMean=np.mean(rawSignal)





print(len(signalMean))
print(len(rawSignal))


processSignal=rawSignal-signalMean
fig1, ax1 = plt.subplots(2, 1, figsize=(15, 15), sharex=True)
ax1[0].plot(rawSignal)
ax1[0].set_ylabel('Tensión [mV]', fontsize=15)
ax1[0].grid()
ax1[0].legend(loc="upper right", fontsize=15)
ax1[0].set_title('Señal ECG cruda', fontsize=15)
ax1[1].plot(processSignal)
ax1[1].set_ylabel('Tensión [mV]', fontsize=15)
ax1[1].grid()
ax1[1].legend(loc="upper right", fontsize=15)
ax1[1].set_title('Señal ECG procesada', fontsize=15)
plt.show()


    