from IPython.display import display
import matplotlib.pyplot as plt
from scipy import signal
#%matplotlib inline
import numpy as np
import os
import shutil
import posixpath
import funciones_fft
import wfdb
import time
from FixedPoint import FXfamily, FXnum
 
def filtrar(a,b,signal,filreredSignal,i):
    acum=0
    for j in range(0,len(b),1):
            acum+=b[j]*signal[i-j]
    for j in range(1,len(a),1):
            acum+=a[j]*filreredSignal[i-j]
    return acum
    
#42  - 03 (50Hz) - 27 - 51
# 05 - 85  -89       
# Se carga archivo de la base de datos
record = wfdb.rdrecord('ecg-id-database-1.0.0/Person_42/rec_1', physical=False)
record.record_name = 'persona1Rec2'

# Se transforma a una lista para poder procesarla
rawSignal = list(record.d_signal[:, 0])

fs = 500    #frecuencia de muestreo
ts = 1 / fs                     # tiempo de muestreo
f, fft_senialRaw = funciones_fft.fft_mag(rawSignal, fs)

# # # ############### IIR con µModeler ####################
# # Se extraen los coeficientes de numerador y denominador
# pasa banda orden 1 con corte en 5 y 15Hz para quedarme con el complejo QRS con float  -0.8816185923631891
bp_b=[0.05920834393330736, 0, -0.05920834393330736]
bp_a=[1,1.8704724177979164, -0.8816185923631891] #a[0] (no se usa) lo relleno con uno para que los indices se correspondan con la convención de nombre.
bp_a_den=[1,-1.8704724177979164, 0.8816185923631891]
# # ############### FIR con pyfda ####################
# filtro_fir = np.load('filtroFIR.npz')
# # Se extraen los coeficientes de numerador y denominador
# bp_b, bp_a = filtro_fir['ba']     
# bp_a_den=bp_a*1
# bp_a=[]

################################################ FIX ##########################################################

# pasa banda orden 1 con corte en 5 y 15Hz para quedarme con el aritmatica de punto fijo
# 485, 0, -485, 30646, -14444 DF1
# 31042, 0, -31042, 30646, -14444 DF2
# b(0.0592,0.00,-0.0592) - a (1,1.87,-0.882) #en el vídeo no tenía estos coeficiente, si uso estos coeficientes que µModeler reporta si funciona y responde muy similar.
# ------------------------------------coeficientes de µmodeler:---------------------------------------
    #     short filter1_coefficients[5] = 
    # {
    # // Scaled for a 16x16:64 Direct form 1 IIR filter with: 
    # // Feedback shift = 14
    # // Output shift = 14
    # // Input  has a maximum value of 1, scaled by 2^15
    # // Output has a maximum value of 1.1110299657053153, scaled by 2^14

    #     485, 0, -485, 30646, -14444// b0 Q13(0.0592), b1 Q13(0.00), b2 Q13(-0.0592), a1 Q14(1.87), a2 Q14(-0.882)

    # };
# ---------------------------------------------------------------------------
fixed_b=[0.0592,0.00,-0.0592]
fixed_a=[1,1.87,-0.882] #a[0] (no se usa) lo relleno con cero para que los indices se corresponda con la convención de nombre.
intbits = np.ceil(np.log2(max(max(abs(np.array(bp_a))), max(abs(np.array(bp_b))))))
fm = FXfamily(16,intbits)
FXnumv = np.vectorize(FXnum)
floatv = np.vectorize(float)
bp_b = floatv(FXnumv(fixed_b,fm))
bp_a = floatv(FXnumv(fixed_a,fm))
bp_a_den=bp_a*-1
bp_a_den[0]=1


######################## RESPUESTA DEL FILTRO ##########################
f_iir, h_iir = signal.freqz(bp_b, bp_a_den, worN=f, fs=fs)




# stop band de 49 a 51Hz
sp_b=[0.9875889380903246, -1.5980786463154453, 0.9875889380903248]
sp_a=[1,1.598078646315445, -0.9751778761806491] 


aux = 0
maximo = 0
tamanio = range(len(rawSignal))
deteccion = []
stop_detection =False
filSignal=list(np.ones(len(bp_b)-1,))
filtrada=list(np.ones(len(bp_b)-1))
filtrada2=list(np.ones(len(bp_b)+2))
periodo=[]
variabilidad=[]
contador=0
refractario=0
compensacion=0
for i in tamanio:
    start=time.time()
    refractario+=1
    if i> len(bp_b)-2:
        filtrada.append(filtrar(bp_a,bp_b,rawSignal,filSignal,i))# guardo copia solo filtrada para graficarla
        filSignal.append(filtrar(bp_a,bp_b,rawSignal,filSignal,i)) # Esta es la señal que se procesa para detectar el R
        if i>len(bp_b)+1:
            filtrada2.append(filtrar(sp_a,sp_b,filtrada,filtrada2,i))# filtrado rechaza banda 50Hz (no fue util y se dejo este paso para mostrar la gráfica)
        filSignal[i-3]=(np.power(filSignal[i-3],3))#elevo para realzar R
        filSignal[i-3] = filSignal[i-3] - 3800 #paso a negativo el ruido cerca de la linea de base
        if filSignal[i-3]<0:
            filSignal[i-3]=0
        if refractario>90:
            stop_detection = False
        else:
            stop_detection = True
        if filSignal[i-3] != 0 and stop_detection==False:
            # contador=0
            aux = filSignal[i-3]
            if aux > maximo:
                maximo = aux
                evento = i-6
         
        if aux < maximo and maximo!=0:            
        # if filSignal[i-3] == 0 and maximo!=0:
            contador+=1
            if contador==220:
                contador=0
                deteccion.append(evento)
                refractario=0
                maximo = 0
                if len(deteccion)>1:
                    periodo.append(deteccion[len(deteccion)-1]-deteccion[len(deteccion)-2])
                if len(periodo)>1:
                    variabilidad.append((periodo[len(periodo)-1]-periodo[len(periodo)-2])/fs*1000)# Variabilidad en ms
                    print('Hay un latido')
                    print((variabilidad[len(variabilidad)-1]),'[ms]')
        demora=time.time()-start
        #######simulación de espera del período de muestreo##############
        # if demora>ts:
        #     compensacion+=demora
        #     # print(demora)
        # else:
        #     if (demora+compensacion)>ts:
        #         compensacion=demora+compensacion-ts
        #     else:
        #         time.sleep(ts-demora-compensacion) 
        #         compensacion=0
    
    
        



print(np.array(deteccion))

f, fft_senialFil2 = funciones_fft.fft_mag(filtrada2, fs)
f, fft_senialFil = funciones_fft.fft_mag(filtrada, fs)

fig, ax1 = plt.subplots(1, 1, figsize=(15, 15), sharex=True)
ax1.scatter(deteccion,np.array(rawSignal)[deteccion], label='detecciones', color='orange')
ax1.plot(rawSignal)




fig, ax2 = plt.subplots(4, 1, figsize=(15, 15), sharex=True)
ax2[0].plot(filSignal, label='señal filtrada', color='orange')
ax2[1].plot(filtrada)
ax2[2].plot(rawSignal)
ax2[3].plot(filtrada2)

fig, ax2 = plt.subplots(3, 1, figsize=(15, 15), sharex=True)
ax2[0].plot(f, fft_senialRaw)
factorDeEscala=np.max(abs(fft_senialRaw))/np.max(abs(h_iir))
ax2[0].plot(f_iir, factorDeEscala*abs(h_iir), label='Filtro', color='orange')
ax2[1].plot(f,fft_senialFil, label='señal filtrada', color='orange')
ax2[2].plot(f, fft_senialFil2)



plt.show()



