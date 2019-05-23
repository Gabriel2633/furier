#Ejercicio4
# 'y' es una senal en funcion del tiempo 't' con las unidades descritas en el codigo.
# a. Grafique la senal en funcion del tiempo en la figura 'senal.png' ('y' vs. 't')
# b. Calule la transformada de Fourier (sin utilizar funciones de fast fourier transform) y
# grafique la norma de la transformada en funcion de la frecuencia (figura 'fourier.png')
# c. Lleve a cero los coeficientes de Fourier con frecuencias mayores que 10000 Hz y calcule 
# la transformada inversa para graficar la nueva senal (figura 'filtro.png')

import numpy as np
import matplotlib.pyplot as plt

#DAdo
n = 512 # number of point in the whole interval
f = 200.0 #  frequency in Hz
dt = 1 / (f * 128 ) #128 samples per frequency
t = np.linspace( 0, (n-1)*dt, n) 
y = np.sin(2 * np.pi * f * t) + np.cos(2 * np.pi * f * t * t)
noise = 1.4*(np.random.rand(n)+0.7)
y = y + noise
#DAdo


#Transdormada de Fourier
def TF(x):
    N = len(x) #Numero de datos
    C = np.zeros(N, dtype=complex) #Coeficientes de Fourier
    W = np.zeros(N) #Frecuencias

    for k in range(N): #Itera sobre las frecuencias
        C[k] = 0.0j
        
        if k<N/2: #
            W[k]=f*k
        else: #
            W[k]=f*(-n+k)
        
        for n in range(N): #Calcula el Coeficiente como la suma de la integral
            C[k] += x[n] * np.exp(-2.0 * np.pi * 1.0j / N ) ** (k * n) 
        
    return W,C

def IFT(x):
    N = len(x)
    X = np.zeros(N, dtype=complex)

    for k in range(N):
        X[k] = 0.0j
        for n in range(N):
            X[k] += x[n] * np.exp(2.0 * np.pi * 1.0j / N ) ** (k * n) 
        
    return X/N

w,y_w=TF(y) #calcula la transformada
y_wf=y_w.copy()
y_wf[abs(w)>10000]=0 #filtra las frecuwencias
yf=IFT(y_wf) #calcula la transdormada inversa de los coeficientes filtrados

plt.figure()
plt.plot(t,y)
plt.show()

plt.figure()
plt.scatter(w,abs(y_w))
plt.show()

plt.figure()
plt.scatter(w,abs(y_wf))
plt.show()

plt.figure()
plt.plot(t,yf)
plt.show()