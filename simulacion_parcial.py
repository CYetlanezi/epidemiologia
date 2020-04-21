import scipy.integrate as spi
import numpy as np
import pylab as pl

betaY=; #Tasa de contacto efectiva de subpoblaciones sintomaticas a susceptibles
phi=0.55; #Probabilidad de transmision relativa
betaA=phi*betaY; #Tasa de contacto efectivo de subpoblacion asintomatica a susceptible
gammaINV=5.1; #Periodo de incubacion
alpha=0.86; #Probabilidad de volverse asintomatico tras la infeccion
lambdaYRinv= np.array([8, 37]); #Tiempo promedio de infeccion de sintomaticos
lambdaARinv= np.array([8, 37]); #Tiempo promedio de infeccion de asintomaticos
lambdaAY=;#transicion sintomatica
lmbdaRS=;#Tasa de recaída
delta=0.032*(1-alpha);#Tasa de mortalidad inducida por la enfermedad
Gamma=; #parámetro demografico
Xi=;#parámetro demografico
S0=;
E0=;
I0=;
Y0=;
A0=;
R0=1-S0-E0-Y0-A0
tV=30*365;
ND=MaxTime=100*365;
TS=1.0


INPUT = np.hstack((S0,E0,I0,Y0,A0,R0))

def diff_eqs(INP,t):  
	'''ecuaciones principales'''
	Y=np.zeros((4))
	V = INP   
	Y[0]= Gamma+(lambdaRS)V[5]-(betaY*(V[3]/N)+betaA*(V[4]/N)+Xi)*V[0] #S
	Y[1]=(betaY*(V[3]/N)+betaA*(V[4]/N)+Xi)*V[0]-(1/gammaINV+Xi)*V[1] #E
	Y[2]=(1/gammaINV)*(1-alpha)*V[1]-(Xi+delta+1/lambdaYRinv)*V[3]+lambdaAY*V[4]#I
    Y[3]=(1/gammaINV)*alpha*V[1]-(1/lambdaYRinv + lambdaAY + Xi)*V[4] #A
    Y[4]=(1/lambdaAYinv)*V[4]+(1/lambdaYRinv)*V[3]-(lambdaRS+Xi)V[5] #R
	return Y   # For odeint

t_start = 0.0; t_end = tV; t_inc = TS
t_range1 = np.arange(t_start, t_end+t_inc, t_inc)
t_start = tV; t_end = ND; t_inc = TS
t_range2 = np.arange(tV, t_end+t_inc, t_inc)
T = np.hstack((t_range1, t_range2))
v=0
RES1 = spi.odeint(diff_eqs,INPUT,t_range1)
v=v0
RES2 = spi.odeint(diff_eqs,RES1[-1],t_range2)
print (RES2)

S = np.hstack((RES1[:,0],RES2[:,0]))
I = np.hstack((RES1[:,1],RES2[:,1]))
R = np.hstack((RES1[:,2],RES2[:,2]))

pl.subplot(311)
pl.plot(T[1:,]/365.0, S[1:,], '-g')
ll=pl.ylim()
tVV=np.repeat([tV/365.],len(ll))
pl.plot(tVV, ll, '--k')
pl.ylabel('Susceptibles')
pl.subplot(312)
pl.plot(T[1:,]/365.0, I[1:,], '-r')
ll=pl.ylim()
tVV=np.repeat([tV/365.],len(ll))
pl.plot(tVV, ll, '--k')
pl.ylabel('Infectious')
pl.subplot(313)
pl.plot(T[1:,]/365.0, R[1:,], '-k')
ll=pl.ylim()
tVV=np.repeat([tV/365.],len(ll))
pl.plot(tVV, ll, '--k')
pl.ylabel('Recovereds/Vaccinated')
pl.xlabel('Time (Years)')

pl.show()



