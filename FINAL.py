# -*- coding: utf-8 -*-
"""
Created on Sun Oct 25 16:33:15 2020

@author: igna
intento en python del tp de computacional para el final

"""
import numpy as np
import matplotlib.pyplot as plt
import numba
import time 

# funciones exactas
def Energia_exacta(R, m):
    rho_0 = 0.037
    p_f = (6*rho_0*R*np.pi**2)**(1/3)
    e_f = p_f**2/(2*m)
    tau = 0.05
    resultado = (3/5)*e_f*(1 + (5/12)*(np.pi**2)*(tau**2))
    return resultado


def Dist_momentos_exacta(q):
    tau = 0.05
    #p_f = (6*rho*np.pi**2)**(1/3)
    #p = 2*(np.random.rand()-0.5)*p_f/np.cbrt(3)
    #q = p/p_f
    numerador = 3*q**2
    denominador = np.exp((1/tau)*(q**2 -1 + (np.pi**2*tau**2/12))) + 1
    return numerador/denominador


@numba.jit(nopython = True)
def p_Fermi(rho):
    return (6*rho*np.pi**2)**(1/3)


def distribuir_xv(rho, N):
    # calculo el tamaño de la caja como función de rho y N
    L = (N/rho)**(1/3)
    # calculo el impulso de Fermi
    p_f = p_Fermi(rho)
    # reparto las partículas de forma random
    x = np.random.rand(N, 3)*L
    v = 2*(np.random.rand(N, 3)-0.5)*p_f
    return x, v


@numba.jit(nopython = True)
def modulo_p(m, v):
    p2 = (m*v)**2
    p = (np.sum(p2, axis = 1))**(1/2)
    return p
    

@numba.jit(nopython = True)
def cinetica(m, p, N):
    """
    Calcula la energía cinética libre por partícula
    """
    e_cin = 1/(2*m)*np.sum(p**2)
    return e_cin/N
    

@numba.jit(nopython = True)
def theta(rho, p, N):
    """
    Calcula la contribución del término de theta en el
    pseudopotencial de Pauli, por partícula
    """
    # parámetros del paper
    rho_0 = 0.037 # parámetro dado en el paper y fijo
    alpha_c = 0.831 # parámetro dado en el paper y fijo
    v_th = 3.560 # parámetro dado en el paper y fijo
    eta = 30 # parámetro dado en el paper y fijo
    # calculo las magnitudes que necesito
    p_f = p_Fermi(rho)
    q = p/p_f
    VC = v_th*(rho/rho_0)**alpha_c
    denominador = 1 + np.exp(-eta*(q**2 - 1))
    Theta = VC/denominador
    V_theta = np.sum(Theta)
    return V_theta/N


@numba.jit(nopython = True)  # para acelerar un poco las cosas
def Vpauli(rho, x, v, N):
    """
    Calcula el valor del "potencial" (ficticio) para emular la
    repulsion de Pauli, por partícula.
    """
    # defino constantes
    rho_0 = 0.037 # parámetro dado en el paper y fijo
    Vq = 13.517 # parámetro dado en el paper y fijo
    Vp = 1.260 # parámetro dado en el paper y fijo
    r0 = 0.845/p_Fermi(rho)
    p0 = 0.193*p_Fermi(rho)
    alpha_a = 0.629 # parámetro dado en el paper y fijo
    alpha_b = 0.665 # parámetro dado en el paper y fijo
    VA = Vq*(rho/rho_0)**alpha_a
    VB = Vp*(rho/rho_0)**alpha_b
    # inicializo la energia
    r_ij = 0.0
    v_ij = 0.0
    E_potencial = 0.0
    for i in range(N):
        for j in range(i+1,N):
            for k in range(3):
                # calculo las distancias entre partículas (tanto x como v)
                r_ij += (x[i, k] - x[j, k])**2
                v_ij += (v[i, k] - v[j, k])**2
            r_ij = r_ij**(1/2)
            v_ij = v_ij**(1/2)
            # uso eso para calcular el pseudopotencial de exclusion
            # reinicio las variables
            E_potencial += VA*np.exp(-r_ij/r0) + VB*np.exp(-v_ij/p0)
            r_ij = 0.0
            v_ij = 0.0
    return E_potencial/N


@numba.jit(nopython = True)
def energia(rho, m, x, v, p, N):
    """
    Suma las contribuciones a la energia de
    Energia partícula libre por partícula
    Energia del "potencial" theta por partícula
    Energia del "potencial" de repulsión de pauli por partícula
    """
    E_cin = cinetica(m, p, N)
    V_th = theta(rho, p, N)
    V_pauli = Vpauli(rho, x, v, N)
    Energia = E_cin + V_th + V_pauli
    return Energia


"""
def montecarlo(rho, m, x, v, p, N):
    # calculo la energía inicial, esto se hace por unica vez al correr esta
    # funcion
    E_inicial = energia(rho, m, x, v, p, N)
    # calculo y defino algunos parámetros
    p_f = p_Fermi(rho)
    T_f = p_f**2/(2*m)
    T = 0.05*T_f
    beta = 1/T
    L = (N/rho)**(1/3)
    # defino el salto random máximo para las coordenadas x y v
    salto_x = L
    salto_v = p_f/m
    # defino los vectores a modificar
    x_modif = x
    v_modif = v
    p_modif = p
    # loop principal
    for i in range(N):
        # modifico las coordenadas de 1 partícula
        x_modif[i,:] = x[i,:] + (np.random.rand(3)-0.5)*salto_x
        v_modif[i,:] = v[i,:] + (np.random.rand(3)-0.5)*salto_v
        p_modif = modulo_p(m, v_modif)
        # chequeo condiciones de contorno periodicas
        x_modif[np.where(x_modif < 0)] = x_modif[np.where(x_modif < 0)] + L
        x_modif[np.where(x_modif > L)] = x_modif[np.where(x_modif > L)] - L
        # calculo la energía una vez que modifiqué la primer partícula
        E_final = energia(rho, m, x_modif, v_modif, p_modif, N)
        Delta_E = E_final - E_inicial
        if Delta_E > 0:
            # defino la probabilidad de aceptación
            p_aceptacion = np.exp(-Delta_E*beta)
            rand = np.random.rand()    
            if p_aceptacion < rand:
                # si no acepto el cambio, entonces vuelvo atrás
                # la modificación al vector y además redefino
                # la energía final como la inicial (porque no acepté)
                x_modif[i,:] = x[i,:]
                v_modif[i,:] = v[i,:]
                p_modif[i] = p[i]
                E_final = E_inicial
    return E_final, x_modif, v_modif, p_modif
"""


def montecarlojit(rho, m, x, v, p, N):
    """
    Lleva a cabo el algoritmo de montecarlo:
        Calcula la energía inicial y la guarda como E_i
        Luego perturba una sola partícula de todo el sistema y vuelve
        a calcular la energía, la guarda como E_f
        Mide la diferencia de energia Delta_E = E_f - E_i:
            Si Delta_E < 0, se acepta el estado (es decir, no se hace
            nada nuevo porque ya modifiqué la posición y momento de 
            la partícula además de su energia)
            Si Delta_E > 0, se arma la probabilidad de aceptación y se
            la compara con un numero random:
                Si p_acep > rand, se acepta el estado (es decir, no se
                hace nada nuevo)
                Si p_acep < rand, no se acepta el estado, entonces hay 
                que volver al estado inicial. Esto es borrar la
                modificación a las coordenadas de posición y momento.
    return: E_final, x, v, p
    """
    # calculo la energía inicial, esto se hace por unica vez al correr esta
    # funcion
    E_inicial = energia(rho, m, x, v, p, N)
    # calculo y defino algunos parámetros
    p_f = p_Fermi(rho)
    T_f = (p_f**2)/(2*m)
    T = 0.05*T_f # segun paper, T/T_f = 0.05 fijo
    beta = 1/T
    L = (N/rho)**(1/3)
    # defino el salto random máximo para las coordenadas x y v
    salto_x = L
    salto_v = p_f/m
    # defino los vectores a modificar
    # los igual a su np.array() (lo cual NO ES REDUNDANTE) para evitar
    # el aliasing
    x_modif = np.array(x)
    v_modif = np.array(v)
    p_modif = np.array(p)
    # loop principal
    for i in range(N):
        for k in range(3):
            # modifico las coordenadas de 1 partícula
            x_modif[i, k] = x[i, k] + (np.random.rand()-0.5)*salto_x
            v_modif[i, k] = v[i, k] + (np.random.rand()-0.5)*salto_v
            # chequeo condiciones de contorno periodicas
            if x_modif[i, k] < 0:
                x_modif[i, k] = x_modif[i, k] + L
            elif x_modif[i, k] > L:
                x_modif[i, k] = x_modif[i, k] - L
        p_modif = modulo_p(m, v_modif)
        # calculo la energía una vez que modifiqué la primer partícula
        E_final = energia(rho, m, x_modif, v_modif, p_modif, N)
        Delta_E = E_final - E_inicial
        if Delta_E > 0:
            # defino la probabilidad de aceptación
            p_aceptacion = np.exp(-Delta_E*beta)
            rand = np.random.rand()    
            if p_aceptacion < rand:
                # si no acepto el cambio, entonces vuelvo atrás
                # la modificación al vector y además redefino
                # la energía final como la inicial (porque no acepté)
                x_modif[i,:] = x[i, :]
                v_modif[i,:] = v[i, :]
                p_modif[i] = p[i]
                E_final = E_inicial

    return E_final, x_modif, v_modif, p_modif


#%%
def main():
    # parámetros para la simu
    N = 100
    m = 1
    Rho = np.arange(0.002,0.082,0.01)
    iteraciones = 4000
    # inicializo un contador k    
    k = 0
    # objetos de datos
    matriz_energia = np.zeros([len(Rho),iteraciones])
    
    for rho in Rho:
        # distribuyo x,v,p para cada rho
        x,v = distribuir_xv(rho, N)
        p = modulo_p(m, v)
        # printeo por donde voy
        print("rho = %.3f" %rho)
        t1 = time.time()    
        for i in range(iteraciones):
            ener, x, v, p = montecarlojit(rho, m, x, v, p, N)
            matriz_energia[k, i] = ener
        t2 = time.time() - t1
        print(t2)
        # le sumo 1
        k += 1
    return matriz_energia, Rho
#%%
matriz_energia, Rho = main()
energias = np.zeros(len(Rho))
longitud = len(matriz_energia[0,:])
inicio = longitud//2
for i in range(len(Rho)):
    energias[i] = np.sum(matriz_energia[i, inicio:])/len(matriz_energia[i, inicio:])
#%%
plt.figure(figsize = (8,6))
#plt.title("Simulación en Python")
plt.plot(Rho/0.037, energias,
         '.',
         markersize = 15,
         label = "Simulada N = 100")
R = np.linspace(0,2,100)
plt.plot(R, Energia_exacta(R, .0238), label = "Exacta m = 0.024")
plt.ylabel("$E/N$",fontsize = 15)
plt.xlabel(r"$\rho/\rho_{0} = R$", fontsize = 15)
plt.legend(fontsize = 15)
plt.grid()
#plt.savefig("SimulacionE_ypaper.pdf")
#%%

import numpy as np
import matplotlib.pyplot as plt


rho = 0.082
N = 64
m = 1
L = (N/rho)**(1/3)
x_vec, v_vec = distribuir_xv(rho, N)
p = modulo_p(m, v_vec)
fig = plt.figure(figsize = (10,10))

for i in range(200):
    e, x_vec, v_vec, p = montecarlojit(rho, m, x_vec, v_vec, p, N)
    cx = np.array(x_vec[0,0])
    cy = np.array(x_vec[0,1])
    cz = np.array(x_vec[0,2])
    
    x = np.array(x_vec[1:,0])
    y = np.array(x_vec[1:,1])
    z = np.array(x_vec[1:,2])
    ax = plt.axes(projection = "3d")
    ax.scatter3D(x,y,z, zdir = "z", s = 60)
    ax.scatter3D(cx,cy,cz, zdir = "z", s = 150)
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_zlim(0, L)
    ax.text(0, 0, 0, "i = %i, E = %f" %(i, e), zdir = "x")
    plt.pause(.001)
#%%
"""
En esta celda lo que quiero probar es si con C y con python, al distribuir
random en cada iteración, obtengo en promedio los mismos valores. Que me den 
los mismos valores medios podría significar:
    1) coincidencia
    2) que la generación de numeros random son semejantes a la vez que los
    términos de la energía dan los mismos resultados
"""
# acá defino los parámetros iniciales
m = 1
rho = 0.05
N = 1000

# acá defino
cin_lista = []
th_lista = []
pauli_lista = []
t1 = time.time()
for i in range(1000):
    x, v = distribuir_xv(rho, N)
    p = modulo_p(m, v)
    cin = cinetica(m, p, N)
    cin_lista.append(cin)
    th = theta(rho, p, N)
    th_lista.append(th)
    Vp = Vpauli(rho, x, v, N)
    pauli_lista.append(Vp)
    # print("cin = %f, th = %f, Vp = %f" %(cin, th, Vp))
t2 = time.time() - t1
cin = np.array(cin_lista)
th = np.array(th_lista)
pauli = np.array(pauli_lista)

##
# Ahora importo lo que viene de usar C
##
Dir = r"C:\Users\ignag\Desktop\Igna\Facultad\Fisica computacional\FINAL\Codigo\src\Tests\test.txt"
Simu_C = np.loadtxt(Dir)

C_cin = np.array(Simu_C[:,0])
C_th = np.array(Simu_C[:,1])
C_pauli = np.array(Simu_C[:,2])
print(t2)
print("<cin..C> = %f <th..C> = %f <pauli..C> = %f"
      %(C_cin.mean(),C_th.mean(),C_pauli.mean()))
print("<cin.py> = %f <th.py> = %f <pauli.py> = %f" 
      %(cin.mean(),th.mean(),pauli.mean()))
print("std..C = %f std..C = %f std..C = %f" 
      %(C_cin.std(),C_th.std(),C_pauli.std()))
print("std.py = %f std.py = %f std.py = %f" 
      %(cin.std(),th.std(),pauli.std()))

del x, v, rho, i, m, N, p, Dir, Simu_C, t1, t2
# siendo las 22:53 del 29/10/2020 veo que tanto C como python hacen lo mismo
# en cuanto al cálculo de los términos de la energía y la distribución de
# partículas. (probe con rho = 0.002, 0.05 N = 100,
# rho 0.05 y N 1000)
# el archivo para este test lo borré.
#%%
"""
Ya sé que C y Python hacen lo mismo en cuanto a la distribución de coordenadas
y calculo de los términos de la energía. Ahora quiero ver como se manejan cada
una cuando corro el montecarlo.
Asi que lo que quiero hacer es algo que me guarde en cada paso de montecarlo
los valores individuales de eneria cinetica, theta, pauli y total, a ver como
van difiriendo de las de python.

Esta se puede correr sin problema que no tarda mucho (70 segundos para 4000 
iteraciones)
"""
# parámetros de la simulación
N = 100
m = 1
rho = 0.002
iteraciones = 4000
# condiciones iniciales de la simulación
x, v = distribuir_xv(rho, N)
p = modulo_p(m, v)
# listas para llenar con los datos de montecarlo
vec_cin = []
vec_th = []
vec_pauli = []
vec_E = []
t1 = time.time()
for i in range(iteraciones):
    E_f, x, v, p = montecarlojit(rho, m, x, v, p, N)
    vec_cin.append(cinetica(m, p, N))
    vec_th.append(theta(rho, p ,N))
    vec_pauli.append(Vpauli(rho, x, v, N))
    vec_E.append(E_f)
    print(i) if i%100 == 0 else None
t2 = time.time() - t1
print(t2)

plt.plot(vec_E, label = "$E_{t}$")
plt.plot(vec_cin, label = "$E_{cin}$")
plt.plot(vec_th, label = r"$E_{\theta}$")
plt.plot(vec_pauli, label = "$E_{Pauli}$")
plt.axis([1,4050,-.2,2.2])
plt.xscale("log")
plt.legend()

del N, m, rho, iteraciones, x, v, p, t1, t2, i, E_f

Dir = r"C:\Users\ignag\Desktop\Igna\Facultad\Fisica computacional"\
       "\FINAL\Codigo\src\Tests\\test_montecarlo.txt"
Datos = np.loadtxt(Dir)
Cvec_E = np.array(Datos[:,0])
Cvec_cin = np.array(Datos[:,1])
Cvec_th = np.array(Datos[:,2])
Cvec_pauli = np.array(Datos[:,3])
plt.figure()
plt.plot(Cvec_E, label = "C$E_{t}$")
plt.plot(Cvec_cin, label = "C$E_{cin}$")
plt.plot(Cvec_th, label = r"C$E_{\theta}$")
plt.plot(Cvec_pauli, label = "C$E_{Pauli}$")
plt.axis([1,4050,-.2,2.2])
plt.xscale("log")
plt.legend()
del Dir, Datos
# Veo que son diferentes todas las energias, más que nada la cinética (y por
# ende la total), pero para V_pauli veo que son bastante parecidas. Sin embargo
# las energías más fáciles de calcular (cinética) son las que más distinto dan
# no sé por qué
#%%
"""
Acá voy a probar qué tal la distribución de momentos tanto de C como de python
esta parte es very difficul porque hay que hacer el histograma y no me estaría
saliendo el guacho
"""
N = 100
m = 1
rho = 0.042
x, v = distribuir_xv(rho, N)
p = modulo_p(m, v)
p_f = p_Fermi(rho/0.037)
Dir = r"C:\Users\ignag\Desktop\Igna\Facultad\Fisica computacional"\
       "\FINAL\Codigo\src\Tests\\test_distribucion_p1.txt"
Datos = np.loadtxt(Dir)
CE = Datos[:,0]
Ccin = Datos[:,1]
Cq = ((Ccin*2)**(1/2))/p_f
#Cq = Datos[:,2]

# plt.plot(CE)
# plt.plot(Ccin)
# plt.plot(Cq)

bines = np.linspace(0, 1.2, 50)

hist, bin_edges = np.histogram(Cq,
                               bins = bines,
                               density = True)
plt.figure(figsize = (8,6))
plt.plot(bin_edges[1:],hist, '.-', label = "Simulacion N = 100")
#distribución exacta
q = np.linspace(0, 1.2, 200)
plt.plot(q, Dist_momentos_exacta(q), label = "Exacta")
plt.xlabel("$q = p/p_{F}$",fontsize = 15)
plt.ylabel(r"$f(q,\tau)$",fontsize = 15)
plt.grid()
plt.legend(fontsize = 15)
del N, m, rho, x, v, p, Dir, Datos 
plt.savefig("Distribucion_q.pdf")
