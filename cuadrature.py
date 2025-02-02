#!/usr/bin/env python3

from scipy.special import legendre
import numpy as np
import math 

def gaussxw(N):

    a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)

    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N, dtype = float)

        p1 = np.copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        dp = (N + 1) * (p0 - x * p1) / (1 - x * x)
        dx = p1 / dp
        x -= dx
        delta = np.max(np.abs(dx))

    w = 2 * (N + 1) * (N + 1)/(N * N * (1 - x * x) * dp * dp)

    return x,w

def gaussxwab(a, b, x, w):

    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

puntos2_pesos2 = gaussxw(2) # aquí usamos "gaussxw" para optener una tupla puntos = x , pesos = w. para N=2

puntos3_pesos3 = gaussxw(3) # exactamente lo mismo a arriba pero para N=3

puntos4_pesos4 = gaussxw(4) # exactamente lo mismo a arriba pero para N=4

puntos5_pesos5 = gaussxw(5) # exactamente lo mismo a arriba pero para N=5

puntos2_pesos2_escalado = gaussxwab(1, 3, puntos2_pesos2[0], puntos2_pesos2[1]) 

puntos3_pesos3_escalado = gaussxwab(1, 3, puntos3_pesos3[0], puntos3_pesos3[1])

puntos4_pesos4_escalado = gaussxwab(1, 3, puntos4_pesos4[0], puntos4_pesos4[1])

puntos5_pesos5_escalado = gaussxwab(1, 3, puntos5_pesos5[0], puntos5_pesos5[1])

def funcion(x): # función para retornar la función que nos dan

    return (x**6) - (x**2) * np.sin(2 * x) # restornamos la función en general sin evaluar, para luego llamarla con su debido parametro

def sumatoria(puntos_pesos_escalado, funcion): # Función para generar el resultado de la integral con N=2

    resultado_N = np.sum(funcion(puntos_pesos_escalado[0])*puntos_pesos_escalado[1]) # variable a la que se le asigna el valor, usamos la función sum de numpy.

    return resultado_N # retornamos el valor

resultado_integral_N2 = sumatoria(puntos2_pesos2_escalado, funcion) # llamamos a la función ya con los parametros

print(f"Resultado para N=2 {resultado_integral_N2}") # imprimimos el resultado

resultado_integral_N3 = sumatoria(puntos3_pesos3_escalado, funcion) # llamamos a la función ya con los parametros

print(f"Resultado para N=3 {resultado_integral_N3}") # imprimimos el resultado

resultado_integral_N4 = sumatoria(puntos4_pesos4_escalado, funcion) # llamamos a la función ya con los parametros

print(f"Resultado para N=4 {resultado_integral_N4}") # imprimimos el resultado

resultado_integral_N5 = sumatoria(puntos5_pesos5_escalado, funcion) # llamamos a la función ya con los parametros

print(f"Resultado para N=5 {resultado_integral_N5}") # imprimimos el resultado
