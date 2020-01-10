# Importação de bibliotecas
from matplotlib import pyplot as plt
import math
import numpy as np

# Tamanho da grade de simulação
Tamanho_Grid = 500

# Escolha da variação espacial
Delta_Espaço = 2

# Cálculo da velocidade da luz no vácuo
mi=4*math.pi*10**-7
eps0=1*(10**-9/(36*math.pi))
v0=1/(math.sqrt(mi*eps0))

# Definição da variação temporal a partir da condição de Courant
# S = sqrt(N), onde N = número de dimensões da simulação
S = 1
Delta_Tempo = Delta_Espaço/(v0*S)

# Definição do número de pontos para a simulação
Numero_de_Pontos = math.floor(Tamanho_Grid/Delta_Espaço)

# Características da gaussiana: Posição do pico (T1), largura (L) e amplitude (A)
T1=int(Numero_de_Pontos/10)
L=T1/2
A=1

# Define o tempo total de simulaçao
T2=2*Numero_de_Pontos

# Cria-se matrizes de zeros para o campo elétrico e magnético
Ex = np.zeros((T2 + 2, 2*Numero_de_Pontos +2))
Hy = np.zeros((T2 + 2, 2*Numero_de_Pontos +2))

# Loop para geração e propagação da gaussiana
for tempo in range(1, T2, 2):
    Ex[tempo][1] = A* (math.exp(-(1 / (L ** 2))* ((tempo - T1)** 2)))
    for espaco in range(1, Numero_de_Pontos, 2):
        Hy[tempo+1][espaco+1] = -((1 / mi) * (1 / v0*S)) * (Ex[tempo][espaco+2] - Ex[tempo][espaco]) + Hy[tempo-1][espaco+1]
        Ex[tempo+2][espaco] = -((1 / eps0) * (1 / v0*S)) * (Hy[tempo + 1][espaco + 1] - Hy[tempo + 1][espaco - 1]) + Ex[tempo][espaco]

# Plotagem da onda
plt.ion()
fig = plt.figure()
se = np.round(range(1,Numero_de_Pontos,2))
espaço = []
for i in range(int(Delta_Espaço), Tamanho_Grid + 1, int(2*Delta_Espaço)):
    espaço.append(i)
i=0
for te in range(1, T2, 2):
     plt.clf()
     plt.plot(espaço, Ex[te][se])
     axes = plt.gca()
     axes.set_xlim([0, Tamanho_Grid])
     axes.set_ylim([-A*1.5, A*1.5])
     plt.xlabel("Distância (m)")
     plt.ylabel("Campo Elétrico (V/m)")
     fig.suptitle("Tempo percorrido: {}s".format(tempo))
     plt.savefig("Simulação de numero {}.png".format(i))
     fig.canvas.draw()
     fig.canvas.flush_events()
     tempo = te*Delta_Tempo
     i = i+1
