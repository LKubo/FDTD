# Importação de bibliotecas
from matplotlib import pyplot as plt
import math
import numpy as np

# Tamanho da grade de simulação
Tamanho_Grid = 2000

# Escolha do passo espacial
Delta_Espaço = 2

# Cálculo da velocidade da luz no vacuo
mi = 4*math.pi*10**-7
eps0 = 1*(10**-9/(36*math.pi))
v0 = 1/(math.sqrt(mi*eps0))

# Definição da variação temporal, a partir da condição de Courant
# S=sqrt(N), onde N=1
S = 1
Delta_Tempo = Delta_Espaço / (S*v0)

# Escolha das permissividades elétricas relativas dos meios (eps1 = N1*eps0)
N1 = 5
N2 = 1
n1 = math.sqrt(N1)
n2 = math.sqrt(N2)
# Cálculo dos novos valores de epslon dos meios
eps1 = N1*eps0
eps2 = N2*eps0

v1 = 1/(math.sqrt(mi*eps2))

# Número de pontos a serem simulados
Numero_Pontos = math.floor(Tamanho_Grid / Delta_Espaço)

# Características da gaussiana: Posição do pico (T1), largura (L) e amplitude (A)
T1 = Numero_Pontos / 10
L = T1/2
A = 400

# Tempo total de simulação T3
T2 = 3 * Numero_Pontos

# Criação de vetores para os campos elétricos e magnéticos
Ex = np.zeros((T2 + 2, 2 * Numero_Pontos))
Hy = np.zeros((T2 + 2, 2 * Numero_Pontos))

Max1 = 0
Max2 = 0
Min1 = 0
Min2 = 0
# Loop para geração da gaussiana
for tempo in range(1, T2, 2):
    Ex[tempo][1] = A * (math.exp(-(1 / (L ** 2)) * ((tempo - T1) ** 2)))
    for espaco in range(1, int(Numero_Pontos / 2), 2):
        Hy[tempo+1][espaco+1] = -((1 / mi) * (Delta_Tempo / Delta_Espaço)) * (Ex[tempo][espaco + 2] - Ex[tempo][espaco]) + Hy[tempo - 1][espaco + 1]
        Ex[tempo+2][espaco] = -((1 / eps1) * (Delta_Tempo / Delta_Espaço)) * (Hy[tempo + 1][espaco + 1] - Hy[tempo + 1][espaco - 1]) + Ex[tempo][espaco]
        if Ex[tempo][espaco]<390:
            if Ex[tempo][espaco]>Max1:
                Max1 = Ex[tempo][espaco]
            if Ex[tempo][espaco]<Min1:
                Min1 = Ex[tempo][espaco]
    for espaco in range(int(Numero_Pontos / 2) + 1, Numero_Pontos - 2, 2):
        Hy[tempo+1][espaco+1] = -((1 / mi) * (Delta_Tempo / Delta_Espaço)) * (Ex[tempo][espaco + 2] - Ex[tempo][espaco]) + Hy[tempo - 1][espaco + 1]
        Ex[tempo+2][espaco] = -((1 / eps2) * (Delta_Tempo / Delta_Espaço)) * (Hy[tempo + 1][espaco + 1] - Hy[tempo + 1][espaco - 1]) + Ex[tempo][espaco]
        if Ex[tempo][espaco]>Max2:
            Max2 = Ex[tempo][espaco]
        if Ex[tempo][espaco]<Min2:
            Min2 = Ex[tempo][espaco]

Coef = Max2/A
print("Amplitude Máxima da onda incidente: {}".format(A))
print("Amplitude Máxima da onda refratada: {}".format(Max2))
print("Amplitude Máxima da onda refletida: {}".format(Max2 - A))
print("Relação de amplitudes (Amplitude refratada/Amplitude incidente): {}".format(Coef))
print("Coeficiente de reflexão teórico: {}".format((n2-n1)/(n2+n1)))
# Plotagem da onda
plt.ion()
fig = plt.figure()
se = np.round(range(1, Numero_Pontos, 2))
espaço = []
for i in range(2 * Delta_Espaço, Tamanho_Grid + 1, 2 * Delta_Espaço):
    espaço.append(i)
i = 0
for te in range(1, T2, 8):
     plt.clf()
     plt.plot(espaço, Ex[te][se])
     plt.axvline(x=1000, c="green", linewidth=1)
     axes = plt.gca()
     axes.set_xlim([0, Tamanho_Grid])
     axes.set_ylim([-A*1.5, A*1.5])
     plt.xlabel("Distância (m)")
     plt.ylabel("Campo Elétrico (V/m)")
     fig.suptitle("Tempo percorrido: {}s".format(tempo))
     # plt.savefig("Simulação de numero {}.png".format(i))
     fig.canvas.draw()
     fig.canvas.flush_events()
     i = i+1
     tempo = te * Delta_Tempo
