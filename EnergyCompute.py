from DelayMin import Params, Get_Edge, Get_Rlk
import gurobipy as gp
import pandas as pd
import numpy as np
import math

## Parameters
J = 10
K = [i for i in range(1,10)]
# Struc = 'Serial'
Struc = 'Parallel'
# Struc = 'General'

## Compute Energy


E = Get_Edge(Struc)

# !能耗的计算还存在一些问题。
for K_i in K:
    R = Get_Rlk(num_servers=K_i+1)
    df_A = pd.read_csv(f'A_matrix_J_{J}_K_{K_i}_{Struc}.csv', index_col=0)
    A = np.array(df_A).tolist()
    Energy_K =  [0 for _ in range(K_i+1)]
    for k_i in range(K_i+1):
        Energytemp1 = gp.quicksum(A[j][k_i]*Params.kappa*(Params.c_j)*((Params.f_k/(10**9))**2) for j in range(J))
        Energytemp2 = gp.quicksum(A[j][k_i]*gp.quicksum(E[j][i]*Params.d_j*A[i][l]*(10**(Params.Pk/10)/1000)/R[k_i][l] for l in range(K_i+1) for i in range(j+1,J)) for j in range(J))
        Energy_K[k_i] = Energytemp1.getConstant() + Energytemp2.getConstant() # in Watts
    # print(sum(Energy_K))
    with open('./EnergyData.txt','a+') as f:
        f.write(f'J={J}K={K_i}_{Struc}, {sum(Energy_K)}\n')
    f.close()