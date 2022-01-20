import gurobipy as gp
from gurobipy import GRB
import math

## Basic Parameters
class Params:
    BW            = 50*(10**6)              # transmit bandwidth
    BW_db         = 10*math.log10(BW)       # BW in dB
    P1            = 10                      # transmit power P_1 = 10 dbm
    Pk            = 32                      # transmit power 2<=k<=K+1, 32dbm
    R_wired       = 100*(10**6)             # wired transmit rate, 100mbps
    nk            = 1                       # available VMs for server k
    kappa         = 10**(-11)               # constant kappa
    d_j           = 10*(10**(6))            # result data size
    c_j           = 100*(10**(6))           # required workload
    Pathloss      = 57                      # Path loss 57 dB
    sigma2        = -174                    # noise power density -174dBm/Hz
    # Parameter Control
    J             = 10                       # num of tasks
    K             = 9                       # num of servers
    f_k           = 3.2*(10**9)             # CPU frequency
    slots         = 500                     # time slots
    DAG_structure = 'General'               # DAG Structure, 'Serial', 'Parallel' or 'General' 
    # DAG_structure = 'Parallel'              # DAG Structure, 'Serial', 'Parallel' or 'General' 
    # DAG_structure = 'Serial'                # DAG Structure, 'Serial', 'Parallel' or 'General' 
    time_window   = 1*10**(-3)              # time window 1ms

## Model
model = gp.Model("OptimalOffloading")


## Variables Initialize

# Compute e_{i,j}
def Get_Edge(Structure_select):
    e  = [[0 for _ in range(Params.J)] for _ in range(Params.J)]
    if Structure_select == 'Serial':
        for i1 in range(Params.J):
            for i2 in range(Params.J):
                if i2-i1 == 1:
                    e[i1][i2] = 1
    
    elif Structure_select == 'Parallel':
        for i1 in range(Params.J):
            for i2 in range(Params.J):
                if (i1 == 0) and (i2 > i1):
                    if i2 == Params.J-1:
                        e[i1][i2]=0
                    else:
                        e[i1][i2] = 1
                elif (i2 == Params.J - 1) and (i2 > i1):
                    e[i1][i2] = 1

    elif Structure_select == 'General':
        for i1 in range(Params.J):
            for i2 in range(Params.J):
                if i2 > i1:
                    if i1 == 0 and i2 <= 2:
                        e[i1][i2]=1
                    elif i1 == 1 and 3<=i2<=4:
                        e[i1][i2]=1
                    elif i1 == 2 and 5<=i2<=6:
                        e[i1][i2]=1
                    elif i2 == 7 and 3<=i1<=4:
                        e[i1][i2]=1
                    elif i2 == 8 and 5<=i1<=6:
                        e[i1][i2]=1
                    elif i2 == 9 and 7<=i1<=8:
                        e[i1][i2]=1
    return e

# Compute R_{lk}
def Get_Rlk():
    R = [[0 for _ in range(Params.K + 1)]for _ in range(Params.K + 1)]
    for i1 in range(Params.K + 1):
        for i2 in range(Params.K + 1):
            if i1 == i2:
                R[i1][i2] = float('inf')
            else:
                if i1 == 0:
                    R[i1][i2] = Params.BW*math.log2(1+10**((Params.P1-Params.Pathloss-Params.BW_db-Params.sigma2)/10-3))
                elif i2 == 0:
                    R[i1][i2] = Params.BW*math.log2(1+10**((Params.Pk-Params.Pathloss-Params.BW_db-Params.sigma2)/10-3))
                else:
                    R[i1][i2] = Params.R_wired
    return R

E = Get_Edge(Params.DAG_structure)

print('--------------------------Edge Matrix---------------------------')
for i in range(Params.J):
    print(f'{E[i]}\n')

Rlk = Get_Rlk()

# Define variables
x  = model.addVars(Params.J, Params.K+1, Params.slots, vtype = GRB.BINARY)
a  = model.addVars(Params.J, Params.K+1, vtype = GRB.BINARY)
tj = model.addVars(Params.J, lb=0, ub=Params.slots, vtype = GRB.INTEGER)
z  = model.addVars(Params.J, lb=0, vtype = GRB.INTEGER) #T_j
model.update()

# time slots
t = [ _ for _ in range(Params.slots)]

# a & tj
model.addConstrs(a[j,k]==gp.quicksum(x[j,k,_] for _ in range(Params.slots)) for j in range(Params.J) for k in range(Params.K+1))
model.addConstrs(tj[j]==gp.quicksum(t[i]*x[j,k,i] for k in range(Params.K+1) for i in range(Params.slots)) for j in range(Params.J))

## Add Constraints
# Add constraint T_j = max_{i<j}...
for j1 in range(Params.J):
    for i1 in range(j1):
        model.addConstr(z[j1]>=E[i1][j1]*tj[i1]+E[i1][j1]*gp.quicksum([E[i1][j1]*a[j1,k]*a[i1,l]*math.ceil((Params.d_j/Rlk[l][k])/Params.time_window) for l in range(Params.K+1) for k in range(Params.K+1)]))

# Compute t_j^c
t_jc={}
for i1 in range(Params.J):
    t_jc[i1] = gp.quicksum(a[i1,k]*math.ceil((Params.c_j/Params.f_k)/Params.time_window) for k in range(Params.K+1))

# Add constraint T_j+t_j^c<=t_j
model.addConstrs(z[j]+t_jc[j] <= tj[j] for j in range(Params.J))

# Add constraint sum_k a_{jk} = 1
for j in range(Params.J):
    model.addConstr(gp.quicksum(a[j,k] for k in range(Params.K+1)) == 1)

# a_00=a_J1=1
model.addConstr(a[0, 0]==1)
model.addConstr(a[Params.J-1, 0]==1)

# Add constraint ...<n_k
Compute_slots = math.ceil((Params.c_j/Params.f_k)/Params.time_window)

for k in range(Params.K+1):
    for ts in range(Params.slots-Compute_slots+1):
        model.addConstr(gp.quicksum(x[j,k,t] for t in range(ts, ts+Compute_slots-1) for j in range(Params.J))<= Params.nk)



model.setObjective(tj[Params.J-1], GRB.MINIMIZE)

model.optimize()

print('--------------------------Var Values---------------------------')
if model.Status==GRB.OPTIMAL:
    print('# Values of a_jk(Assignment matrix) #')
    for j in range(Params.J):
        ltemp = []
        for k in range(Params.K+1):
            ltemp.append(a[j,k].X)
        print(ltemp)
    print('# Values of tj(Completion time of tasks) #')
    for _ in range(Params.J):
        print(tj[_].X)
    print('# Values of T_j(Ready time of tasks) #')
    for _ in range(Params.J):
        print(z[_].X)
    # print('# n_k in time slot check #')
    # for j in range(Params.J):
    #     print(f'----------- Task{j} -----------')
    #     for k in range(Params.K+1):
    #         print(f'----------- Task{j} on Server{k} -----------')
    #         for t in range(Params.slots):
    #             # print(x[j,k,t].X)
    #             if x[j,k,t].X == 1:
    #                 print(f'{t}')