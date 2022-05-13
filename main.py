# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 12:40:12 2019

@author: renat
"""


import time
import sys
import os
from dateutil.relativedelta import relativedelta
from tree import tree 
from mpi4py import MPI

from param import *
from forward_AR import forward
from backward import backward
from modelo import modelo

rd.seed(seed)
simu=0;
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

def sendcorte(pacote):
   N=len(reservatorios)+len(h_agreg)
   sendcountes=[len(data[n])*(N*2+1) for n in range(len(data))]
   offsets = np.zeros(nprocs) 
   offsets[1:]=np.cumsum(sendcountes)[:-1]
   pacotao = np.zeros((sum(len(data[n]) for n in range(len(data))),2*N+1), dtype = 'd')
   comm.Allgatherv([pacote, MPI.DOUBLE], [pacotao, sendcountes, offsets, MPI.DOUBLE], )
   
   return pacotao
def cut_selection_func(corrected_list,cortes,vetor_v,vetor_i):
    for i in sorted(cortes.keys()-corrected_list):
        a_ma = np.array([]);
        for k in cortes.keys():
            if k < i:
                aux_vetor = cortes[i][-1] + np.inner(cortes[k][len(reservatorios):2*len(reservatorios)],cortes[i][:len(reservatorios)])\
                            + sum(np.inner(cortes[k][2*len(reservatorios)+(q+ordem)*len(h_agreg):2*len(reservatorios)+(q+ordem+1)*len(h_agreg)],cortes[i][2*len(reservatorios)+q*len(h_agreg):2*len(reservatorios)+(q+1)*len(h_agreg)]) for q in range(ordem))
                if vetor_v[k]<aux_vetor:
                    vetor_v[k]=aux_vetor
                    vetor_i[k]=i
            a_ma = np.append(a_ma,cortes[k][-1] + np.inner(cortes[i][len(reservatorios):2*len(reservatorios)],cortes[k][:len(reservatorios)])\
                   + sum(np.inner(cortes[i][2*len(reservatorios)+(q+ordem)*len(h_agreg):2*len(reservatorios)+(q+ordem+1)*len(h_agreg)],cortes[k][2*len(reservatorios)+q*len(h_agreg):2*len(reservatorios)+(q+1)*len(h_agreg)]) for q in range(ordem)))
        vetor_v[i]= np.amax(a_ma)
        vetor_i[i] = [*cortes][np.argmax(a_ma)]
    corrected_list = set(val for val in vetor_i.values())  
    return corrected_list,vetor_v,vetor_i



stages = len(g);
stage,aux_mc = modelo(mc,g,f,0,0)
nos_per=[abert**(i) for i in range(5)]
nos_per_aux=[];aux = 0
for i in nos_per:
    aux += i
    nos_per_aux.append(aux) 

start_program = time.time();resolve={}
o=0;time_stop = time.time() - start_program;time_ot = [];passo = 5;o_aux = 0;NC_PDDE=1
custo_inf = [];para = 1;tol = []; custo_ant = 0;o=0;tempo = {};custo_sup=[]
k=0;v_estado = {};y_estado = {};num_cut={};corte_coef={}
vetor_v = {};vetor_i = {};num_ant = {}
corrected_list={};itera = {}
cuts = {};tempo_selecao = {};tempo_selecao_adicao={};cortes={};cortes_armazen = {};tempo_fase = {}
for j in range(1,len(g)):
  cortes[j] = {};vetor_v[j]={};vetor_i[j]={};itera[j]=np.array([]);num_ant[j]=0;corrected_list[j]=[]
  ini_c = nos_per_aux[g.loc[j][0]-1]-abert**(g.loc[j][0]-1)+1;
  fini_c = nos_per_aux[g.loc[j][0]-1]+1
  cortes_armazen[j] = {};num_cut[j]=0;cuts[j]={};
  for a in range(ini_c,fini_c):
        cuts[j][a]=[];

# while time_stop<7200:
for o in range(5):
    if para==0 and o>1:break
    # vsult_s[o] = {}; custo_s[o] = {};ysult_s[o] = {};
    sampled_scenarios = []
    for s in range(NC_PDDE):
        sampled_scenarios.append([])
        for j in range(1,len(g)+1):
            ini = nos_per_aux[g.loc[j][0]-1]-abert**(g.loc[j][0]-1)+1
            fini = nos_per_aux[g.loc[j][0]-1]+1
            if g.loc[j][0]==1:ini=0;fini=abert
            for t in range(g.loc[j][0]):
                sampled_scenarios[s] = sampled_scenarios[s]+rd.sample(range(ini,fini),1)

    cen = {}
    for s in range(NC_PDDE):
        cen[s]={};
        cen[s]['no'] = sampled_scenarios[s]
        cen[s]['y']=[]
        aux = 0;aux1 = 1
        for j in range(len(cen[s]['no'])-1):
            ac=0;
            if int(g.loc[aux1][0])==1:cen[s]['y'].append(cen[s]['no'][j])
            else:
                ini = nos_per_aux[g.loc[aux1][0]-1]-abert**(g.loc[aux1][0]-1)+1
                fini = nos_per_aux[g.loc[aux1][0]-1]+1
                for n in range(ini,cen[s]['no'][j]):
                    ac+=1; 
                    if ac>=abert: ac=0  
                cen[s]['y'].append(ac)
            aux+=1
            if aux==g.loc[aux1][0]:aux1+=1;aux=0

    if rank == 0:
        data = np.arange(len(cen))
    
        # determine the size of each sub-task
        ave, av = divmod(data.size, nprocs)
        counts = [ave + 1 if p < av else ave for p in range(nprocs)]
    
        # determine the starting and ending indices of each sub-task
        starts = [sum(counts[:p]) for p in range(nprocs)]
        ends = [sum(counts[:p+1]) for p in range(nprocs)]
    
        # converts data into a list of arrays 
        data = [data[starts[p]:ends[p]] for p in range(nprocs)]

    else:
        data = None

    # sends data to all processes
    data = comm.bcast(data, root=0)
    
    
    # defines the scenarios that each process will solve
    if data[rank].size!=0:
        cen = dict((key,value) for key, value in cen.items() if key in data[rank])
    else: cen = {}
    
    tempo_fase[o]={};temp = time.time()
    # runs a forward phase
    cen, stage,list_st,tempo_aux,tempo_aux_aux = forward(stage,res,stages,cen,g,det,aux_mc,simu)
    tempo_fase[o][0] = time.time() - temp
    
    runtime_limite = 10*tempo_aux_aux.mean()
    if runtime_limite==0: runtime_limit = 0.001

    ## lower bound calculation

    if rank ==0:
        custo_inf.append(abs(cen[data[rank][0]]['custo'][1][0]))
        
        tol.append((abs(cen[data[rank][0]]['custo'][1][0])-custo_ant)/(abs(cen[data[rank][0]]['custo'][1][0])+0.00000001))

        custo_ant = copy.deepcopy(abs(cen[data[rank][0]]['custo'][1][0]))
    
    ## Backward
    tempo_selecao[o]={};tempo_selecao_adicao[o]={};
    for j in range(1,len(g)+1):
        tempo_selecao[o][j]={};tempo_selecao_adicao[o][j]={};
        for s in cen:
            tempo_aux[j][s][1]={} ## 1 indica backward
    t=start_month+relativedelta(months=T)
    sub_resolve = {};
    for j in range(2,stages+1):
        sub_resolve[j] = np.zeros(abert)
    tempo_fase[o][1] = 0
    for aux in range(1,len(g)):
        j = stages - aux + 1; t-=relativedelta(months=g.loc[j][0]);
        temp = time.time()
        PIs,PIys,vsults,ysults,custos,tempo_aux,sub_resolve = backward(t,j,tempo_aux,stage,res,stages,cen,g,runtime_limite,sub_resolve)
        comm.barrier()
        tempo_fase[o][1] += time.time() - temp
        pacote = np.array([np.hstack((PIs[s],vsults[s],PIys[s],ysults[s],custos[s])) for s in cen.keys()])
        pacotao = sendcorte(pacote)
        
        itera[j-1]=np.append(itera[j-1],len(pacotao)*[o])
        for nuc in pacotao:
            num_cut[j-1]+=1
            cortes[j-1][num_cut[j-1]] = nuc
            cortes_armazen[j-1][num_cut[j-1]] = nuc
        ini_c = nos_per_aux[g.loc[j-1][0]-1]-abert**(g.loc[j-1][0]-1)+1;
        fini_c = nos_per_aux[g.loc[j-1][0]-1]+1
        if cut_selection==1:
            tempo_selecao_aux = time.time()
            corrected_list[j-1],vetor_v[j-1],vetor_i[j-1] = cut_selection_func(corrected_list[j-1],cortes[j-1],vetor_v[j-1],vetor_i[j-1]) 
            print(corrected_list)
            tempo_selecao[o][j][rank] = time.time()-tempo_selecao_aux
            tempo_selecao_aux2 = time.time()
            stage[j-1][f'm{j-1}'].remove(cuts[j-1])
            cuts[j-1] = stage[j-1][f'm{j-1}'].addConstrs((stage[j-1][f'alfa_{a}']\
                        - gb.quicksum(stage[j-1][f'v{reservatorios[h]}_{a}']*cortes[j-1][i][h] for h in range(len(reservatorios)))\
                        - gb.quicksum(gb.quicksum(stage[j-1][f'yn{h}_{a}']*cortes[j-1][i][h-1+2*len(reservatorios)+q*len(h_agreg)] for h in h_agreg) for q in range(ordem))>=\
                        cortes[j-1][i][-1]) for a in range(ini_c,fini_c) for i in corrected_list[j-1])                    
            cortes[j-1] = {k: v for k, v in cortes[j-1].items() if k in corrected_list[j-1]}
            vetor_v[j-1] = {k: v for k, v in vetor_v[j-1].items() if k in corrected_list[j-1]}
            vetor_i[j-1] = {k: v for k, v in vetor_i[j-1].items() if k in corrected_list[j-1]}
            tempo_selecao_adicao[o][j][rank] = time.time()-tempo_selecao_aux2 
        else:
            cuts[j-1] = stage[j-1][f'm{j-1}'].addConstrs((stage[j-1][f'alfa_{a}']\
                        - gb.quicksum(stage[j-1][f'v{reservatorios[h]}_{a}']*cortes[j-1][i][h] for h in range(len(reservatorios)))\
                        - gb.quicksum(gb.quicksum(stage[j-1][f'yn{h}_{a}']*cortes[j-1][i][h-1+2*len(reservatorios)+q*len(h_agreg)] for h in h_agreg) for q in range(ordem))>=\
                        cortes[j-1][i][-1]) for a in range(ini_c,fini_c) for i in range(num_ant[j-1]+1,len(cortes[j-1])+1))                    
           
        num_ant[j-1] = num_cut[j-1]

        stage[j-1][f'm{j-1}'].update();resolve[o] = sub_resolve
        
        stage[j-1][f'm{j-1}'].write(f'modelo{j-1}_3.lp')
    # tempo_fase[o][1] = time.time() - temp
    if rank >= 1:
        comm.send(tempo_aux, dest=0, tag=rank)
    else:
        for i in range(1,nprocs):
            tempo_aux1 = comm.recv(source=i, tag=i)
            for j in range(1,stages+1):
                tempo_aux[j].update(tempo_aux1[j])
        tempo[o] = tempo_aux

    time_ot.append((time.time() - start_program)/3600)    
    time_stop = time.time() - start_program;o+=1;

end_program = time.time() - start_program;
df = {'o':[],'j':[],'s':[],'tempo':[]}
for k in range(o):
    for j in tempo_selecao[k]:
        if j==1:continue
        for s in tempo_selecao[k][j]:
            df['o'].append(k)
            df['j'].append(j)
            df['s'].append(s)
            df['tempo'].append(tempo_selecao[k][j][s])
df = pd.DataFrame(df)
df.to_csv(os.path.join(path1,f'tempo_selecao{rank}.csv'),sep = ';')
df = {'o':[],'j':[],'s':[],'tempo':[]}
for k in range(o):
    for j in tempo_selecao[k]:
        if j==1:continue
        for s in tempo_selecao[k][j]:
            df['o'].append(k)
            df['j'].append(j)
            df['s'].append(s)
            df['tempo'].append(tempo_selecao_adicao[k][j][s])
df = pd.DataFrame(df)
df.to_csv(os.path.join(path1,f'tempo_selecao_adicao{rank}.csv'),sep = ';')
df = {'j':[],'num cortes':[],'cortes':[]}
for j in corrected_list:
        if j==1:continue
        df['j'].append(j)
        df['cortes'].append(corrected_list[j])
        df['num cortes'].append(len(corrected_list[j]))
df = pd.DataFrame(df)
df.to_csv(os.path.join(path1,f'selecao{rank}.csv'),sep = ';')


if rank==0:
    t1 = pd.DataFrame([end_program], columns = ['tempo total'])
    t2 = pd.DataFrame([o], columns = ['iteracoes'])
    t3 = pd.DataFrame(custo_inf, columns = ['limite_inferior'])
    t4 = pd.DataFrame(time_ot, columns = ['tempo por iteracao'])
    t5 = pd.DataFrame(tol, columns = ['tolerancia'])
    t6 = pd.DataFrame([custo_inf[-1]], columns = ['custo otimo'])
    # t7 = pd.DataFrame([custo_sup[-1]], columns = ['custo superior'])
    t = pd.concat([t1,t2,t3,t4,t5,t6],axis=1)
    # t.to_csv(os.path.join(path1,f'evolution_{stages}_{seed}_{phi}_{abert}.csv'),sep = ';')
    t.to_csv(os.path.join(path1,f'evolution_{stages}_{seed}_{abert}.csv'),sep = ';')


    df = {'o':[],'tempo':[],'j':[],'custo':[]}
    for h in reservatorios:
          df[f'PI{h}']=[];df[f'vsult{h}']=[]
    for h in h_agreg:
        for q in range(ordem):
            df[f'PIy{h}{q}']=[];df[f'ysult{h}{q}']=[]
    for j in cortes_armazen:
        for k in range(1,len(cortes_armazen[j])+1):
            df['o'].append(int(itera[j][k-1]))
            df['j'].append(j)
            df['tempo'].append(time_ot[int(itera[j][k-1])])
            df['custo'].append(cortes_armazen[j][k][-1])
            for h in range(len(reservatorios)):
                # print(j,h,h+len(reservatorios),cortes[j][k][h+len(reservatorios)])
                df[f'PI{reservatorios[h]}'].append(cortes_armazen[j][k][h])
                df[f'vsult{reservatorios[h]}'].append(cortes_armazen[j][k][h+len(reservatorios)])
            for h in h_agreg:
                for q in range(ordem):
                    df[f'PIy{h}{q}'].append(cortes_armazen[j][k][h+2*len(reservatorios)-1])
                    df[f'ysult{h}{q}'].append(cortes_armazen[j][k][h+2*len(reservatorios)+len(h_agreg)-1])
                  
    df = pd.DataFrame(df)
    cols=df.columns[4:]
    df.to_csv(os.path.join(path1,f'cortes_{stages}_{seed}_v3.csv'),sep = ';')  
    df = {'o':[],'fase':[],'tempo':[],}
    for k in range(o):
        for fase in [0,1]:
            df['o'].append(k)
            df['fase'].append(fase)
            df['tempo'].append(tempo_fase[k][fase])
    df = pd.DataFrame(df)
    df.to_csv(os.path.join(path1,f'tempo_fase_{stages}_{seed}.csv'),sep = ';')   
    df = {'o':[],'j':[],'s':[],'fase':[],'abertura':[],'tempo':[]}
    for k in range(o):
        for j in tempo[k]:
            for s in tempo[k][j]:
                for fase in [0,1]:
                    if fase==1 and j==1: continue
                    if fase==0 and j==len(g): continue                    
                    if fase == 0:
                        df['o'].append(k)
                        df['j'].append(j)
                        df['s'].append(s)
                        df['fase'].append(fase)
                        df['abertura'].append(0)
                        df['tempo'].append(tempo[k][j][s][fase])
                    else:
                        for i in range(abert):
                            df['o'].append(k)
                            df['j'].append(j)
                            df['s'].append(s)
                            df['fase'].append(fase)
                            df['abertura'].append(i)
                            df['tempo'].append(tempo[k][j][s][fase][i])
    df = pd.DataFrame(df)
    df.to_csv(os.path.join(path1,f'tempo_{stages}_{seed}.csv'),sep = ';')
    df = {'o':[],'j':[],'s':[],'fase':[],'abertura':[],'tempo':[],'restricoes':[],'variaveis':[]}
    for k in range(o):
        for j in tempo[k]:
            for s in tempo[k][j]:
                for fase in [0,1]:
                    if fase==1 and j==1: continue
                    if fase==0 and j==len(g): continue                    
                    if fase == 0:
                        df['o'].append(k)
                        df['j'].append(j)
                        df['s'].append(s)
                        df['fase'].append(fase)
                        df['abertura'].append(0)
                        df['tempo'].append(tempo[k][j][s][fase][0])
                        df['restricoes'].append(tempo[k][j][s][fase][1])
                        df['variaveis'].append(tempo[k][j][s][fase][2])
                    else:
                        for i in range(abert):
                            df['o'].append(k)
                            df['j'].append(j)
                            df['s'].append(s)
                            df['fase'].append(fase)
                            df['abertura'].append(i)
                            df['tempo'].append(tempo[k][j][s][fase][i][0])
                            df['restricoes'].append(tempo[k][j][s][fase][i][1])
                            df['variaveis'].append(tempo[k][j][s][fase][i][2])
    df = pd.DataFrame(df)
    df.to_csv(os.path.join(path1,f'tempo_{stages}_{seed}.csv'),sep = ';')
    
    df = {'o':[],'j':[],'abertura':[],'resolve':[]}
    for k in range(o):
        for j in resolve[k]:
            # for i in range(abert):
            df['o'] = df['o'] + abert*[k]
            df['j'] = df['j'] + abert*[j]
            df['abertura'] = df['abertura'] + list(np.arange(abert))
            df['resolve'] = df['resolve'] + list(resolve[k][j])
    df = pd.DataFrame(df)
    df.to_csv(os.path.join(path1,f'subproblemas_resolvidos_{stages}_{seed}.csv'),sep = ';')
