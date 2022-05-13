# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 17:26:11 2021

@author: rehpe
"""
from param import *

def backward(t,j,tempo_aux,stage,res,stages,cen,g,runtime_limite,sub_resolve):
    from tree import tree 
    import time
    from dateutil.relativedelta import relativedelta
    tree1 = tree(int(g.loc[j][0]),abert);
    PI = {};PIy = {};pi_matrix={};b_matrix = {};flag_time = np.zeros(abert);custo = np.zeros(abert);
    PIs={};PIys={};vsults = {};custos = {};ysults={}
    for s in cen:
        PIys[s]=np.array([])
        prob = np.zeros(abert)
    #     # custo = np.array([]);
        for i in range(abert):
            for h in h_agreg:
                if bacia==0:h_aux = hidros.COMPATIBILIDADE_SPT.loc[h]
                else: h_aux = h
                for q in range(ordem):
                    stage[j][f'yp{h}_{q}'].LB = cen[s][f'ysult{h}_{q}'][j-2]
                    stage[j][f'yp{h}_{q}'].UB = cen[s][f'ysult{h}_{q}'][j-2]
                stage[j][f'ctr{h}_1'].RHS = res[(res.Iteradores==t)][f'{i+1}'].loc[h_aux]

            for h in reservatorios:
                    stage[j][f'ctr_{h}'].RHS = cen[s][f'vsult{h}'][j-2]
            stage[j][f'm{j}'].update()
            stage[j][f'm{j}'].Params.timeLimit=80
            # stage[j][f'm{j}'].Params.timeLimit=runtime_limite
            # stage[j][f'm{j}'].Params.timeLimit=0.0009
            stage[j][f'm{j}'].optimize()
            
            # stage[j][f'm{j}'].write(f'SDDP_back{j}{i}{sev}.lp')
            # stage[j][f'm{j}'].write(f'SDDP_back{j}{i}{sev}.sol')
            stage[j][f'm{j}'].write(f'SDDP_back_{j}{i}_v4.lp')
            stage[j][f'm{j}'].write(f'SDDP_back_{j}{i}_v4.sol')
            tempo_aux[j][s][1][i] = np.zeros(3)
            # tempo_aux[j][s][1][i] = stage[j][f'm{j}'].runtime
            tempo_aux[j][s][1][i][0] = stage[j][f'm{j}'].runtime
            tempo_aux[j][s][1][i][1] = stage[j][f'm{j}'].NumConstrs
            tempo_aux[j][s][1][i][2] = stage[j][f'm{j}'].NumVars - stage[j][f'm{j}'].NumBinVars
            status = stage[j][f'm{j}'].status
            if status == gb.GRB.INFEASIBLE:
                stage[j][f'm{j}'].computeIIS()
                stage[j][f'm{j}'].write("model.ilp")
                print('backward infeasible')
            elif status == 9:
                # print('passou do tempo',stage[j][f'm{j}'].runtime)
                flag_time[i] = 1
                # b_matrix[i] = stage[j][f'm{j}'].getAttr('RHS', stage[j][f'm{j}'].getConstrs())
                # b_matrix[i] = np.array([stage[j][f'ctr{b}_1'].RHS for b in hidros.BACIA.unique()])
                b_matrix[i] = np.array([stage[j][f'ctr{h}_1'].RHS for h in h_agreg])
            elif status != 2:
                print("não é ótimo")
            else: 
                PI[i] = np.array([stage[j][f'ctr_{h}'].pi for h in reservatorios])
                custo[i] = np.array(stage[j][f'm{j}'].objVal)
                PIy[i] = np.array([stage[j][f'yp{h}_{q}'].RC for h in h_agreg for q in range(ordem)])
                pi_matrix[i] =  np.array([stage[j][f'ctr{h}_1'].pi for h in h_agreg]) 
                b_matrix[i] = np.array([stage[j][f'ctr{h}_1'].RHS for h in h_agreg])
                sub_resolve[j][i] = 1
                # stage[j][f'm{j}'].write(f'SDDP_back{j}sinc.lp')
                # stage[j][f'm{j}'].write(f'SDDP_back{j}sinc.sol')
        
        idx = np.where(flag_time==1)[0]
        flag0 = np.where(flag_time==0)[0]
        
        if flag0.size==0: 
            flag_time = np.zeros(abert)
            sorteados = rd.sample(range(abert),int(0.5*abert))
            for i in sorteados:
                for h in h_agreg:
                    if bacia==0:h_aux = hidros.COMPATIBILIDADE_SPT.loc[h]
                    else: h_aux = h
                    for q in range(ordem):
                        stage[j][f'yp{h}_{q}'].LB = cen[s][f'ysult{h}_{q}'][j-2]
                        stage[j][f'yp{h}_{q}'].UB = cen[s][f'ysult{h}_{q}'][j-2]
                    stage[j][f'ctr{h}_1'].RHS = res[(res.Iteradores==t)][f'{i+1}'].loc[h_aux]
                for h in reservatorios:
                    stage[j][f'ctr_{h}'].RHS = cen[s][f'vsult{h}'][j-2]

                stage[j][f'm{j}'].update()
                stage[j][f'm{j}'].Params.timeLimit=80
                # stage[j][f'm{j}'].Params.timeLimit=0.0009
                stage[j][f'm{j}'].optimize()
                sub_resolve[j][i] = 1
                # stage[j][f'm{j}'].write(f'SDDP_back{j}{i}v4.lp')
                # stage[j][f'm{j}'].write(f'SDDP_back{j}{i}v4.sol')
                tempo_aux[j][s][1][i] = np.zeros(3)
                tempo_aux[j][s][1][i][0] = stage[j][f'm{j}'].runtime
                tempo_aux[j][s][1][i][1] = stage[j][f'm{j}'].NumConstrs
                tempo_aux[j][s][1][i][2] = stage[j][f'm{j}'].NumVars - stage[j][f'm{j}'].NumBinVars
                # tempo_aux[j][s][1][i] = stage[j][f'm{j}'].runtime
                status = stage[j][f'm{j}'].status
                if status == gb.GRB.INFEASIBLE:
                    stage[j][f'm{j}'].computeIIS()
                    stage[j][f'm{j}'].write("model.ilp")
                    print('backward infeasible')
                elif status != 2:
                    print("não é ótimo")
                else: 
                    PI[i] = np.array([stage[j][f'ctr_{h}'].pi for h in reservatorios])
                    PIy[i] = np.array([stage[j][f'yp{h}_{q}'].RC for h in h_agreg for q in range(ordem)])
                    custo[i] = np.array(stage[j][f'm{j}'].objVal)
                    pi_matrix[i] =  np.array([stage[j][f'ctr{h}_1'].pi for h in h_agreg]) 
                    b_matrix[i] = np.array([stage[j][f'ctr{h}_1'].RHS for h in h_agreg])
            for i in [dix for dix in list(np.arange(abert)) if dix not in sorteados]:
                flag_time[i] = 1
                
        #tempu = time.time()
        idx = np.where(flag_time==1)[0]
        flag0 = np.where(flag_time==0)[0]      
        menores = np.argsort(custo)[:mt.floor((risk)*abert)]
        maiores = np.argsort(custo)[mt.floor((risk)*abert):]      
        for i in idx:
            custo_linha = [custo[ix]-np.inner(b_matrix[ix],pi_matrix[ix])+np.inner(b_matrix[i],pi_matrix[ix]) for ix in flag0]
            PI[i] = PI[flag0[np.argmax(custo_linha)]]
            PIy[i] = PIy[flag0[np.argmax(custo_linha)]]
            custo[i] = np.max(custo_linha)
        for k in maiores: prob[k]=prob2
        for k in menores: prob[k]=prob1
        
        if mc==1:
            PIs[s]={}
            for i in range(abert): 
                PIs[s][i] = np.multiply(PI[i],prob[i])  
            vsults[s] = np.array(cen[s][f'vsult{h}'][j-2] for h in reservatorios) 
            custos[s][i] = np.multiply(custo[i],prob[i]) - np.inner(vsults[s],PIs[s][i])
        else:
            aux = np.zeros(len(reservatorios));
            aux2 = np.zeros(len(h_agreg))
            
            for i in range(abert): 
                aux += np.multiply(PI[i],prob[i]) 
                aux2 += np.multiply(PIy[i],prob[i]) 
            PIys[s] = aux2; PIs[s] = aux; 
            
            vsults[s] = np.array([cen[s][f'vsult{h}'][j-2] for h in reservatorios])
            ysults[s] = np.array([cen[s][f'ysult{h}_{q}'][j-2-q] for h in h_agreg for q in range(ordem)])
            custos[s] = np.inner(custo,prob) - np.inner(vsults[s],PIs[s])- sum(np.inner(ysults[s][q*len(h_agreg):(q+1)*len(h_agreg)],PIys[s][q*len(h_agreg):(q+1)*len(h_agreg)]) for q in range(ordem))
        #tap = time.time() - tempu 
    return PIs,PIys,vsults,ysults,custos,tempo_aux,sub_resolve#,tap