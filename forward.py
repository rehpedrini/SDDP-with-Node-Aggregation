# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 14:51:22 2020

@author: rehpe
"""
from param import *
import time
def salva_var(cen,stage,t2,dur,nos_per_aux):
    from dateutil.relativedelta import relativedelta
    for per in range(dur):
        t2+=relativedelta(months=+per)
        if dur == 1: ini=1;fini=2
        else:ini = nos_per_aux[per-1]+1; fini  = nos_per_aux[per]+1
        cen['phsult'].append(sum(sum(sum(duracao_patamar[(duracao_patamar.Iteradores == t2)][f'{l}'].iloc[0]*stage[f'ph{h}_{l}_{n}'].x for l in range(1,PAT+1)) for h in hidros.index) for n in \
                                range(ini,fini))/abert**(per))
        cen['vsult'].append(sum(sum(stage[f'v{h}_{n}'].x for h in reservatorios)\
                                   for n in range(ini,fini))/abert**(per))
        cen['gtsult'].append(sum(sum(sum(stage[f'gt{i}_{l}_{n}'].x for l in range(1,PAT+1)) for i in termos.index)\
                                   for n in range(ini,fini))/abert**(per))
        cen['CMO'].append(sum(sum(sum(duracao_patamar[(duracao_patamar.Iteradores == t2)][f'{l}'].iloc[0]*stage[f'CMO_{sub}_{l}_{n}'].pi/horas for sub in range(1,SUBS+1)) for l in range(1,PAT+1)) for n in range(ini,fini)))
        cen['dfc'].append(sum(sum(sum(stage[f'd{sub}_{l}_{n}'].x for sub in range(1,SUBS+1)) for l in range(1,PAT+1)) for n in range(ini,fini))/abert**(per))    
        for h in reservatorios:
            cen[f'vsult_{h}'].append(sum(stage[f'v{h}_{n}'].x for n in range(ini,fini))/abert**(per))
        for h in h_agreg:
            cen[f'ysult_{h}'].append(sum(stage[f'y{h}_{n}'].x for n in range(ini,fini))/abert**(per))
        for l in range(1,PAT+1):
            cen[f'phsult_{l}'].append(sum(sum(stage[f'ph{h}_{l}_{n}'].x for h in hidros.index) for n in \
                                    range(ini,fini))/abert**(per))
            cen[f'gtsult_{l}'].append(sum(sum(stage[f'gt{i}_{l}_{n}'].x for i in termos.index)\
                                       for n in range(ini,fini))/abert**(per))
            cen[f'CMO_{l}'].append(sum(sum(stage[f'CMO_{sub}_{l}_{n}'].pi/horas for sub in range(1,SUBS+1)) for n in range(ini,fini)))
            cen[f'dfc_{l}'].append(sum(sum(stage[f'd{sub}_{l}_{n}'].x for sub in range(1,SUBS+1)) for n in range(ini,fini))/abert**(per))                          
    return cen,t2
def forward(stage,res,stages,cen,g,det,aux_mc,simu):
    from tree import tree 
    from dateutil.relativedelta import relativedelta
    tempo_aux_aux = np.zeros(stages-1)
    tempo_aux = {};list_st = []
    for j in range(1,len(g)+1):
        tempo_aux[j] = {}
        for s in cen:
            tempo_aux[j][s] = {}
            # tempo_aux[j][s][0]={} ## 0 indica forward
            tempo_aux[j][s][0] = np.zeros(3)
    nos_per=[abert**(i) for i in range(5)]
    nos_per_aux=[];aux = 0
    for i in nos_per:
        aux += i
        nos_per_aux.append(aux) 
    if det ==0:
        j=1;tree1 = tree(int(g.loc[j][0]),abert);#t=start_month
        ini = nos_per_aux[g.loc[j][0]-1]-abert**(g.loc[j][0]-1)+1
        fini = nos_per_aux[g.loc[j][0]-1]+1
    
        # for b in hidros.BACIA.unique():
        for h in h_agreg:
            ac=1;cont=0;t_aux = 1;t=start_month
            if bacia==0:h_aux = hidros.COMPATIBILIDADE_SPT.loc[h]
            else: h_aux = h
            for q in range(ordem):
                stage[j][f'yp{h}_{q}'].UB = y_hist.loc[h_aux, y_hist.columns == t-relativedelta(months=+q+1)][0]
                stage[j][f'yp{h}_{q}'].LB = y_hist.loc[h_aux, y_hist.columns == t-relativedelta(months=+q+1)][0]
            stage[j][f'ctr{h}_1'].RHS = res[(res.Iteradores==t)]['1'].loc[h_aux]

        for h in reservatorios:
            stage[j][f'ctr_{h}'].RHS = hidros.loc[h]['VOLUME_INICIAL']*(hidros.loc[h]['VOLUME_MAXIMO_OPERACIONAL']-hidros.loc[h]['VOLUME_MINIMO_OPERACIONAL'])/100
            
            
    
        stage[j][f'm{j}'].update()
        # stage[j][f'm{j}'].reset()
        stage[j][f'm{j}'].Params.timeLimit=80
        stage[j][f'm{j}'].optimize()
        
        status = stage[j][f'm{j}'].status
        if status == gb.GRB.INFEASIBLE:
            stage[j][f'm{j}'].computeIIS()
            stage[j][f'm{j}'].write("model.ilp")
            print('forward infeasible')
        elif status != 2:
            print("não é ótimo")
        list_st.append(status)
        # stage[j][f'm{j}'].write(f'SDDP_11v4.lp')
        # stage[j][f'm{j}'].write(f'SDDP_11v4.sol')
        for s in cen:
            # tempo_aux[j][s][0] = stage[j][f'm{j}'].runtime
            tempo_aux[j][s][0][0] = stage[j][f'm{j}'].runtime
            tempo_aux[j][s][0][1] = stage[j][f'm{j}'].NumConstrs
            tempo_aux[j][s][0][2] = stage[j][f'm{j}'].NumVars - stage[j][f'm{j}'].NumBinVars
            tempo_aux_aux[j-1] = stage[j][f'm{j}'].runtime
            # print(stage[j][f'm{j}'].runtime)
            # for b in hidros.BACIA.unique():
            #     if g.loc[j][0]>1:
            #         for q in range(ordem):
            #             cen[s][f'ysult{b}_{q}'] = [stage[1]['y{0}_{1}'.format(b,cen[s]['no'][q])].x]    
            #     else: 
            #         for q in range(ordem):
            #             cen[s][f'ysult{b}_{q}'] = [stage[1]['y{0}_{1}'.format(b,1)].x]
            for h in h_agreg:
                if g.loc[j][0]>1:
                    for q in range(ordem):
                        cen[s][f'ysult{h}_{q}'] = [stage[1]['yn{0}_{1}'.format(h,cen[s]['no'][q])].x] 
                else: 
                    for q in range(ordem):
                        cen[s][f'ysult{h}_{q}'] = [stage[1][f'yn{h}_1'].x]
            for h in reservatorios:
                if g.loc[j][0]>1:  cen[s][f'vsult{h}'] = [stage[1]['v{0}_{1}'.format(h,cen[s]['no'][0])].x]
                else: cen[s][f'vsult{h}'] = [stage[1][f'v{h}_1'].x]

            if simu==1:
                cen[s]['phsult'] = []; cen[s]['vsult'] = []; cen[s]['gtsult'] = [];cen[s]['dfc'] = []; cen[s]['CMO'] = []; 
                for h in reservatorios:cen[s][f'vsult_{h}'] = []
                for h in h_agreg:cen[s][f'ysult_{h}'] = []
                for l in range(1,PAT+1):cen[s][f'phsult_{l}'] = [];cen[s][f'gtsult_{l}'] = [];cen[s][f'dfc_{l}'] = [];cen[s][f'CMO_{l}'] = []
                t2=start_month
                cen[s],t2 = salva_var(cen[s],stage[j],t2,g.loc[j][0],nos_per_aux)
                
            cen[s]['corte'] = sum(sum(sum(stage[j][f'd{sub}_{l}_{n}'].x for sub in range(1,SUBS+1)) for l in range(1,PAT+1)) for n in tree1.keys())
            if mc == 1: cen[s]['custo_simu'] = copy.deepcopy(stage[j][f'm{j}'].objVal) - copy.deepcopy(sum(sum(stage[j][f'alfa_{i}_{a}'].x for a in range(abert)) for i in range(ini,fini)))/abert**(g.loc[j][0]-1)
            else: cen[s]['custo_simu'] = copy.deepcopy(stage[j][f'm{j}'].objVal) - copy.deepcopy(sum(stage[j][f'alfa_{i}'].x for i in range(ini,fini)))/abert**(g.loc[j][0]-1)
            cen[s]['custo']={}; cen[s]['custo'][1] = [copy.deepcopy(stage[j][f'm{j}'].objVal)]
            cen[s]['alfa']={}; 
            if mc ==1:cen[s]['alfa'][1] = [copy.deepcopy(sum(sum(stage[j][f'alfa_{i}_{a}'].x for a in range(abert)) for i in range(ini,fini)))/abert**(g.loc[j][0]-1)]
            else:cen[s]['alfa'][1] = [copy.deepcopy(sum(stage[j][f'alfa_{i}'].x for i in range(ini,fini)))/abert**(g.loc[j][0]-1)]
        

        for s in cen:
            inp = g.loc[1][0];t=start_month+relativedelta(months=+g.loc[1][0]);t2 = t
            for j in range(2,stages):
                tree1 = tree(int(g.loc[j][0]),abert);
                ini = nos_per_aux[g.loc[j][0]-1]-abert**(g.loc[j][0]-1)+1
                fini = nos_per_aux[g.loc[j][0]-1]+1
                # for b in hidros.BACIA.unique():
                #     for q in range(ordem):
                #          stage[j][f'ctr_y{b}_{q}'].RHS = cen[s][f'ysult{b}_{q}'][j-2]
                for h in reservatorios:
                    stage[j][f'ctr_{h}'].RHS = cen[s][f'vsult{h}'][j-2]
                for h in h_agreg:
                    for q in range(ordem):
                         stage[j][f'yp{h}_{q}'].UB = cen[s][f'ysult{h}_{q}'][j-2]
                         stage[j][f'yp{h}_{q}'].LB = cen[s][f'ysult{h}_{q}'][j-2]
                    if bacia==0:h_aux = hidros.COMPATIBILIDADE_SPT.loc[h]
                    else: h_aux = h
                    stage[j][f'ctr{h}_1'].RHS = res[(res.Iteradores==t)][f"{cen[s]['y'][inp]+1}"].loc[h_aux]
                inp+=1
                stage[j][f'm{j}'].update()
                # stage[j][f'm{j}'].reset()
                stage[j][f'm{j}'].Params.timeLimit=80
                stage[j][f'm{j}'].optimize()
                # print('Number of constraints: %d' % stage[j][f'm{j}'].NumConstrs);
                # print('Number of continuous variables: %d' % (stage[j][f'm{j}'].NumVars));
                # tempo_aux[j][s][0] = stage[j][f'm{j}'].runtime
                # tempo_aux[j][s][0] = stage[j][f'm{j}'].runtime
                tempo_aux[j][s][0][0] = stage[j][f'm{j}'].runtime
                tempo_aux[j][s][0][1] = stage[j][f'm{j}'].NumConstrs
                tempo_aux[j][s][0][2] = stage[j][f'm{j}'].NumVars - stage[j][f'm{j}'].NumBinVars
                tempo_aux_aux[j-1] = stage[j][f'm{j}'].runtime
                # print(stage[j][f'm{j}'].runtime)
                # tempo_aux.append(time.time()-tempo_ini)
                status = stage[j][f'm{j}'].status
                if status == gb.GRB.INFEASIBLE:
                    stage[j][f'm{j}'].computeIIS()
                    stage[j][f'm{j}'].write("model.ilp")
                    print('forward infeasible')
                list_st.append(status)
                # stage[j][f'm{j}'].write(os.path.join(path1,f'SDDP_{j}{s}v4.lp'))
                # stage[j][f'm{j}'].write(os.path.join(path1,f'SDDP_{j}{s}v4.sol'))
                # print(stage[j]['v2_1'].x)
                # list_vol = [stage[j][f'v{h}_1'].x for h in reservatorios]
                # negvol = [num for num in list_vol if num < 0]
                # print(j,negvol)
                for h in reservatorios:
                    if g.loc[j][0]>1:cen[s][f'vsult{h}'].append(stage[j]['v{0}_{1}'.format(h,cen[s]['no'][inp])].x);
                    else:cen[s][f'vsult{h}'].append(stage[j][f'v{h}_1'].x);
                for h in h_agreg:
                    if g.loc[j][0]>1:
                        cen[s][f'ysult{h}_{q}'].append(stage[j]['yn{0}_{1}'.format(h,cen[s]['no'][inp-q])].x);
                    else:
                        cen[s][f'ysult{h}_{q}'].append(stage[j][f'yn{h}_1'].x);
                            
                # for b in hidros.BACIA.unique():
                #     if g.loc[j][0]>1:
                #         for q in range(ordem):
                #             cen[s][f'ysult{b}_{q}'].append(stage[j]['y{0}_{1}'.format(b,cen[s]['no'][inp-q])].x);
                #     else:
                #         for q in range(ordem):
                #             cen[s][f'ysult{b}_{q}'].append(stage[j]['y{0}_{1}'.format(b,1)].x);
                            
                if simu == 1:
                    cen[s],t2 = salva_var(cen[s],stage[j],t2,g.loc[j][0],nos_per_aux)
                    t2 += relativedelta(months=+1)
                cen[s]['corte'] += sum(sum(sum(stage[j][f'd{sub}_{l}_{n}'].x for sub in range(1,SUBS+1)) for l in range(1,PAT+1)) for n in tree1.keys())
                if mc ==1: cen[s]['custo_simu'] += copy.deepcopy(stage[j][f'm{j}'].objVal) - copy.deepcopy(sum(sum(stage[j][f'alfa_{i}_{a}'].x for a in range(abert)) for i in range(ini,fini)))/abert**(g.loc[j][0]-1)
                else: cen[s]['custo_simu'] += copy.deepcopy(stage[j][f'm{j}'].objVal) - copy.deepcopy(sum(stage[j][f'alfa_{i}'].x for i in range(ini,fini)))/abert**(g.loc[j][0]-1)
                if g.loc[j][0]>1: inp+=1
                t+=relativedelta(months=+g.loc[j][0])

    else:         
        for s in cen:
            cen[s]['phsult'] = [];cen[s]['gtsult'] = [];cen[s]['CMO'] = [];cen[s]['dfc'] = []; cen[s]['vsult'] = []; cen[s]['corte'] = 0; cen[s]['custo_simu'] = 0;
            for h in hidros.index:
                if hidros.loc[h]['USINA_FIO_DAGUA']==0:
                    cen[s][f'vsult{h}'] = []
                cen[s][f'ysult{h}'] = []
        for s in cen:
            inp = g.loc[j][0];t=start_month;
            for j in range(1,stages):
                tree1 = tree(int(g.loc[j][0]),abert);
                ini = nos_per_aux[g.loc[j][0]-1]-abert**(g.loc[j][0]-1)+1
                fini = nos_per_aux[g.loc[j][0]-1]+1
                for n in tree1.keys():
                    
                    for h in hidros.index:
                        if n==1: 
                            if j==1:
                               stage[j][f'ctr{h}_{n}'].RHS = res[h][1].loc[t].iloc[0];
                            else: stage[j][f'ctr{h}_{n}'].RHS = res[h][s+1].loc[t].iloc[0];
                        else: stage[j][f'ctr{h}_{n}'].RHS = res[h][s+1].loc[t].iloc[0];
                    t+=relativedelta(months=+1); 
                if j==1:
                    for h in hidros.index:
                        if hidros.loc[h]['USINA_FIO_DAGUA']==0:
                            stage[j][f'ctr_{h}'].RHS = hidros.loc[h]['VOLUME_INICIAL']*(hidros.loc[h]['VOLUME_MAXIMO_OPERACIONAL']-hidros.loc[h]['VOLUME_MINIMO_OPERACIONAL'])/100
                        for q in range(ordem):
                            stage[j][f'ctr_y{h}_{q}'].RHS =  hist_bacia.iloc[q-1][f'{h}']
                else:
                    for h in hidros.index:
                        for q in range(ordem):
                             stage[j][f'ctr_y{h}_{q}'].RHS = cen[s][f'ysult{h}_{q}'][j-2]
                        if hidros.loc[h]['USINA_FIO_DAGUA']==0:
                            stage[j][f'ctr_{h}'].RHS = cen[s][f'vsult{h}'][j-2] 
                stage[j][f'm{j}'].update()
                stage[j][f'm{j}'].optimize()
                # stage[j]['m{0}'.format(j)].write(os.path.join(path1,f'simu_{j}{s}_{mc}_3.lp'))
                # stage[j]['m{0}'.format(j)].write(os.path.join(path1,f'simu_{j}{s}_{mc}_3.sol'))
                for h in hidros.index:
                    if g.loc[j][0] == 1: aux_in = 1 
                    else: aux_in = g.loc[j][0]-1
                    if hidros.loc[h]['USINA_FIO_DAGUA']==0:
                        cen[s][f'vsult{h}'].append(stage[j]['v{0}_{1}'.format(h,aux_in)].x);
                    for q in range(ordem):
                        cen[s][f'ysult{h}_{q}'].append(stage[j]['y{0}_{1}'.format(h,aux_in)].x);
                if simu == 1:
                    for n in tree1.keys():
                        cen[s]['phsult'].append(sum(stage[j]['ph{0}_1_{1}'.format(h,n)].x for h in hidros.index));
                        cen[s]['gtsult'].append(sum(stage[j]['gt{0}_1_{1}'.format(i,n)].x for i in termos.index)); 
                        cen[s]['dfc'].append(stage[j][f'd1_1_{n}'].x); 
                        cen[s]['CMO'].append(stage[j][f'CMO_1_1_{n}'].pi/horas);
                        cen[s]['vsult'].append(sum(stage[j]['v{0}_{1}'.format(h,n)].x for h in hidros[(hidros.USINA_FIO_DAGUA==0)].index));
                        cen[s]['corte'] += sum(sum(stage[j][f'd{sub}_{l}_{n}'].x for sub in range(1,SUBS+1)) for l in range(1,PAT+1))
                        for h in hidros[(hidros.USINA_FIO_DAGUA==0)].index:
                            cen[s][f'vsult_{h}'].append(stage[j]['v{0}_{1}'.format(h,n)].x);
                    if mc==1: cen[s]['custo_simu'] += copy.deepcopy(stage[j][f'm{j}'].objVal) - copy.deepcopy(sum(sum(stage[j][f'alfa_{i}_{a}'].x for a in range(aux_mc)) for i in range(ini,fini)))/aux_mc**(g.loc[j][0]-1)
                    # else: cen[s]['custo_simu'] += copy.deepcopy(stage[j][f'm{j}'].objVal) - copy.deepcopy(sum(stage[j][f'alfa_{i}'].x for i in range(ini,fini)))/abert**(g.loc[j][0]-1)
                    else: cen[s]['custo_simu'] += copy.deepcopy(stage[j][f'm{j}'].objVal) - copy.deepcopy(sum(stage[j][f'alfa_{i}'].x for i in range(ini,fini)))
                
    return cen, stage,list_st,tempo_aux,tempo_aux_aux

