# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 14:15:43 2020

@author: rehpe
"""
from param import *
# def modelo(mc,g,f,t_simu): 
def modelo(mc,g,f,t_simu,import_cuts): 
    from tree import tree
    from dateutil.relativedelta import relativedelta
    import datetime as dt
    stages = len(g)
    nos_per=[abert**(i) for i in range(int(g.loc[2][0]))]
    nos_per_aux=[];aux = 0
    for i in nos_per:
        aux += i
        nos_per_aux.append(aux) 

    if import_cuts ==1:
        if flag_2per==0:
            cortes = pd.read_csv(os.path.join(path1,f'cortes_{stages}_{seed}_{phi}_{sev}.csv'),sep = ';',index_col = 0)
        else: cortes = pd.read_csv(os.path.join(path1,f'cortes_{T}_{seed}.csv'),sep = ';',index_col = 0)
    
    obj = 0; stage = {};t=start_month;ind_aux=0;cut={}
    for j in range(1,len(g)+1):
        cut[j]=[]
        tree1 = tree(int(g.loc[j][0]),abert)
        stage[j]={}; 
        stage[j][f'm{j}'] = gb.Model(f'PDDE_{j}')
        
        for h in reservatorios:
            stage[j][f'vp{h}'] = stage[j][f'm{j}'].addVar(lb=0,vtype=gb.GRB.CONTINUOUS, name=f'vp{h}')
            stage[j][f'ctr_{h}'] = stage[j][f'm{j}'].addLConstr(stage[j][f'vp{h}'] == 0)
            
        for h in h_agreg:
            for o in range(ordem):
                stage[j][f'yp{h}_{o}'] = stage[j][f'm{j}'].addVar(lb = -gb.GRB.INFINITY,vtype=gb.GRB.CONTINUOUS, name=f'yp{h}_{o}')
            # stage[j][f'yf{h}'] = stage[j][f'm{j}'].addVar(lb=0,vtype=gb.GRB.CONTINUOUS, name=f'yf{h}')

        ## Variáveis
        ini = nos_per_aux[int(g.loc[j][0])-1]-abert**(int(g.loc[j][0])-1)+1
        fini = nos_per_aux[int(g.loc[j][0])-1]+1
        if mc==1:
            if import_cuts==1: aux_mc = len(cortes[(cortes['o']==0) & (cortes['j']-1<=1) & (cortes['s']==0)])
            else: aux_mc = 0
        else: aux_mc=0
        for a in range(ini,fini):
            if mc == 1:
                if det == 1 and simu==1:
                    for i in range(aux_mc):
                        stage[j][f'alfa_{a}_{i}'] = stage[j][f'm{j}'].addVar(lb=0,  vtype=gb.GRB.CONTINUOUS, name=f'alfa_{a}_{i}')
                else:
                    for i in range(abert):
                        stage[j][f'alfa_{a}_{i}'] = stage[j][f'm{j}'].addVar(lb=0,  vtype=gb.GRB.CONTINUOUS, name=f'alfa_{a}_{i}')
                
            else: stage[j][f'alfa_{a}'] = stage[j][f'm{j}'].addVar(lb=0,  vtype=gb.GRB.CONTINUOUS, name=f'alfa_{a}')
        ac=1;cont=0;t_aux = 1
        for n in tree1.keys():
            cont+=1; ac+=1;
            for l in range(1,PAT+1):
                for i in termos.index:
                        stage[j][f'gt{i}_{l}_{n}'] = stage[j][f'm{j}'].addVar(lb=0, ub=termos.loc[i]['POTENCIA_MAXIMA']*duracao_patamar[(duracao_patamar.Iteradores == t)][f'{l}'].iloc[0], vtype=gb.GRB.CONTINUOUS, name=f'gt{i}_{l}_{n}') 
                for sub in range(1,SUBS+1):
                    valor_demanda = demanda[t==demanda.index][f'{sub}'].iloc[0]*duracao_patamar[(duracao_patamar.Iteradores == t)][f'{l}'].iloc[0]
                    stage[j][f'd{sub}_{l}_{n}'] = stage[j][f'm{j}'].addVar(lb=0, ub=valor_demanda, vtype=gb.GRB.CONTINUOUS, name=f'd{sub}_{l}_{n}') 
            
            for h in h_agreg: stage[j][f'yn{h}_{n}'] = stage[j][f'm{j}'].addVar(lb = -gb.GRB.INFINITY, vtype=gb.GRB.CONTINUOUS, name=f'yn{h}_{n}')
            for h in hidros.index:
                stage[j][f'yf{h}_{n}'] = stage[j][f'm{j}'].addVar(lb=0, vtype=gb.GRB.CONTINUOUS, name=f'yf{h}_{n}')
                stage[j][f'y{h}_{n}'] = stage[j][f'm{j}'].addVar(lb = -gb.GRB.INFINITY, vtype=gb.GRB.CONTINUOUS, name=f'y{h}_{n}')
                if hidros.loc[h]['USINA_FIO_DAGUA']==0:
                    stage[j][f'v{h}_{n}'] = stage[j][f'm{j}'].addVar(lb=0, ub=hidros.loc[h]['VOLUME_MAXIMO_OPERACIONAL']-hidros.loc[h]['VOLUME_MINIMO_OPERACIONAL'], vtype=gb.GRB.CONTINUOUS, name=f'v{h}_{n}')
                for l in range(1,PAT+1):
                    stage[j][f's{h}_{l}_{n}'] = stage[j][f'm{j}'].addVar(lb=0,  vtype=gb.GRB.CONTINUOUS, name=f's{h}_{l}_{n}')
                    stage[j][f'ph{h}_{l}_{n}'] = stage[j][f'm{j}'].addVar(lb=0, ub=hidros.loc[h]['POTENCIA_MAXIMA'],vtype=gb.GRB.CONTINUOUS, name=f'ph{h}_{l}_{n}')
                    stage[j][f'q{h}_{l}_{n}'] = stage[j][f'm{j}'].addVar(lb=0, ub=hidros.loc[h]['TURBINAMENTO_MAXIMO'], vtype=gb.GRB.CONTINUOUS, name=f'q{h}_{l}_{n}')
    
                # stage[j][f'yc{h}_{n}'] = stage[j][f'm{j}'].addVar(lb=-1000,ub = 0, vtype=gb.GRB.CONTINUOUS, name=f'yc{h}_{n}')
            ## Restrições
            if n==1:
                for h in h_agreg:
                    if bacia==0:
                        compat = hidros.COMPATIBILIDADE_SPT.loc[h]
                        stage[j][f'ctr{h}_{n}'] = stage[j][f'm{j}'].addLConstr(stage[j][f'yn{h}_{n}'] - gb.quicksum(f[(f.Iteradores==t)].loc[compat][f'{o+1}']*stage[j][f'yp{h}_{o}'] for o in range(ordem))  == 0)
                    else:
                        stage[j][f'ctr{h}_{n}'] = stage[j][f'm{j}'].addLConstr(stage[j][f'yn{h}_{n}'] - gb.quicksum(f[(f.Iteradores==t)].loc[h][f'{o+1}']*stage[j][f'yp{h}_{o}'] for o in range(ordem))  == 0)

                for h in hidros.index:                  
                    if hidros.loc[h]['USINA_FIO_DAGUA']==1:
                        for l in range(1,PAT+1):
                            stage[j][f'm{j}'].addLConstr(stage[j][f'q{h}_{l}_{n}'] + stage[j][f's{h}_{l}_{n}'] - gb.quicksum(stage[j]['q{0}_{1}_{2}'.format(jus,l,n)] + stage[j]['s{0}_{1}_{2}'.format(jus,l,n)] for jus in hidros.JUSANTE_TURBINAS[hidros.JUSANTE_TURBINAS == h].index.tolist()) - stage[j][f'y{h}_{n}']  - stage[j][f'yf{h}_{n}'] == 0)                        
                    else:  stage[j][f'm{j}'].addLConstr(stage[j][f'v{h}_{n}']/cz + gb.quicksum(duracao_patamar[(duracao_patamar.Iteradores == t)][f'{l}'].iloc[0]*(stage[j][f'q{h}_{l}_{n}'] + stage[j][f's{h}_{l}_{n}'] - sum(stage[j]['q{0}_{1}_{2}'.format(jus,l,n)] + stage[j]['s{0}_{1}_{2}'.format(jus,l,n)] for jus in hidros.JUSANTE_TURBINAS[hidros.JUSANTE_TURBINAS == h].index.tolist())) for l in range(1,PAT+1)) - stage[j][f'vp{h}']/cz - stage[j][f'y{h}_{n}'] - stage[j][f'yf{h}_{n}'] == 0)                        
                    ### fph
                    if fph_flag==1:
                        for l in range(1,PAT+1):
                            for i in range(fph.PLANO[h].count()):
                                if hidros.loc[h]['USINA_FIO_DAGUA']==0:
                                      stage[j][f'm{j}'].addLConstr(stage[j][f'ph{h}_{l}_{n}'] - (stage[j][f'v{h}_{n}'] + stage[j][f'vp{h}'])*fph.V[h].iloc[i]/2  - fph.Q[h].iloc[i]*stage[j][f'q{h}_{l}_{n}'] - fph.S[h].iloc[i]*stage[j][f's{h}_{l}_{n}'] <= 0)

                                else: stage[j][f'm{j}'].addLConstr(stage[j][f'ph{h}_{l}_{n}'] - fph.Q[h].iloc[i]*stage[j][f'q{h}_{l}_{n}'] - fph.S[h].iloc[i]*stage[j][f's{h}_{l}_{n}'] <= +fph.RHS[h].iloc[i])
    
            else:
                for h in hidros.index:
                    # stage[j][f'm{j}'].addLConstr(stage[j][f'y{h}_{n}'] - stage[j][f'yn{h}_{n}']  ==  -gl.valor.loc[hidros.COMPATIBILIDADE_SPT.loc[h]])
                    if hidros.loc[h]['USINA_FIO_DAGUA']==1:
                        for l in range(1,PAT+1):
                            stage[j][f'ctr_{h}_{n}'] = stage[j][f'm{j}'].addLConstr(stage[j][f'q{h}_{l}_{n}'] + stage[j][f's{h}_{l}_{n}'] - gb.quicksum(stage[j]['q{0}_{1}_{2}'.format(jus,l,n)] + stage[j]['s{0}_{1}_{2}'.format(jus,l,n)] for jus in hidros.JUSANTE_TURBINAS[hidros.JUSANTE_TURBINAS == h].index.tolist()) - stage[j][f'y{h}_{n}']  - stage[j][f'yf{h}_{n}']== 0)
                    else: stage[j][f'ctr_{h}_{n}'] = stage[j][f'm{j}'].addLConstr(stage[j][f'v{h}_{n}']/cz + gb.quicksum(duracao_patamar[(duracao_patamar.Iteradores == t)][f'{l}'].iloc[0]*(stage[j][f'q{h}_{l}_{n}'] + stage[j][f's{h}_{l}_{n}'] - gb.quicksum(stage[j]['q{0}_{1}_{2}'.format(jus,l,n)] + stage[j]['s{0}_{1}_{2}'.format(jus,l,n)] for jus in hidros.JUSANTE_TURBINAS[hidros.JUSANTE_TURBINAS == h].index.tolist()))for l in range(1,PAT+1)) \
                                                                                   - stage[j]['v{0}_{1}'.format(h,tree1[n]['pais'][0])]/cz - stage[j][f'y{h}_{n}'] - stage[j][f'yf{h}_{n}'] == 0)
                    
                    
                    ### fph
                    if fph_flag==1:
                        for l in range(1,PAT+1):
                            for i in range(fph.PLANO[h].count()):
                                if hidros.loc[h]['USINA_FIO_DAGUA']==1:
                                    stage[j][f'm{j}'].addLConstr(stage[j][f'ph{h}_{l}_{n}'] - fph.Q[h].iloc[i]*stage[j][f'q{h}_{l}_{n}'] - fph.S[h].iloc[i]*stage[j][f's{h}_{l}_{n}'] <= +fph.RHS[h].iloc[i])
                                else: stage[j][f'm{j}'].addLConstr(stage[j][f'ph{h}_{l}_{n}'] - (stage[j][f'v{h}_{n}'] + stage[j]['v{0}_{1}'.format(h,tree1[n]['pais'][0])])*fph.V[h].iloc[i]/2  - fph.Q[h].iloc[i]*stage[j][f'q{h}_{l}_{n}'] - fph.S[h].iloc[i]*stage[j][f's{h}_{l}_{n}'] <= 0)
                                
                for h in h_agreg:
                    if bacia==1:
                        if ordem < g.loc[j][0]:stage[j][f'ctr{h}_{n}'] = stage[j][f'm{j}'].addLConstr(stage[j][f'yn{h}_{n}'] - gb.quicksum(f[(f.Iteradores==t)].loc[h][f'{o+1}']*stage[j]['yn{0}_{1}'.format(h,tree1[n]['pais'][o])] for o in range(ordem))== res[(res.Iteradores==t)][f'{ac}'].loc[h])
                        else: stage[j][f'ctr{h}_{n}'] = stage[j][f'm{j}'].addLConstr(stage[j][f'yn{h}_{n}'] - gb.quicksum(f[(f.Iteradores==t)].loc[h][f'{o+1}']*stage[j]['yn{0}_{1}'.format(h,tree1[n]['pais'][o])] for o in range(g.loc[j][0]))\
                              - gb.quicksum(f[(f.Iteradores==(t-relativedelta(months=+o+1)).month)].loc[h]['{0}'.format(o+1)]*stage[j]['yp{0}_{1}'.format(h,o-g.loc[j][0])] for o in range(ordem-g.loc[j][0]))== res[(res.Iteradores==t)][f'{ac}'].loc[h])
                    else:
                        compat = hidros.COMPATIBILIDADE_SPT.loc[h]
                        if ordem < g.loc[j][0]:stage[j][f'ctr{h}_{n}'] = stage[j][f'm{j}'].addLConstr(stage[j][f'yn{h}_{n}'] - gb.quicksum(f[(f.Iteradores==t)].loc[compat][f'{o+1}']*stage[j]['yn{0}_{1}'.format(h,tree1[n]['pais'][o])] for o in range(ordem))== res[(res.Iteradores==t)][f'{ac}'].loc[compat])
                        else: stage[j][f'ctr{h}_{n}'] = stage[j][f'm{j}'].addLConstr(stage[j][f'yn{h}_{n}'] - gb.quicksum(f[(f.Iteradores==t)].loc[compat][f'{o+1}']*stage[j]['yn{0}_{1}'.format(h,tree1[n]['pais'][o])] for o in range(g.loc[j][0]))\
                              - gb.quicksum(f[(f.Iteradores==(t-relativedelta(months=+o+1)).month)].loc[h]['{0}'.format(o+1)]*stage[j]['yp{0}_{1}'.format(h,o-g.loc[j][0])] for o in range(ordem-g.loc[j][0]))== res[(res.Iteradores==t)][f'{ac}'].loc[compat])
            for h in hidros.index:
                compat = hidros.COMPATIBILIDADE_SPT.loc[h]
                if bacia==0: stage[j][f'm{j}'].addLConstr(stage[j][f'y{h}_{n}'] - stage[j][f'yn{h}_{n}']  ==  -gl.valor.loc[compat])
                else: stage[j][f'm{j}'].addLConstr(stage[j][f'y{h}_{n}'] - stage[j][f'yn{hidros.loc[h].BACIA}_{n}']*coef_part[(coef_part['MES']==t.month)][f'{h}']  ==  -gl.valor.loc[compat])
            for sub in range(1,SUBS+1):
                for l in range(1,PAT+1):
                    valor_demanda = demanda[t==demanda.index][f'{sub}'].iloc[0]*patamar[(patamar.Iteradores == t)][f'{l}'].loc[sub]*duracao_patamar[(duracao_patamar.Iteradores == t)][f'{l}'].iloc[0]
                    stage[j][f'CMO_{sub}_{l}_{n}'] = stage[j][f'm{j}'].addLConstr(gb.quicksum(stage[j][f'gt{i}_{l}_{n}'] for i in termos[(termos.SUBSISTEMA == sub)].index) + duracao_patamar[(duracao_patamar.Iteradores == t)][f'{l}'].iloc[0]*gb.quicksum(stage[j][f'ph{h}_{l}_{n}'] for h in hidros[(hidros.SUBSISTEMA == sub)].index) + stage[j][f'd{sub}_{l}_{n}'] == valor_demanda) 
            # produtibilidade constante
            if fph_flag==0:
                for l in range(1,PAT+1):
                    for h in hidros.index:
                        stage[j][f'm{j}'].addLConstr(stage[j][f'ph{h}_{l}_{n}'] == prod['0'][h]*stage[j][f'q{h}_{l}_{n}']) 
            if ac>=abert: ac=0  
            if cont>=abert**(t_aux-1):cont=0;t+=relativedelta(months=+1);t_aux+=1;ac=0

        if fim_percurso == 1 and j==61:
           cortes_fim_data = pd.read_csv(os.path.join(path,'Sistema Fredo/cortes_fim.csv'),sep = ';',index_col = 0)
           for k in cortes_fim:
               for a in range(ini,fini):
                   stage[j][f'm{j}'].addConstr(stage[j][f'alfa_{a}'] >= cortes_fim_data['custo'][k]\
                          +sum(stage[j]['v{0}_{1}'.format(h,a)]*cortes_fim_data[f'PI{h}'][k] for h in reservatorios))
        if import_cuts == 1:
            if flag_2per == 0:
                if mc==1:
                    for k in cortes.index:
                        if cortes['j'][k]==j: 
                            if cortes['tempo'][k]<t_simu:
                                for a in range(ini,fini):
                                    i=cortes['i'][k]
                                    stage[j][f'm{j}'].addConstr(stage[j][f'alfa_{a}_{i}'] >= cortes['custo'][k]\
                                            +sum(stage[j]['v{0}_{1}'.format(h,a)]*cortes[f'PI{h}'][k] for h in reservatorios))
                                            # +sum(stage[j]['yn{0}_{1}'.format(b,a)]*cortes[f'PIy{b}'][k] for b in hidros.BACIA.unique()))
                else:
                    for k in cortes.index:
                        if cortes['j'][k]==j: 
                            if cortes['tempo'][k]<t_simu:
                                for a in range(ini,fini):
                                    stage[j][f'm{j}'].addConstr(stage[j][f'alfa_{a}'] >= cortes['custo'][k]\
                                            +sum(stage[j][f'v{h}_{a}']*cortes[f'PI{h}'][k] for h in reservatorios)\
                                            +sum(stage[j][f'yn{h}_{a}']*cortes[f'PIy{h}{q}'][k] for q in range(ordem) for h in hidros.index))
            else: 
                ind_aux += g.loc[j][0]
                if mc==1:
                    for k in cortes.index:
                        if cortes['j'][k]==ind_aux: 
                            if cortes['tempo'][k]<t_simu:
                                for a in range(ini,fini):
                                    i=cortes['i'][k]
                                    stage[j][f'm{j}'].addConstr(stage[j][f'alfa_{a}_{i}'] >= cortes['custo'][k]\
                                            +sum(stage[j][f'v{h}_{a}']*cortes[f'PI{h}'][k] for h in reservatorios)\
                                            +sum(stage[j][f'yn{h}_{a}']*cortes[f'PIy{h}{q}'][k] for q in range(ordem) for h in hidros.index))
                else:
                    for k in cortes.index:
                        if cortes['j'][k]==ind_aux: 
                            if cortes['tempo'][k]<t_simu:
                                for a in range(ini,fini):
                                    stage[j][f'm{j}'].addConstr(stage[j][f'alfa_{a}'] >= cortes['custo'][k]\
                                            +sum(stage[j][f'v{h}_{a}']*cortes[f'PI{h}'][k] for h in reservatorios)\
                                            +sum(stage[j][f'yn{h}_{a}']*cortes[f'PIy{h}{q}'][k] for q in range(ordem) for h in hidros.index))
                                        
        obj = 0; cont_2=0; t_aux2=1; prob=1;t2=start_month
        for n in tree1.keys():
            cont_2+=1;
            # print(t_aux2,t2)
            obj += prob*horas*npf.npv(tax,(t_aux2-1)*[0]+[gb.quicksum(gb.quicksum(termos.loc[i]['CUSTO']*stage[j][f'gt{i}_{l}_{n}'] for i in termos.index) \
                    + gb.quicksum(stage[j][f'd{sub}_{l}_{n}']*cd.loc[sub][(cd.Iteradores == t2)][f'{l}'] for sub in range(1,SUBS+1))for l in range(1,PAT+1))\
                    + gb.quicksum(stage[j][f'yf{h}_{n}']*5250 for h in hidros.index)])
            if cont_2>=abert**(t_aux2-1):cont_2=0;t_aux2+=1;prob=prob/abert;t2+=relativedelta(months=+1)
        if mc ==1: 
            if det == 1 and simu==1:obj += npf.npv(tax,(t_aux2-1)*[0]+[gb.quicksum(gb.quicksum(stage[j][f'alfa_{a}_{i}'] for a in range(ini,fini))/aux_mc**(g.loc[j][0]-1) for i in range(aux_mc))])
            else:obj += npf.npv(tax,(t_aux2-1)*[0]+[gb.quicksum(gb.quicksum(stage[j][f'alfa_{i}_{a}'] for i in range(ini,fini))/abert**(g.loc[j][0]-1) for a in range(abert))])
        else:obj += npf.npv(tax,(t_aux2-1)*[0]+[gb.quicksum(stage[j][f'alfa_{i}'] for i in range(ini,fini))/abert**(g.loc[j][0]-1)])
        stage[j][f'm{j}'].setObjective(obj, gb.GRB.MINIMIZE) 
        stage[j][f'm{j}'].update()
        stage[j][f'm{j}'].Params.outputflag=0
        stage[j][f'm{j}'].Params.method=1

         ### runtimelimit
        # stage[j][f'm{j}'].Params.FeasibilityTol=10**(-9)
        # stage[j][f'm{j}'].write(f'modelo{j}v3.lp')
        # stage[j][f'm{j}'].write(os.path.join(path1,f'modelo{j}_v3.lp'))
        # print('Number of constraints: %d' % stage[j][f'm{j}'].NumConstrs);
        # print('Number of continuous variables: %d' % (stage[j][f'm{j}'].NumVars - stage[j][f'm{j}'].NumBinVars));
        # print('Number of binary variables: %d' % stage[j][f'm{j}'].NumBinVars);
    return stage,aux_mc