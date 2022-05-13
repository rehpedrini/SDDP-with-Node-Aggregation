# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 16:12:41 2020

@author: rehpe
"""
import pandas as pd
import gurobipy as gb
import numpy as np
import copy
from heapq import nlargest, nsmallest
import numpy_financial as npf
import math as mt
import timeit
from itertools import product
import random as rd
# from datetime import date
import datetime as dt
# from sorteio_v3 import sorteio
import os

# T = 2
# start_month = dt.datetime(2020,11,1)
start_month = dt.datetime(2020,1,1)
# start_month = dt.datetime(2019,8,1)
# start_month = dt.date(2019,8,1)
horas = 1
sev = 0
tax = 0
# tax = 0.008
abert = 2
risk = 0.5
phi  = 0
prob1 = (1/abert)*(1-phi)
prob2 = (1/abert)*(1-phi)+(1/abert)*(phi/risk)
## system data
SUBS = 1;PAT = 1
cz = 3600*24*30/1000000
ordem = 1;seed = 1;ordem = 1;mc = 0;det = 0;fph_flag = 0;sort_dist=0;flag_2per=0;
lhs_flag=1;flag_2020=0;cut_selection=1;fim_percurso = 0;bacia = 0
# import_cuts = 0;simu = 0
path = 'D:\Dropbox\Renata Pedrini (1)\SDDP-Agregada\Dados'
path1 = 'D:\Dropbox\Renata Pedrini (1)\SDDP-Agregada/Gurobi/novo_sistema'
# path2 = 'D:/Renata/Doutorado/DRO/Dados'
#################### Dados FPH #################################
# linhas = [28,28,28,20]
# usinas = ['Campos Novos','Barra Grande', 'Machadinho', 'Itá']
# coef_fph={}
# for h in range(1,NH+1):
#     coef_fph[h] = pd.read_excel(os.path.join(path,'FPH_reduzida.xlsx'), header = None, sheet_name=usinas[h-1], skiprows = [*range(0,2),*range(linhas[h-1],100+linhas[h-1])], usecols  = [*range(0, 3)])
###############################################################

########################################## Leitura de Dados #######################################################################
# file = os.path.join(path,'Sistema Fredo/HIDRELETRICA_v2.csv')
file = os.path.join(path,'Sistema Fredo/HIDRELETRICA.csv')
data = pd.read_csv(file,sep=";",index_col=0,header=0) 
hidros = data[(data.USINA_EM_OPERACAO == 1)]
reservatorios = hidros[(hidros.USINA_FIO_DAGUA==0)].index
fiodagua = hidros[(hidros.USINA_FIO_DAGUA==1)].index
v0 = {h: hidros.loc[h]['VOLUME_INICIAL']*(hidros.loc[h]['VOLUME_MAXIMO_OPERACIONAL']-hidros.loc[h]['VOLUME_MINIMO_OPERACIONAL'])/100 for h in reservatorios}
### TERMELÉTRICAS ############################################################################
file = os.path.join(path,'Sistema Fredo/TERMELETRICA.csv')
# file = os.path.join(path,'Sistema Fredo/TERMELETRICA_v1.csv')
data = pd.read_csv(file,sep=";",index_col=0,header=0)
termos = data[(data.USINA_EM_OPERACAO == 1)]

########### DEMANDA #################################################
# file = os.path.join(path,'Sistema Fredo/DEMANDA_ESTAGIO_v1.csv')
file = os.path.join(path,'Sistema Fredo/DEMANDA_ESTAGIO_v6.csv')
# file = os.path.join(path,'Sistema Fredo/DEMANDA_ESTAGIO_10anos.csv')
# file = os.path.join(path,'Sistema Fredo/DEMANDA_ESTAGIO_toy.csv')
demanda = pd.read_csv(file,sep=";",index_col=0,header=0)
demanda['DATA'] = [dt.datetime(demanda['ANO'][c],demanda['MES'][c],1) for c in demanda.index]
demanda = demanda.drop(columns=['ANO','MES'])
demanda = demanda.set_index('DATA')

########### PROCESSO ESTOCASTICO #################################################
file = os.path.join(path,'Sistema Fredo/VARIAVEL_ALEATORIA_coeficiente_linear_auto_correlacao.csv')
f = pd.read_csv(file,sep=";",index_col=0,header=0).fillna(0)


########### FPH #################################################
## FPH SPT
# file = os.path.join(path,'Sistema Fredo/HIDRELETRICA_FPH.csv')
# fph = pd.read_csv(file,sep=";",index_col=0,header=0)
## FPH Fredo
file = os.path.join(path,'Sistema Fredo/FUNCAO_PRODUCAO_HIDRELETRICA_FPH1.csv')
fph = pd.read_csv(file,sep=";",index_col=0,header=0)

####### produtividade constante ##################################
# file = os.path.join(path,'Sistema Fredo/produtividade_small_system.csv')
file = os.path.join(path,'Sistema Fredo/produtividade.csv')
prod = pd.read_csv(file,sep=";",index_col=0,header=0)

#### HISTÓRICO AFLUÊNCIA ########################################
file = os.path.join(path,'Sistema Fredo/HISTORICO_AFLUENCIA_MES_ANO.csv')
hist = pd.read_csv(file,sep=";",index_col=0,header=0)
file = os.path.join(path,'Sistema Fredo/HISTORICO_AFLUENCIA_MES_ANO_BACIA.csv')
hist_bacia = pd.read_csv(file,sep=";",index_col=0,header=0)
file = os.path.join(path,'Sistema Fredo/COEFICIENTE_PARTICIPACAO.csv')
coef_part = pd.read_csv(file,sep=";",index_col=0,header=0)
#### PATAMAR DE CARGA ###############################################
# file = os.path.join(path,'Sistema Fredo/PATAMAR_ESTAGIO_v1.csv')
# patamar = pd.read_csv(file,sep=";",index_col=0,header=0)
# patamar['DATA'] = [dt.datetime(patamar['ANO'][c],patamar['MES'][c],1) for c in patamar.index]
# patamar = patamar.drop(columns=['ANO','MES'])
# patamar = patamar.set_index('DATA')
file = os.path.join(path,'Sistema Fredo/SUBMERCADO_AttMatrizPremissa_PorPeriodoPorIdPatamarCarga_sempat.csv')
patamar = pd.read_csv(file,sep=";",index_col=0,header=0)
patamar.Iteradores = pd.to_datetime(patamar.Iteradores)
file = os.path.join(path,'Sistema Fredo/DADOS_AttMatrizOperacional_PorPeriodoPorIdPatamarCarga_sempat.csv')
duracao_patamar = pd.read_csv(file,sep=";",index_col=0,header=0)
duracao_patamar['Iteradores'] = pd.to_datetime(duracao_patamar['Iteradores'])
file = os.path.join(path,'Sistema Fredo/SUBMERCADO_PATAMAR_DEFICIT_AttMatrizOperacional_PorPeriodoPorIdPatamarCarga.csv')
cd =  pd.read_csv(file,sep=";",index_col=0,header=0)
cd = cd[(cd.AttMatriz == 'custo')]
cd.Iteradores = pd.to_datetime(cd.Iteradores)
########### PROCESSO ESTOCASTICO #################################################
if bacia == 0:
    file = os.path.join(path,f'Sistema Fredo/ProcessoEstocasticoHidrologico_{seed}/VARIAVEL_ALEATORIA_coeficiente_linear_auto_correlacao.csv')
    f = pd.read_csv(file,sep=";",index_col=0,usecols=[0,2,3],header=0).fillna(0)
    f.Iteradores = pd.to_datetime(f.Iteradores)
    f = f[(f.index.isin(hidros.COMPATIBILIDADE_SPT.values))]
    file = os.path.join(path,f'Sistema Fredo/ProcessoEstocasticoHidrologico_{seed}/VARIAVEL_ALEATORIA_INTERNA_grau_liberdade.csv')
    gl = pd.read_csv(file,sep=";",index_col=0,header=0)
    gl= gl[(gl.index.isin(hidros.COMPATIBILIDADE_SPT.values))]
    file = os.path.join(path,f'Sistema Fredo/ProcessoEstocasticoHidrologico_{seed}/VARIAVEL_ALEATORIA_INTERNA_tendencia_temporal.csv')
    y_hist = pd.read_csv(file,sep=";",index_col=0,header=0)
    y_hist.columns = y_hist.columns[:3].tolist() + pd.to_datetime(y_hist.columns[3:]).tolist()
    y_hist = y_hist[(y_hist.index.isin(hidros.COMPATIBILIDADE_SPT.values))]
    # for i in hidros.COMPATIBILIDADE_SPT.values:
    #     y_hist.loc[i,'novo index'] = hidros[hidros.COMPATIBILIDADE_SPT==i].index[0]
    #     gl.loc[i,'novo index'] = hidros[hidros.COMPATIBILIDADE_SPT==i].index[0]
    #     f.loc[i,'novo index'] = hidros[hidros.COMPATIBILIDADE_SPT==i].index[0]
    # y_hist.reset_index(inplace=True)
    # y_hist = y_hist.set_index('novo index')
    # gl.reset_index(inplace=True)
    # gl = gl.set_index('novo index')
    # f.reset_index(inplace=True)
    # f = f.set_index('novo index')
    h_agreg = hidros.index
else:
    file = os.path.join(path,'Sistema Fredo/ProcessoEstocasticoHidrologico_bacia/VARIAVEL_ALEATORIA_coeficiente_linear_auto_correlacao.csv')
    f = pd.read_csv(file,sep=";",index_col=0,header=0).fillna(0)
    f.Iteradores = pd.to_datetime(f.Iteradores)
    file = os.path.join(path,'Sistema Fredo/ProcessoEstocasticoHidrologico_bacia/VARIAVEL_ALEATORIA_INTERNA_grau_liberdade.csv')
    gl = pd.read_csv(file,sep=";",index_col=0,header=0)
    gl_aux = gl.sum(level=0)
    gl = gl.set_index('idVariavelAleatoriaInterna')
    # gl.columns = gl.columns[:2].tolist() + pd.to_datetime(gl.columns[2:]).tolist()
    file = os.path.join(path,'Sistema Fredo/ProcessoEstocasticoHidrologico_bacia/VARIAVEL_ALEATORIA_INTERNA_tendencia_temporal.csv')
    y_hist = pd.read_csv(file,sep=";",index_col=0,header=0)
    y_hist.columns = y_hist.columns[:3].tolist() + pd.to_datetime(y_hist.columns[3:]).tolist()
    y_hist = y_hist.sum(level=0)
    y_hist = y_hist.drop(y_hist.columns[:2], axis=1)
    # h_agreg = hidros.BACIA.unique()
    h_agreg = np.sort(hidros.BACIA.unique())
    for t in y_hist.columns:
        y_hist[t]=y_hist[t]-gl_aux.valor
######### PERÍODO - ESTÁGIO #############################################
file = os.path.join(path,'Sistema Fredo/DADOS_periodo_estagio.csv')
g = pd.read_csv(file,sep=";",index_col=0,header=0)
T = g.duracao.sum()
# T = 4
# g = [1,2]
# num_d = len(g);
#### RESÍDUOS ###################################################
if bacia == 0:
    file = os.path.join(path,f'Sistema Fredo/ProcessoEstocasticoHidrologico_{seed}/VARIAVEL_ALEATORIA_residuo_espaco_amostral.csv')
    res = pd.read_csv(file,sep=";",index_col=0,header=0)
    res.Iteradores = pd.to_datetime(res.Iteradores)
    res = res[(res.index.isin(hidros.COMPATIBILIDADE_SPT.values))]
    # for i in hidros.COMPATIBILIDADE_SPT.values:
    #     res.loc[i,'novo index'] = hidros[hidros.COMPATIBILIDADE_SPT==i].index[0]
    # res.reset_index(inplace=True)
    # res = res.set_index('novo index')
else:
    file = os.path.join(path,'Sistema Fredo/ProcessoEstocasticoHidrologico_bacia/VARIAVEL_ALEATORIA_residuo_espaco_amostral.csv')
    res = pd.read_csv(file,sep=";",index_col=0,header=0)
    res.Iteradores = pd.to_datetime(res.Iteradores)
    
#########################################################################
# cortes_fim = np.array([4096, 4097, 4098, 4099, 4100, 4101, 4102, 4103, 4104, 4105, 4106, 4107, 4108, 4109, 4110, 4111, 4112, 4113, 4114, 4115, 4116, 4117, 4118, 4119, 4120, 3021, 3057, 3176, 3204, 3212, 3218, 3226, 3246, 3248, 3259, 3273, 3284, 3292, 3299, 3322, 3334, 3349, 3358, 3362, 3363, 3366, 3378, 3380, 3387, 3400, 3413, 3441, 3443, 3451, 3454, 3455, 3458, 3463, 3464, 3466, 3469, 3470, 3472, 3475, 3483, 3485, 3487, 3491, 3492, 3493, 3497, 3502, 3503, 3504, 3505, 3509, 3510, 3512, 3513, 3514, 3517, 3519, 3520, 3522, 3526, 3528, 3531, 3532, 3536, 3537, 3540, 3547, 3553, 3554, 3556, 3559, 3561, 3562, 3567, 3572, 3575, 3578, 3581, 3583, 3584, 3587, 3588, 3594, 3602, 3603, 3604, 3607, 3608, 3611, 3615, 3616, 3619, 3620, 3622, 3623, 3625, 3629, 3632, 3633, 3636, 3637, 3638, 3640, 3642, 3643, 3645, 3646, 3648, 3650, 3652, 3653, 3655, 3657, 3658, 3660, 3662, 3663, 3664, 3665, 3667, 3668, 3671, 3672, 3673, 3674, 3675, 3676, 3678, 3680, 3682, 3683, 3684, 3686, 3689, 3690, 3692, 3693, 3694, 3697, 3699, 3700, 3702, 3703, 3704, 3706, 3712, 3713, 3714, 3715, 3717, 3718, 3721, 3722, 3723, 3724, 3726, 3727, 3729, 3731, 3733, 3734, 3736, 3737, 3738, 3739, 3740, 3742, 3743, 3744, 3747, 3748, 3749, 3750, 3751, 3752, 3753, 3754, 3755, 3756, 3759, 3760, 3761, 3762, 3763, 3765, 3767, 3771, 3772, 3774, 3776, 3777, 3778, 3779, 3780, 3782, 3783, 3784, 3786, 3787, 3788, 3789, 3790, 3791, 3792, 3794, 3795, 3797, 3799, 3800, 3801, 3803, 3804, 3805, 3807, 3808, 3809, 3810, 3811, 3812, 3813, 3816, 3817, 3818, 3819, 3820, 3821, 3822, 3823, 3824, 3825, 3826, 3827, 3828, 3830, 3831, 3833, 3834, 3835, 3836, 3837, 3838, 3839, 3840, 3843, 3844, 3845, 3847, 3848, 3849, 3850, 3851, 3852, 3853, 3854, 3856, 3857, 3858, 3859, 3860, 3861, 3862, 3863, 3864, 3865, 3866, 3867, 3869, 3870, 3871, 3872, 3873, 3874, 3875, 3876, 3877, 3878, 3879, 3880, 3881, 3882, 3883, 3884, 3885, 3886, 3887, 3888, 3889, 3890, 3891, 3892, 3893, 3894, 3895, 3896, 3897, 3898, 3900, 3902, 3903, 3904, 3905, 3906, 3907, 3908, 3909, 3910, 3911, 3912, 3913, 3914, 3915, 3916, 3917, 3918, 3919, 3920, 3921, 3922, 3923, 3924, 3925, 3926, 3927, 3928, 3929, 3930, 3931, 3932, 3933, 3934, 3935, 3936, 3937, 3938, 3939, 3940, 3941, 3942, 3943, 3944, 3945, 3946, 3947, 3948, 3949, 3950, 3951, 3952, 3953, 3954, 3955, 3956, 3957, 3958, 3959, 3960, 3961, 3962, 3963, 3964, 3965, 3966, 3967, 3968, 3969, 3970, 3971, 3972, 3973, 3974, 3975, 3976, 3977, 3978, 3979, 3980, 3981, 3982, 3983, 3984, 3985, 3986, 3987, 3988, 3989, 3990, 3991, 3992, 3993, 3994, 3995, 3996, 3997, 3998, 3999, 4000, 4001, 4002, 4003, 4004, 4005, 4006, 4007, 4008, 4009, 4010, 4011, 4012, 4013, 4014, 4015, 4016, 4017, 4018, 4019, 4020, 4021, 4022, 4023, 4024, 4025, 4026, 4027, 4028, 4029, 4030, 4031, 4032, 4033, 4034, 4035, 4036, 4037, 4038, 4039, 4040, 4041, 4042, 4043, 4044, 4045, 4046, 4047, 4048, 4049, 4050, 4051, 4052, 4053, 4054, 4055, 4056, 4057, 4058, 4059, 4060, 4061, 4062, 4063, 4064, 4065, 4066, 4067, 4068, 4069, 4070, 4071, 4072, 4073, 4074, 4075, 4076, 4077, 4078, 4079, 4080, 4081, 4082, 4083, 4084, 4085, 4086, 4087, 4088, 4089, 4090, 4091, 4092, 4093, 4094, 4095])