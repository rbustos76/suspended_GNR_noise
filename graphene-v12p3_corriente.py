#
#  Para contruir todo se tomo "graphene-ac_v8_corriente.py" y "graphene-zz_v8_corriente.py"
#  desde: "Dropbox/SiFeRa/Cosas-listas-Chile/iLmax-vs-eps-vs-L(corregido)/eps.eq.0p100/L=26"
#  como base.
#
#  La versión 9 es lo mismo pero se le agregó el bombeo neto por ciclo.
#  La versión 10 es lo mismo pero se le agregó maker desde afuera (sirve par zz y ac)
#  La versión 11 se le agregó tarjeta de entrada y se corrigio error con xmin
#    (ahora da igual que versión 9)
#  La versión 12p3 "Under construction"
#     Se corrigio error con Anderson
#            RandE=(1.-random.random())*noise_delta*V0
#     Se agregó dTdA y \Sum_n T_n (1-T_n) para calcular ruido de Termico y shot noise.
#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Este programa calcula la corriente de bombeo máxima para una nanocinta de
grafeno con borde zz o ac.'''

import numpy as np
import kwant 
import matplotlib.pyplot as plt
import time
from GNR_functions_v0p6 import GNR
import sys as Opsys


if __name__ == "__main__":
# Parámetros fijos del sistema
 E0 = 0.      # Energías de sitio de los leads
 V0 = -2.66   # Hopping
 dEdz = 1.    # Factor dE/dz
 q = 0.       # Factor q
 h = 0.0001   # Incremento para derivadas
 q_el = 1.    # Carga del electrón
 cte = q_el*(2.*np.pi)**(-1)  # Constante para corrientes
 acc = 1./np.sqrt(3.)  # Distancia C-C
 namein='graphene-v12p3_corriente.inp'

# Lectura tarjeta de entrada
 data=GNR.reader(namein,3,' ',True)
 for i in range(len(data)):
    if len(data[i]) != 10:
       print('Mal lectura. Parámetros de entrada = ',len(data[i]), ' para linea ',i+4)
       for j in range(len(data[i])):        
          print('i = ',[i],' j = ',[j],' data[i][j] = ', data[i][j])
       Opsys.exit()
    else:
       print(data[i][0],data[i][1],data[i][2],data[i][3],data[i][4] \
       ,data[i][5],data[i][6],data[i][7],data[i][8],data[i][9])
 print('Leidos ',len(data),' conjunto de parámetros de la tarjeta de entrada: ',namein)
 
# Parámetros que cambian
#tipo Lx Ly nx ny eps nameout    #orden de los parámetros  
 for i in range(len(data)):
   print()
   tipo=data[i][0]
   Lx=int(data[i][1])
   Ly=int(data[i][2])
   nx_max=int(data[i][3])
   ny_max=int(data[i][4])
   eps=float(data[i][5])
   noise_type=data[i][6]
   noise_delta=float(data[i][7])
   noise_seed=int(data[i][8])
   nameout=data[i][9]


   p_ran,v_ran=GNR.set_ranE(tipo,Lx,Ly,noise_type,noise_seed,noise_delta,V0)
#   p_ran=[NLx,NLy,DLx,DLy,Lxmin,Lymin]
#   out01 = open('test-01.dat','w')
#   out01.write(f' NLx = {p_ran[0]:<8.4G} NLy = {p_ran[1]:<8.4G} \
#DLx = {p_ran[2]:<8.4G} DLy = {p_ran[3]:<8.4G} \
#Lxmin = {p_ran[4]:<8.4G} Lymin = {p_ran[5]:<8.4G} \n')
#   out01.write(' m      l       xmin     xmax     ymin     ymax     v_ran \n')
#  for m in range(0,p_ran[0]):
#      xmin=p_ran[4]+m*p_ran[2]
#      xmax=p_ran[4]+(m+1)*p_ran[2]
#      for l in range(0,p_ran[1]):
#         ymin=p_ran[5]+l*p_ran[3]
#         ymax=p_ran[4]+(l+1)*p_ran[3]
#         out01.write(f' {m:<6d} {l:<6d} {xmin:<8.4G} {xmax:<8.4G} {ymin:<8.4G} {ymax:<8.4G} {v_ran[m][l]:<8.4G} \n')
#   out01.close()
 #out03 = open('iLmax-ac_L='+str(Lx)+'_eps=0p1_v9.dat','w')

   if nameout=='default':
      nameout='iLmax-v12p3-'+tipo+'-'+str(Lx)+'-'+str(Ly)+'-'+str(nx_max)+'-'+str(ny_max) \
          +'-'+str(eps)+'-'+noise_type+'-'+str(noise_delta)+'-'+str(noise_seed)+'.dat'
   print(nameout)

#   paraelcodigo

   out03 = open(nameout,'w')
   out03.write(f'# tipo={tipo:<3s}| Lx={Lx:<5d}| Ly={Ly:<5d}| nx_max={nx_max:<3d}| ny_max={nx_max:<3d}| eps={eps:<8.4G} \n')
   out03.write(f'# noise_type={noise_type:<9s}| noise_seed={noise_seed:<10d}| noise_delta={noise_delta:<8.4G} \n')
#   out03.write(f'#  E0 = {E0:<8.4G}| V0 = {V0:<8.4G} | dEdz = {dEdz:<8.4G}| q = {q:<8.4G}| h = {h:<8.4G} \n')
   out03.write(f'#{"nx":<4s} {"ny":<4s} {"iL_max ":>13s} {"iL_pump ":>13s} ')
   out03.write(f'{"Transmittance ":>13s} {"Reflectance ":>13s} {"dTdq":>13s} {"sum T*(1-T)":>13s}\n')
   xmin=-0.1
   xmax=Lx+0.1
   if tipo=='ac':
      ymin=-0.1
      ymax=Ly+0.1
   elif tipo=='zz':
      ymin=-0.6
      ymax=Ly+0.6
   else:
      print('mal tipo en main')
      Opsys.exit()   
 #kwant.plot(GNR.maker(E0,V0,dEdz,0.1,1,0,Lx,Ly,xmax,xmin,ymax,ymin,tipo),dpi=300)
 #Opsys.exit()
 # Cálculo de la corriente máxima
   print('Energía:',eps)
   print('Tamaño del sistema:')
   print('    Lx =',Lx,'a = ',Lx*2.46,'Angstroms = ',Lx*2.46/10000.,'micras')
   print('    Ly =',Ly,'a = ',Ly*2.46,'Angstroms = ',Ly*2.46/10000.,'micras')
 # Barrido en modos
   for nx in range(1,nx_max+1):
     for ny in range(0,ny_max+1):
       start = time.time()
       param_maker=[E0,V0,dEdz,q,nx,ny,Lx,Ly,xmax,xmin,ymax,ymin,tipo,p_ran,v_ran]
       S0 = kwant.smatrix(GNR.maker(param_maker),eps)
       param_maker[3]=q+h
       Sq_f = kwant.smatrix(GNR.maker(param_maker),eps)
       param_maker[3]=q-h
       Sq_i = kwant.smatrix(GNR.maker(param_maker),eps)
       param_maker[3]=q+h
       param_maker[4]=0
       param_maker[5]=0
       SE_f = kwant.smatrix(GNR.maker(param_maker),eps)
       param_maker[3]=q-h
       param_maker[4]=0
       param_maker[5]=0
       SE_i = kwant.smatrix(GNR.maker(param_maker),eps)
       iL_max=0.
       iL_pump=0.
       modS=0.
       Trans=0.
       T1mT=0.
       dTdq=0.
       Refle=0.
       # Barrido en número de leads
       for k in range(2):
       # Matrices de las cuales necesitamos los elementos
          sub = S0.submatrix(0,k) # Lk
          subq_f=Sq_f.submatrix(0,k)
          subq_i=Sq_i.submatrix(0,k)
          dSdqsub = (subq_f-subq_i)/(2*h)
          dSdEsub = (SE_f.submatrix(0,k)-SE_i.submatrix(0,k))/(2*h)
          # Barrido en las componentes internas de las submatrices

#	calculo transmitancia!
          for i in range(len(sub)):
             for j in range(len(np.transpose(sub))): # Suponiendo que todos los leads tienen el mismo número de canales
             # Transmitancia y reflectancia
                if k==1:  #solo cuenta el sub-bloque 0-1
                   Tn=np.real(sub[i,j]*np.conj(sub[i,j]))
                   Trans += Tn
                   T1mT+= Tn*(1.-Tn)
                   Tnq_f=np.real(subq_f[i,j]*np.conj(subq_f[i,j]))
                   Tnq_i=np.real(subq_i[i,j]*np.conj(subq_i[i,j]))
                   dTdq  += (Tnq_f-Tnq_i)/(2*h)
                else:
                   Refle += np.real(sub[i,j]*np.conj(sub[i,j]))
             # Valor máximo de la corriente de orden 1
                iL_max  += cte*np.imag((dSdqsub[i,j]*np.conj( sub[i,j])))
             # Corriente de bombeo
                iL_pump += cte*np.imag((dSdqsub[i,j]*np.conj( dSdEsub[i,j] )))
#                modS += cte*abs(dSdqsub[i,j])
       end = time.time()
       print(f'nx = {nx:<4d} | ny = {ny:<4d} | iLmax = {iL_max:<13.6G} | iLpump = {iL_pump:<13.6G} | Trans={Trans:<13.6G} | tiempo = {end-start:<13.6G}')
       out03.write(f'{nx:^4d} {ny:^4d} {iL_max:>13.6G} {iL_pump:>13.6G} ')
       out03.write(f'{Trans:>13.6G} {Refle:>13.6G} {dTdq:>13.6G} {T1mT:>13.6G}\n')
   out03.close()
