# Para que se vea bien en gnuplot
# splot 'GNR-plotter-v0p0_zz_L=4_modo_1-0.dat' u 3:4:($6==1?$5:1/0)  palette pt 7
# plot 'GNR-plotter-v0p0_zz_L=4_modo_1-0.dat' u 3:4:($6==1?$5:1/0) w image
#
#  Con esta version de maker N interesantes N=(3*l+2)   N=5,8,11,16,17,20,23,26,50,101
#  Con version GNR_plotter_v0p3.py N interesantes N=(3*l+1)   N=4,7,10,...,25
#
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''Este programa calcula la corriente de bombeo máxima para una nanocinta de
grafeno con borde zigzag.'''

import numpy as np
import csv
import kwant 
import matplotlib.pyplot as plt
import matplotlib.backends
import time
import math
import sys as Opsys
import random



class GNR:
#%% Funciones Construcción del sistema
 def frequency(nx,ny,Lx,Ly):
#    '''Devuelve la frecuencia asociada a la combinación de modos nx y ny para
#    una membrana rectangular de lados Lx y Ly.'''
    v = 1.  # Velocidad de propagación en la membrana
    return v*(np.sqrt((nx*np.pi/Lx)**2 + (ny*np.pi/Ly)**2))

#
#  Se agregó la definición de graphene porque era peligroso definirlo por fuera
#  ya que la definición de qué cosa es "graphene" está entrando por afuera de
#  las varibles de entrada. Lo mismo pasa con "a" y "b"
#  Se agregó la posibilidad de que solo modifique energías de sition en un dado parche
#  Se agregó la posibilidad de energías de sitio aleatorias
#
 def maker(param_maker):
   E0,V0,dEdz,q,nx,ny,Lx,Ly,xmax,xmin,ymax,ymin,tipo,p_ran,v_ran=param_maker
   NLx,NLy,DLx,DLy,Lxmin,Lymin=p_ran
#  Región rectangular de scattering
   def region_zz(pos):  # Se añade 0.5 para obtener la cantidad de sitios correcta
       return (0 <= pos[0] <= Lx) and (-0.5 <= pos[1] <= Ly+0.5)  #para zz
   def region_ac(pos):
       return (0 <= pos[0] <= Lx) and (0 <= pos[1] <= Ly)        #para ac
   if tipo == 'zz':
#     Parámetros de red de grafeno_zz:
      graphene = kwant.lattice.general([(1,0),(0.5,np.sqrt(3)*0.5)],  # Vectores primitivos
                                    [(0,0),(0,1/np.sqrt(3))])      # Vectores del motivo
      a,b = graphene.sublattices  # La sublattice "a" se forma moviendo el C en (0,0)
                            # La sublattice "b" se forma moviendo el C en (0,1/sqrt(3))
#     Construccion systema grafeno_zz
      sys = kwant.Builder()
      sys[graphene.shape(region_zz,(Lx/2,Ly/2))] = 0
      sys[kwant.builder.HoppingKind((0,0), a,b)] = V0
      sys[kwant.builder.HoppingKind((0,1), a,b)] = V0
      sys[kwant.builder.HoppingKind((-1,1),a,b)] = V0
#      sys.eradicate_dangling()
      #print('Sitios del sistema antes de añadir leads:  ',len(sites))
      #kwant.plot(sys)  # Ploteo del sistema para controlar
# Construccion Leads grafeno_zz
      sym = kwant.TranslationalSymmetry(graphene.vec((1,0)))
      sym.add_site_family(graphene.sublattices[0],other_vectors = [(-1,2)])
      sym.add_site_family(graphene.sublattices[1],other_vectors = [(-1,2)])
      lead0 = kwant.Builder(sym)
      lead0[graphene.shape(region_zz,(0,0))] = E0
      lead0[kwant.builder.HoppingKind((0,0), a,b)] = V0
      lead0[kwant.builder.HoppingKind((0,1), a,b)] = V0
      lead0[kwant.builder.HoppingKind((-1,1),a,b)] = V0
      lead0.eradicate_dangling()
      lead1 = lead0.reversed()
      lead1.eradicate_dangling()
   elif tipo=='ac':
#     Parámetros de red de grafeno_ac:
      graphene = kwant.lattice.general([(np.sqrt(3)*0.5,0.5),(0,1)],  # Vectores primitivos
                                    [(0,0), (1/np.sqrt(3),0)])      # Vectores del motivo
      a,b = graphene.sublattices  # La sublattice "a" se forma moviendo el C en (0,0)
                            # La sublattice "b" se forma moviendo el C en (0,1/sqrt(3))
#     Construccion systema grafeno_ac
      sys = kwant.Builder()
      sys[graphene.shape(region_ac,(Lx/2,Ly/2))] = 0
      sys[kwant.builder.HoppingKind((0,0), a,b)] = V0
      sys[kwant.builder.HoppingKind((1,-1),a,b)] = V0
      sys[kwant.builder.HoppingKind((1,0), a,b)] = V0
#      sys.eradicate_dangling()
      #print('Sitios del sistema antes de añadir leads:  ',len(sites))
      #kwant.plot(sys)  # Ploteo del sistema para controlar
# Construccion Leads grafeno_ac
      sym = kwant.TranslationalSymmetry(graphene.vec((2,-1)))
      lead0 = kwant.Builder(sym)
      lead0[graphene.shape(region_ac,(0,0))] = E0
      lead0[kwant.builder.HoppingKind((0,0), a,b)] = V0
      lead0[kwant.builder.HoppingKind((1,-1),a,b)] = V0
      lead0[kwant.builder.HoppingKind((1,0), a,b)] = V0
      lead0.eradicate_dangling()
      lead1 = lead0.reversed()
      lead1.eradicate_dangling()
   else:
      print('mal tipo!')
#
# Asignación energías de sitios
#   out02 = open('test-02.dat','w')
   sites=list(sys.sites())
   for site in sites:
       fam = site.family
       pos = site.pos
       tag = site.tag
#   Busca el microparche, la energía aleatoria del mismo y lo suma a la energía de sitio
       m=int((pos[0]-Lxmin)/DLx)
       l=int((pos[1]-Lymin)/DLy)
       if m>=0 and m<NLx and l>=0 and l<NLy:
          RandE=v_ran[m][l]
       else:
          RandE=0
       E_site = E0+RandE
#   Si el átomo está en el parche de entrada cambia las energía de sitio
       if (xmin < pos[0] and pos[0] <= xmax and ymin < pos[1] and pos[1] <= ymax):
         if nx==0 and ny==0:  # Sube todas las energías de sitio por igual
           E_site += dEdz*q
         else:                # Sube las energías de sitio como una onda
           E_site += dEdz*np.sin(nx*np.pi*pos[0]/Lx)*np.cos(ny*np.pi*pos[1]/Ly)*q
#   Ahora sí asigna la energía de sitio correspondiente
       sys[fam(tag[0],tag[1])] = E_site

#       out02.write(f' {m:<6d} {l:<6d} {pos[0]:<8.4G} {pos[1]:<8.4G} {RandE:<8.4G} \n')

   sys.attach_lead(lead0)
   sys.attach_lead(lead1)
   #sites2=list(sys.sites())
   #print('Sitios del sistema después de añadir leads:',len(sites2))
   #kwant.plot(sys)  # Ploteo del sistema+leads para controlar

#   out02.close()
#   Opsys.exit()	# para en cualquier momento
   return sys.finalized()

#  Condition es una función ad-doc que predice cuando SI HAY átomos en el parche (control=True)
#   a veces en el último sitio en "y2 predice que va a haber un átomo que no está
#	igual eso no afecta el resultado (el problema sería alrevés)
 def condition(x,y,Lx,Ly,l,m,tipo):
   if tipo=='zz':
     if x < 0:
        control=False
     elif x > Lx:
        control=False
     elif y < -0.5:
        control=False
     elif y > Ly+0.5:
        control=False
     elif (l-0)%4==0 and (m+1)%4==0:
        control=True
     elif (l-1)%4==0 and (m-1)%4==0:
        control=True
     elif (l-2)%4==0 and (m-1)%4==0:
        control=True
     elif (l-3)%4==0 and (m+1)%4==0:
        control=True
     else:
        control=False
   elif tipo=='ac':
     if x < 0:
        control=False
     elif x > Lx:
        control=False
     elif y < 0:
        control=False
     elif y > Ly:
        control=False
     elif (m-1)%4==0 and (l-2)%4==0:
        control=True
     elif (m-2)%4==0 and (l-2)%4==0:
        control=True
     elif (m-3)%4==0 and (l-4)%4==0:
        control=True
     elif (m-4)%4==0 and (l-4)%4==0:
        control=True
     else:
        control=False
   else:
     print('Mal tipo en condition')
     Opsys.exit()
   return control 

#
#	Función que sirve para leer archivos con cualquier formato..
#
 def reader(name,nskip,string,nonull):
   with open(name) as file:
      for i in range(nskip):
         next(file)
      csvreader = csv.reader(file,delimiter=string)	#elementos separados por string
      datos=[]
      for row in csvreader:
         row_clean = [x for x in row if x != '']	#solo guarda elementos no nulos
         datos.append(row_clean)
   if nonull:
      datos_clean = [ele for ele in datos if ele != []]
      return datos_clean # Devuelve una lista de listas sin elementos nulos, con solo strings.
   else:
      return datos # Devuelve una lista de listas, con solo strings.


#
#  Me genera arreglos de parches donde solo entra un átomo
#  me dice la coordenada del parche
#  me dice si hay átomo o no en el parche
#  Si hay átomo le asigna una energía de sitio aleatoria
#
 def set_ranE(tipo,Lx,Ly,noise_type,noise_seed,noise_delta,V0):
   random.seed(noise_seed)
   if tipo=='zz':	#calcula parche optimo para GNR_zz
      DLy= 0.866025404/2.
      DLx = 0.5/2.
      Lymin=-DLy*(1+1/2.)
      NLy=math.ceil(Ly/DLy)+2
      Lxmin=-DLx*(1+1/2.)
      NLx=math.ceil(Lx/DLx)+2
   else:		#calcula parche optimo para GNR_ac
      DLy = 0.5/2.
      DLx= 0.866025404/2.
      Lymin=-DLy*(2+1/2.)
      NLy=math.ceil(Ly/DLy)+6
      Lxmin=-DLx*(1+1/2.)
      NLx=math.ceil(Lx/DLx)+2
   p_ran=[NLx,NLy,DLx,DLy,Lxmin,Lymin]
   v_ran=[ [ '' for l in range(0,NLy)] for m in range(0,NLx)]
   for m in range(0,NLx):
#      xmin=Lxmin+m*DLx
#      xmax=Lxmin+(m+1)*DLx
#      x=(xmax+xmin)/2.
      for l in range(0,NLy):
#         ymin=Lymin+l*DLy
#         ymax=Lymin+(l+1)*DLy
#         y=(ymax+ymin)/2.
         RandE=0.0
         if noise_type=='Anderson':
           # Anderson´s noise added with W/V0 = noise_delta
           RandE=(1.-random.random())*noise_delta*V0
         elif noise_type=='Vacancies':
           #  A vacancy is the same as taking the site energy very very large
           if random.random()<=noise_delta:
              RandE=1000*V0 
         v_ran[m][l]=RandE
   return p_ran,v_ran


#   if noise_type=='Anderson':
#     RandE=(random.random()*noise_delta*V0
#     print(f'Anderson´s noise added with W/V0 = {noise_delta:<6g} and seed = {noise_seed:<10d}')
#   elif noise_type=='Vacancies':
#  A vacancy is the same as taking the site energy very very large
#     RandE=1000*V0 if random.random()<=noise_delta  else 0
#     print(f'{noise_delta:<6.2G} % of vacancies added with seed = {noise_seed:<10d}')
#   else:
#     RandE=0.0
#     print('No noise added.')





if __name__ == "__main__":
   print('Usted no debería ver esto.')
else:
   print('Funciones tomadas de GNR_functions_v0p6')
#   print('Ojo que Lx da distinto N que en GNR_plotter_v0p3.py !')
