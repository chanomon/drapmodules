#!/usr/bin/env python

#
# NOAAa DRAP model implementation
# 
# Copyright (C) 2021 Author:Elizandro Huipe Domratcheva (hdomeli@gmail.com)
#                      
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import math
import scipy.stats as st
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('TkAgg')  # Usa el backend TkAgg (interactivo)
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ephem
from datetime import datetime
#import turtle
import geopy.distance
import matplotlib.cm as cm
import pandas as pd
from tqdm import tqdm
import glob
import time as TIME
mpl.rcParams["savefig.directory"] = './'
nearside = ccrs.NearsidePerspective(central_latitude=0)#sun_lat)
carree = ccrs.PlateCarree(central_longitude=180.0)
#global parameters
RT = geopy.distance.EARTH_RADIUS #[km]  earth radius
proj = ccrs.NearsidePerspective()

def subsolarpoint(datetime):#subsolar point from datetime  object
    greenwich = ephem.Observer()
    greenwich.lat = "0"
    greenwich.lon = "0"
    greenwich.date = datetime
    sun = ephem.Sun(greenwich)
    sun.compute(greenwich.date)
    sun_lon = math.degrees(sun.ra - greenwich.sidereal_time() )
    if sun_lon < -180.0 :
      sun_lon = 360.0 + sun_lon 
    elif sun_lon > 180.0 :
      sun_lon = sun_lon - 360.0
    sun_lat = math.degrees(sun.dec)
    return(sun_lon,sun_lat)

def calc_param(flux,time):        #calculate parameters of solar afectation 
    freqs = np.arange(5,40)#[5,10,15,25,30,35]    #(AFradius for frequencies and subsolarpoints) 
    lon_sun,lat_sun = subsolarpoint(time)
    fdic = {}
    HAF = (10*math.log10(flux)+65)
    for f in freqs:
        #argument = (f/HAF)**(4/3)
        #chirad = np.arccos((f/HAF)**(4/3)) ##zenith angle in radians
        #if not math.isnan(chirad):
        #    r_sun = chirad*RT                                       ##radius of afectation
        #else:
        #    r_sun = 0.0
        ratio = (f / HAF) ** (4 / 3)
        
        try:
            chirad = np.arccos(ratio)  # Intenta calcular el arco coseno
            r_sun = chirad * RT  # Calcula el radio de afectación
        except ValueError:  # Captura el error si ratio está fuera del rango [-1, 1]
            r_sun = 0.0  # Asigna 0.0 si no se puede calcular chirad
        fdic[str(f)]=[r_sun]
    for f in range(len(freqs)):
        if f == len(freqs)-1:
            fdic[str(freqs[f])].append(0.0)
        else:
            fdic[str(freqs[f])].append(fdic[str(freqs[f+1])][0])
            
    return lon_sun,lat_sun,fdic

def returnpointsAF3(lonlist, latlist, subsolarlon,subsolarlat,RA,dictoAF): #return affected frequencies 
                                                                           #for points given
    t = []
    AFlon = []
    AFlat = []
    subsolarcoords = (subsolarlat,subsolarlon)
    for lt,ln in zip(latlist,lonlist):
        MaxdistAF = dictoAF['5'][0]   ######Maximum distance of Affected Freqs
        dist = geopy.distance.distance((lt,ln),subsolarcoords).km
        AFlon.append(ln)
        AFlat.append(lt)
        if dist > MaxdistAF:
            t.append(0)
        else:
            for f in dictoAF:
                if dictoAF[f][0] != 0.:
                    if dist <= dictoAF[f][0] and dist > dictoAF[f][1]:
                        t.append(int(f))
                
    return AFlon,AFlat,t

def Drapmap3(flux,time):#ploting DRAP with NOAA points
    lons = np.arange(-178,180,2)
    lats = np.arange(-89,91,2)
    
    lon = []
    lat=[]
    for n in lons:
        for t in lats:
            lon.append(n)
            lat.append(t)
    params = calc_param(flux,time)
    RA =  params[2]['5'][0]  #[km]
    subslon = params[0]
    subslat = params[1]
    dictoAF = params[2]
    AFlon,AFlat,t=returnpointsAF3(lon,lat,subslon,subslat,RA,dictoAF )
    plt.figure(figsize=(6, 3))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0))
    ax.set_extent([-180, 180, -90, 90])#, crs=ccrs.LambertCylindrical())
    plt.title('DRAP  '+ 'flux='+format(flux, ".2e")+'  '+time.strftime("%Y-%m-%d %H:%M"))
    ax.coastlines(resolution='110m')
    AF = ax.scatter(AFlon,AFlat,c=t, cmap="jet",vmin=0, vmax=35,alpha=0.3)#,transform=carree)#,transform=ccrs.PlateCarree())
    ax.scatter([subslon],[subslat],marker='*',color='orange',label='subsolar point')
    plt.colorbar(AF, label='HAF by 1dB Absorption [MHz]')
    plt.rcParams['figure.dpi'] = 200
    #plt.legend()
    plt.show()


def subsolarpointdistribution(fluxes,datetimes): #input: a list of xrs fluxes in the long band and a list of their respective t_max in datetime object.
    dx = 5                                       #output: a finite distribution of subsolarpoints over an earths map.
    dy = 5
    size = len(fluxes)
    fraction = 1/size
    print(size)
    lon = range(180,-180,-dx)#[str(i) for i in range(+180,-180,-1)]
    lat = range(90,-90,-dy)#[str(j) for j in range(+90,-90,-1)]
    geomatrix = pd.DataFrame(0, index=lat, columns=lon)
    #flag = 0
    #flag2= 0
    for f,d in zip(fluxes,datetimes):
        lo,la = subsolarpoint(d)
        for index in geomatrix:
            #print(index)
            if la<index and la>index-dy:
                y = index
                #flag += 1
                #print(index)
                break
            #else:
            #    pass
        for column in geomatrix:
            #print(column)
            if lo<column and lo>column-dx:
                x = column
                #flag2+= 1
                #print(column)
                break
            #else:
            #    pass
        old_value = geomatrix.loc[y,x]
        #print(x,y,f*1e8)
        geomatrix.at[y,x]=old_value+fraction
        #print(geomatrix.at[y,x])
    #print(flag,flag2)
    #plt.figure(figsize=(6, 3))
    #ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0))
    #ax.set_extent([-180, 180, -90, 90])#, crs=ccrs.LambertCylindrical())
    #plt.title('Solar Flare Energy distribution')
    #ax.coastlines(resolution='110m')
    #extent = [-180, 180, -90, 90]
    #plt.imshow(heatmap.T, extent=extent, origin='lower')
    #plt.imshow(geomatrix,extent=extent, cmap ="jet")
    #plt.hist2d(geomatrix.columns,geomatrix.index,geomatrix, cmap='jet')
    #plt.xticks(range(len(geomatrix.columns)), geomatrix.columns)
    #plt.yticks(range(len(geomatrix.index)), geomatrix.index)
    #plt.rcParams['figure.dpi'] = 200
    #plt.colorbar(label='Xray flux (x1e8) in 1975-2017 ')
    #plt.show()

    plt.figure(figsize=(8, 4))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0.0))
    ax.set_extent([-180, 180, -90, 90])#, crs=ccrs.LambertCylindrical())
    ax.coastlines(resolution='110m')
    extent = [-180, 180, -90, 90]
    ax.set_xticks(np.arange(-180, 181, 60), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(-90, 91, 30), crs=ccrs.PlateCarree())
    ax.set_xlabel('Lon')
    ax.set_ylabel('Lat')
    plt.imshow(geomatrix,extent=extent, cmap ="jet")
    plt.rcParams['figure.dpi'] = 200
    plt.colorbar(shrink=0.6,label='Percentage of subsolar points')#fraction=0.0235
    plt.show()


#if __name__ == "__main__":
##### Definición de parámatros de fulguración
####como ejemplo se pone el evento X1 del 28 de Octubre
flux = 1e-4   ##flujo de rayos X [W/m^2 s] ver como referencia la escala de la NOAA de fulguraciones
time = datetime(2021, 10 , 28, hour=15, minute=35)	## Hora y fecha en tiempo universal del evento en formato YYYY,MM,DD hour=hh, minute=mm
							## no usar cifras con un cero antes, por ejemplo si ocurrió
							## un evento en la hora "01:08" poner "hour=1, minute=8"

Drapmap3(flux,time)

