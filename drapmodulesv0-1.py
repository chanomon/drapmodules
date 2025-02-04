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
        chirad = np.arccos((f/HAF)**(4/3)) ##zenith angle in radians
        if not math.isnan(chirad):
            r_sun = chirad*RT                                       ##radius of afectation
        else:
            r_sun = 0.0
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
#############################################
##### Definición de parámatros de fulguración
#############################################
####como ejemplo se pone el evento X1 del 28 de Octubre
flux = 1e-4   ##flujo de rayos X [W/m^2 s] ver como referencia la escala de la NOAA de fulguraciones
time = datetime(2021, 10 , 28, hour=15, minute=35)	## Hora y fecha en tiempo universal del evento en formato YYYY,MM,DD hour=hh, minute=mm
							## no usar cifras con un cero antes, por ejemplo si ocurrió
							## un evento en la hora "01:08" poner "hour=1, minute=8"

Drapmap3(flux,time)

