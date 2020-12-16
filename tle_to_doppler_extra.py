# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 20:51:29 2020

@author: Mindaugas
"""

import numpy as np
from scipy.spatial.transform import Rotation as rot
import datetime

tle_file = "GNSS_2020_10_07.txt" #name of the TLE file
with open(tle_file,'r') as f:
    tle_data = f.read()

tle_data = tle_data.split("\n") #split the file into lines
for i in range(len(tle_data)): #split the lines into words and remove empty elements from the list
    if(i%3 != 0):    
        tle_data[i] = tle_data[i].split(" ")
        while("" in tle_data[i]):
            tle_data[i].remove("")

sat = [] #list of all objects in the TLE file

for i in range(0, len(tle_data), 3): #transfer of data from raw TLE into the object list
    sat.append({
        "name":tle_data[i],
        "year":2000+int(tle_data[i+1][3][:2]),
        "day":float(tle_data[i+1][3][1:]),
        "ballistic_coe":float(tle_data[i+1][4]),
        "inclination":float(tle_data[i+2][2]),
        "ascension":float(tle_data[i+2][3]),
        "eccentricity":float("0."+tle_data[i+2][4]),
        "argument_per":float(tle_data[i+2][5]),
        "mean_anomaly":float(tle_data[i+2][6]),
        "mean_motion":float(tle_data[i+2][7])
        })

#physical constant definitions
m_earth = 5.9722e24 #mass of central body
grav = 6.67430e-11 #gravitatinal constant
r_earth = 6371e3 #radius of central body
loc_lat = (54.+53./60.+50./3600.)/360*2*np.pi #positive north, latitude of observation point
loc_lon = (23.+53./60.+10./3600.)/360*2*np.pi #positive east, longitude of observation point
loc_rot = rot.from_euler('z', 0) #transformation of coordinates to sidereal ref
loc_rot = loc_rot * rot.from_euler('y', -loc_lat)
loc_v = loc_rot.apply([r_earth, 0, 0])
axis_dist_loc = (loc_v[0]**2 + loc_v[1]**2)**0.5 #location distance from earth axis
vel_loc = axis_dist_loc*2*np.pi/24/3600 #local velocity
vel_loc_v = loc_rot.apply([0, vel_loc, 0]) #local velocity rotated
sid_day = datetime.timedelta(hours = 23.9344696) #sidereal day length

phase_steps = 500
substeps = 2

#shortened attribute name declaration
name = "name"
year = "year"
eday = "day"
incd = "inclination"
ascd = "ascension"
ecce = "eccentricity"
arpd = "argument_per"
mand = "mean_anomaly"
mavr = "mean_angular_vel"
mmtd = "mean_motion"
smjl = "semimajor_axis"
smnl = "semiminor_axis"
apol = "apogee"
perl = "perigee"
prtq = "plane_rot"
pnmv = "plane_n"
sfdl = "semifocal_dist"
lper = "longitude_per"
trtq = "total_rot"
datt = "tle_date"
orbt = "full_orb_time"
tiar = "tle_anomaly_rad"
cvlv = "current_vel"
cpsv = "current_pos"
dopf = "doppler_freq"
arpq = "argument_per_rot"

progress_stps = 10

#physical calculations
for i in range(len(sat)):
    sat[i][mavr] = sat[i][mmtd]/24/3600*2*np.pi #mean angular spead in radians/s
    sat[i][smjl] = (grav*m_earth/(sat[i][mavr]**2))**(1./3) #semi major axis
    sat[i][apol] = sat[i][smjl]*(1+sat[i][ecce]) #apogee
    sat[i][perl] = sat[i][smjl]*(1-sat[i][ecce]) #perigee
    inc_rot = rot.from_euler('z', sat[i][ascd]*np.pi/180) #rotation to the ascending node
    inclination_vector = inc_rot.apply(np.array([1.0, 0.0, 0.0])) #rotation vector for plane inclination
    sat[i][prtq] = rot.from_rotvec(sat[i][incd]*np.pi/180 * inclination_vector) #plane inclication rotation
    sat[i][pnmv] = sat[i][prtq].apply(np.array([0.0, 0.0, 1.0])) #orbital plane normal vector
    sat[i][arpq] = rot.from_rotvec(sat[i][arpd]*np.pi/180 * sat[i][pnmv]) #rotation to argument of perigee
    sat[i][trtq] = sat[i][arpq] * sat[i][prtq] * inc_rot #total rotation of orbit
    sat[i][sfdl] = (sat[i][apol]-sat[i][perl])/2 #semi focal distance
    sat[i][smnl] = (sat[i][smjl]**2 - sat[i][sfdl]**2)**0.5 #semi minor axis
    sat[i][datt] = datetime.datetime(sat[i][year]-1, 12, 31, 0, 0, 0, 0, datetime.timezone.utc) + datetime.timedelta(days=sat[i][eday]) #date object of tle data
    sat[i][tiar] = sat[i][mand]*np.pi/180 #mean anomaly in radians
    sat[i][orbt] = 2*np.pi/sat[i][mavr] #orbital period in seconds
    sat[i][cpsv] = [np.array([sat[i][perl], 0, 0])] #first positional point in orbit (perigee), on the X axis
    sat[i][cvlv] = [np.array([0, (grav*m_earth*(2/sat[i][perl]-1/sat[i][smjl]))**0.5, 0])] #velocity on the first point in orbit (perigee), on the Y axis
    sat[i][dopf] = [] #list of doppler shifts
    time_step = sat[i][orbt]/phase_steps #time step of physical calculations
    curr_orb_time = time_step #time iterator
    j = 1 #array iterator
    
    while (curr_orb_time <= sat[i][orbt]): #while time iterator is less than orbital period
        grav_dir = -sat[i][cpsv][j-1]/np.linalg.norm(sat[i][cpsv][j-1]) #gravitational force direction vector
        grav_accel = grav_dir*grav*m_earth/np.linalg.norm(sat[i][cpsv][j-1])**2 #gravitational acceleration
        sat[i][cpsv].append(sat[i][cpsv][j-1]+sat[i][cvlv][j-1]*time_step+grav_accel/2*time_step**2) #new position
        sat[i][cvlv].append(sat[i][cvlv][j-1]+time_step*grav_accel) #new velocity
        for k in range(substeps): #iteration to make the new position and velocity more precise
            grav_dir = -sat[i][cpsv][j]/np.linalg.norm(sat[i][cpsv][j])
            grav_accel_n = grav_dir*grav*m_earth/np.linalg.norm(sat[i][cpsv][j])**2
            grav_accel_n = (grav_accel+grav_accel_n)/2
            sat[i][cpsv][j] = sat[i][cpsv][j-1]+sat[i][cvlv][j-1]*time_step+grav_accel_n/2*time_step**2
            sat[i][cvlv][j] = sat[i][cvlv][j-1]+time_step*grav_accel_n
        j += 1 #add to iterator
        curr_orb_time += time_step #add to time iterator
        
    for j in range(len(sat[i][cpsv])): #loop to transfer all orbital coordinates to global coordinate system
        sat[i][cpsv][j] = sat[i][trtq].apply(sat[i][cpsv][j]) #position
        sat[i][cvlv][j] = sat[i][trtq].apply(sat[i][cvlv][j]) #velocity
        
    if(i/len(sat)*progress_stps%1 < 1/len(sat)*progress_stps): #progress bar - 10 symbols until finished
        print("=", end = " ")
print("\n finished orbital calculations")
    
date_calc = datetime.datetime(2020, 10, 6, 0, 0, 0, 0, datetime.timezone(datetime.timedelta(hours = +3),"EEST")) #start time of doppler calculations
calc_period = datetime.timedelta(days = 0.5) #period for calculation
calc_step = datetime.timedelta(minutes = 20) #calculation timestep
base_freq = 1575.42e6 #base signal frequency
lightspeed = 3e8 #speed of light
curr_time = date_calc #time iterator
end_time = date_calc + calc_period #end time of calculation
plotx = [] #object to store x values of the plot (hours since start of calculation)

def GMST(CT): #Greenwich mean sidereal time calculation (in hours)
    CT_d = (CT - datetime.datetime(2000, 1, 1, 12, 0, 0, 0, datetime.timezone(datetime.timedelta(hours = +3))))/datetime.timedelta(days = 1)
    CT_mn = (CT_d - CT_d%1)+0.5
    CT_h = (CT_d-CT_mn)*24
    CT_c = CT_d / 36525
    GMST = 6.697374558 + 0.06570982441908*CT_mn + 1.00273790935*CT_h + 0.000026*CT_c**2
    return GMST%24

while(curr_time <= end_time): #iteration
    phi = GMST(curr_time)*2*np.pi/24-loc_lon #calculation of ground location rotation angle from sidereal coordinate origin
    curr_loc_v = rot.from_euler('z', phi).apply(loc_v) #applying ground rotation to location
    curr_vel_loc_v = rot.from_euler('z', phi).apply(vel_loc_v) #and velocity
    for i in range(len(sat)): #iteration through satellite list
        tle_dtime = curr_time - sat[i][datt] #date difference from tle data
        tle_dtime /= datetime.timedelta(days = 1) #time difference in days
        tle_dtime *= sat[i][mmtd] #time difference in orbits
        tle_dtime += sat[i][mand]/360 #mean anomaly compensation
        curr_phase = tle_dtime % 1 #truncating to fraction of one orbit
        j = int(curr_phase*phase_steps) #calculating the nearest index on the satellite position list
        curr_v = sat[i][cvlv][j] #velocity at location
        curr_p = sat[i][cpsv][j] #location
        vis_v = curr_p - curr_loc_v #direction vector from ground location to satellite
        vis_bool = np.dot(curr_loc_v, vis_v)>=0. #visibility check
        if(vis_bool == False): #if invisible
            sat[i][dopf].append(None) #append an empty object to plot data
        else: #if visible
            proj_dir = curr_p - curr_loc_v #direction vector from ground to satellite
            rel_vel = curr_v - curr_vel_loc_v #relative velocity
            proj_vel = -np.dot(rel_vel, proj_dir)/np.linalg.norm(proj_dir) #relative velocity projected to direction
            doppler = (proj_vel/lightspeed) * base_freq #doppler shift
            sat[i][dopf].append(doppler) #append shift to plot data
    plotx.append(-(date_calc-curr_time)/datetime.timedelta(hours = 1)) #append time to plot x coordinate
    curr_time += calc_step #iterate time

#ground truth data
h_p_px = 0.0083045
hz_p_px = 19.802
base_data_x = np.array([71, 121, 164, 203, 232, 265, 297, 327, 344, 1032, 1078, 1107, 1128, 1155, 1175, 1194, 1219])*h_p_px
base_data_y = np.array([319, 315, 296, 263, 227, 173, 108, 40, 0, 332, 293, 253, 216, 161, 113, 64, 0])*hz_p_px
    
import matplotlib.pyplot as plt

#create plot
plt.figure(figsize=(10, 5), dpi=149)
for i in range(len(sat)):
    plt.plot(plotx, sat[i][dopf], linewidth = 0.3, label = sat[i][name])
    
#plt.plot(base_data_x, base_data_y, linewidth = 0.3, label = "base_data") #ground truth plot
plt.grid()
plt.xticks(np.arange(0, 12, step=1))

# plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc = "lower center", mode = "expand", ncol=3)
plt.xlim(0, 12)
plt.xlabel("time, h (" + date_calc.strftime("%Y-%m-%d") + ")")
plt.ylabel("doppler, Hz")
plt.savefig("plot.png")

#extra calculation of specific time
curr_time = datetime.datetime(2020, 10, 6, 7, 28, 0, 0, datetime.timezone(datetime.timedelta(hours = +3),"EEST")) #specific time of checking
phi = GMST(curr_time)*2*np.pi/24-loc_lon #calculation of ground location rotation angle from sidereal coordinate origin
curr_loc_v = rot.from_euler('z', phi).apply(loc_v) #applying ground rotation to location
curr_vel_loc_v = rot.from_euler('z', phi).apply(vel_loc_v) #and velocity
for i in range(len(sat)): #iteration through satellite list
    tle_dtime = curr_time - sat[i][datt] #date difference from tle data
    tle_dtime /= datetime.timedelta(days = 1) #time difference in days
    tle_dtime *= sat[i][mmtd] #time difference in orbits
    tle_dtime += sat[i][mand]/360 #mean anomaly compensation
    curr_phase = tle_dtime % 1 #truncating to fraction of one orbit
    j = int(curr_phase*phase_steps) #calculating the nearest index on the satellite position list
    curr_v = sat[i][cvlv][j] #velocity at location
    curr_p = sat[i][cpsv][j] #location
    vis_v = curr_p - curr_loc_v #direction vector from ground location to satellite
    vis_bool = np.dot(curr_loc_v, vis_v)/np.linalg.norm(curr_loc_v)/np.linalg.norm(vis_v)>=0.15 #visibility check
    if(vis_bool == True): #if visible
        proj_dir = curr_p - curr_loc_v #direction vector from ground to satellite
        rel_vel = curr_v - curr_vel_loc_v #relative velocity
        proj_vel = -np.dot(rel_vel, proj_dir)/np.linalg.norm(proj_dir) #relative velocity projected to direction
        doppler = str((proj_vel/lightspeed) * base_freq) #doppler shift
        output = "Satellite name: " + sat[i][name] + "  Doppler shift: " + doppler
        print(output)