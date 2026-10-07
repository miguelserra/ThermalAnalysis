#############################################################
from math import sin, cos, pi, log
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
#############################################################


def standard_fire_curve(double t):
    #t in minutes
    return 20+345*log(8*t+1,10)

    
def density_concrete(double T):
    #T in centigrade

    cdef double ro = 2400
    
    if T > 115 and T <= 200:
        ro *= 1 - 0.02 * (T - 115) / 85
    elif T > 200 and T <= 400:
        ro *= 0.98 - 0.03 * (T - 200) / 200
    elif T > 400 and T <= 1200:
        ro *= 0.95 - 0.07 * (T - 400) / 800
    else:
        ro *= 1
    
    return ro
        

def specific_heat_concrete(double T, double moisture_percent):
    #T in centigrade
    #moisture_percent in %
    
    cdef double Cmoist= 373.33 * moisture_percent + 903.33

    if T > 100 and T < 115:
        return 900 + (T - 100)
    elif T >= 100 and T <= 115:
        return Cmoist
    elif T >= 100 and T <= 115:
        return (1000 - Cmoist) / 85 * T + 40 / 17 * (Cmoist - 575)
    elif T > 200 and T <= 400:
        return 1000 + (T - 200) / 2
    elif T > 400 and T <= 1200:
        return 1100
    else:
        return 900

    
def thermal_conductivity_concrete(double T):
    #T in centigrade
    return 1.36 - 0.136 * (T / 100) + 0.0057 * (T / 100)**2
    
    
def heat_flux_convection(double Tgas, double Tsurface, double alpha_c = 25):
    #Tgas and Tsurface in centigrade
    #alpha_C [W/m2K] - coefficient of heat transfer by convection
    return alpha_c * (Tgas - Tsurface)


def heat_flux_radiation(double Tgas,
                        double Tsurface,
                        double phi      = 1.0 ,
                        double epson_m  = 0.7 ,
                        double epson_f  = 1.0 , 
                        double sigma    = 5.67e-8 ):

    """
    Tgas and Tsurface in centigrade
    phi     = 1         #configuration factor
    epson_m = 0.7       #surface emissivity of the member - 0.7 for CONCRETE!!!
    epson_f = 1.0       #emissivity of the fire
    sigma   = 5.67e-8   #[W/m2K4] - Stephan Boltzmann constant
    """
    
    return phi*epson_m*epson_f*sigma*((Tgas+273)**4-(Tsurface+273)**4)

    
    
def myprint(title,array):
    result = "%s\n" % title
    for line in array:
        for val in line:
            result += "%.2f\t" % val
        result += "\n"
    print(result)


# class TimeMe:
    # def __init__(self,label):
        # self.label = label
        # self.sumTime = 0
    # def START(self):
        # self.time = time.time()
    # def STOP(self): 
        # self.sumTime += time.time()-self.time
    # def REPORT(self):
        # print("%s: %.2es" % (self.label,self.sumTime))




def RUN(double H, double B, double t_final, int divX = 100, int divY = 100, int material=1):  
    """
    > Materials:
        Concrete    = 1
        Steel       = 2
        Composite   = 3

    """

    cdef double max_dt, time_print_every

    max_dt = .075            #min

    time_print_every = 5    #analysis min

    #############################################################   
    """ SETUP """
    cdef int i,j,k,u,v

    cdef double dx = B*1e-3/divX
    cdef double dy = H*1e-3/divY

    cdef int vx, vy, dir, maxT
    cdef double Tgas, Tij, T0, T1, NetCHF, NetRHF, sumConduct, sumHF, vol_specific_heat
    
    cdef double[:,:]    Temperatures
    cdef double[:,:]    coef_x     
    cdef double[:,:]    coef_y     
    cdef double[:,:,:]  Conductance
    cdef double[:,:,:]  HeatFlows  
    
    
    
    Directions  = (0,1,2,3)
    vec_x       = (0,1,0,-1)
    vec_y       = (-1,0,1,0)




    cdef double moisture_percent
    if material == 1:
        
        #constants
        moisture_percent = 1.5
        
        #functions
        specific_heat        = lambda T: specific_heat_concrete(T, moisture_percent)
        density              = density_concrete
        thermal_conductivity = thermal_conductivity_concrete
        
        
    #############################################################



    #temperature gets 1 extra row/column on each side to store gas temperatures    
    Temperatures = np.array(  [ [  float(20) for i in range(0,divX+2) ] for j in range(0,divY+2) ]                      , dtype = np.float64)
    coef_x       = np.array(  [ [  float( 0) for i in range(0,divX+2) ] for j in range(0,divY+2) ]                      , dtype = np.float64)
    coef_y       = np.array(  [ [  float( 0) for i in range(0,divX+2) ] for j in range(0,divY+2) ]                      , dtype = np.float64)
    Conductance  = np.array([ [ [  float( 0) for i in range(0,divX+2) ] for j in range(0,divY+2) ] for k in Directions ], dtype = np.float64)
    HeatFlows    = np.array([ [ [  float( 0) for i in range(0,divX+2) ] for j in range(0,divY+2) ] for k in Directions ], dtype = np.float64)



    cdef double dt = 1e-8
    cdef double curr_time = dt
    cdef int count_time_print = 1

    while curr_time <= t_final:
        if curr_time > time_print_every*count_time_print:
            print("### TIME = %.3f min\tCurrent dt = %.3f" % (curr_time,dt))
            count_time_print += 1
            
        #assuming here temperature is the same on all four sides !!!!!!!!!!!!!!!!
        Tgas = standard_fire_curve(curr_time)

        for i in range(0,divX+2):
            for j in range(0,divY+2):
                if (i in [0, divX+1] ) or ( j in [0, divY+1]  ):
                    Temperatures[i,j] = Tgas


        #this is constant throughout the step
        for i in range(1,divX+1):
            for j in range(1,divY+1):
                Tij = Temperatures[i,j]
                coef_x[i,j] = dx/(2*thermal_conductivity(Tij))
                coef_y[i,j] = dy/(2*thermal_conductivity(Tij))

        
        for dir in Directions:
            
            vx, vy = vec_x[dir], vec_y[dir]

            
            for i in range(1,divX+1):
            
                for j in range(1,divY+1):

                    u,v = i+vx,j+vy
                
                    T0 = Temperatures[i,j]
                    T1 = Temperatures[u,v]
                      
                    if dir % 2 == 0: #assess direction of analysis
                        Conductance[dir,i,j] = dy/(coef_x[i,j] + coef_x[u,v])
                    else:
                        Conductance[dir,i,j] = dx/(coef_y[i,j] + coef_y[u,v])

                        
                    if ( i in [1, divX] ) or ( j in [1, divY] ): 
                        #heat flows from boundary
                        NetCHF = heat_flux_convection(T1,T0)
                        NetRHF =  heat_flux_radiation(T1,T0)
                        
                        if dir % 2 == 0:
                            HeatFlows[dir,i,j] =  dy*(NetCHF + NetRHF)
                        else:
                            HeatFlows[dir,i,j] =  dx*(NetCHF + NetRHF)
                            
                                
                    else:
                        #heat flows from internal elements
                        HeatFlows[dir,i,j] = Conductance[dir,i,j] * (T1-T0)

                        
        lim_dt = 1e3
        for j in range(1,divY+1):
            
            for i in range(1,divY+1):
            
                T0 = Temperatures[i,j]
                
                sumConduct = 0
                sumHF      = 0
                
                for dir in Directions:
                    sumConduct += Conductance[dir,i,j]
                    sumHF      +=   HeatFlows[dir,i,j]


                vol_specific_heat = specific_heat(T0) * density(T0)
                
                lim_dt = min(vol_specific_heat*dx*dy/sumConduct/60,lim_dt)
                
                Temperatures[i][j] +=  (dt*60)/(dx*dy*vol_specific_heat)*sumHF


        dt = min(lim_dt, max_dt)
        curr_time += dt
        

    #adicionar scatter com a posicao das armaduras


    #eliminate temperature lines related to gas temperatures
    Temperatures   = np.array([ [ val for val in line[1:-1] ] for line in Temperatures[1:-1]])
    
    x = np.linspace(-B/2,B/2,divX)
    y = np.linspace(-H/2,H/2,divY)
    #x_grid = np.linspace(0,B,300)
    #y_grid = np.linspace(0,H,300)
    #Z     = griddata((x, y), Temperatures, (x_grid[None,:], y_grid[:,None]), method='cubic')
    X, Y  = np.meshgrid(x,y)
    
    
    maxT = int(np.amax(Temperatures)/100)*100+100
    levels   = [ i for i in range(0,maxT+1, 5 ) ]
    isolines = [ i for i in range(0,maxT+1,100) ]
    CS1 = plt.contour(X, Y, Temperatures, levels=isolines, colors='k')
    CS2 = plt.contourf(X, Y, Temperatures, levels)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.tick_params(direction='out')
    plt.colorbar(ticks = isolines)
    plt.clabel(CS1, inline=1, fontsize=7.5, fmt='%i', manual = True)
    
    plt.show()
