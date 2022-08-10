import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

##########################################################################
#Arbitrary number of cells in a triangular lattice stimulated by a square 
#wave. Each cells has a randomly choosen hmax value
#The number of peaks past threshold are counted
#The number of edges past between nearest-neighbors are counted. The peaks 
#between nearest-neighbors needs to be at least five time unit apart.

##########################################################################
#Number of trials to average over for a given g value
trials = 100

#Total number of cells in lattice
row_size = 25 
number_of_cells = row_size**2
k = 1 #k does nothing, args needs to be a tuple in odeint


N = 100
theta = 0
hmin_value = -.25 #Monostable low regime
f = .01

# a needs to be defined within the function
b = f + 1/(N*(theta + 1))
c = 1/(3*(N**2)*(theta + 1))
epsilon = .1
d = epsilon
e = epsilon

g = np.array([.0001, .00025, .0005, .001, .0025, .005, .01, .025, .05, .1, .25, .5, 1])
gi = 0

I = 50
t = np.linspace(0,3*I,(10*(3*I) + 1))


ave_number_of_edges_list = np.array([])
std_number_of_edges_list = np.array([])

ave_number_of_peaks_list = np.array([])
std_number_of_peaks_list = np.array([])



#Starting loop for a given g value
while gi < np.size(g):

    #Arrays to be filled for a given g value
    number_of_edges_list = np.array([])
    number_of_peaks_list = np.array([])
    
    i = 0
    while i < trials:

        #Each cell has a random normally distributed hmax
        #The boundary between mono stable low and excitable it about -0.128
        #The boundary for the excitable regime is -0.068
        mu_hmax = -.127 #.05 std from boundary between mono low and excitable regimes, this ensures ~52% of cells are placed in excitable regime and ~48% stay in mono low
        sigma_hmax = .0198
        
        #hmax = np.random.normal(mu_hmax, sigma_hmax, number_of_cells)
        hmax = np.random.normal(mu_hmax, sigma_hmax, number_of_cells)
        #h is bounded below by -1/3 and the excitable-oscillatory regime is around -.0068
        ii = 0
        while ii < np.size(hmax):
            if hmax[ii] < hmin_value:
                hmax[ii] = hmin_value
            elif hmax[ii] > -.068:
                hmax[ii] = -.068
                
            ii += 1
        ii = 0
        hmin = np.array([])
        while ii < np.size(hmax):
            hmin= np.append(hmin, hmin_value)
            ii += 1
        ii = 0
        
        
        #Creating array of initial conditions. Note: they are all 0
        def IC(i):
            IC_list = np.array([])
            j = 0
            while j < i:
                IC_list = np.append(IC_list, 0)
                j += 1
            return IC_list
        
        
        
        def MM(state,t, i, k):
            ODEs = np.array([])

            if (t // I) % 2 != 0:
                h = hmax     
            else:
                h = hmin 
            
            j = 1
            while j < i + 1:
                x = state[2*(j-1)]
                y = state[2*(j-1) + 1] 
                
                
                
                #x_r = state[2*(j-1) + 2]
                #x_l = state[2*(j-1) - 2]
                
                #For odd rows
                #x_otr = state[2*((j-(row_size))-1)]
                #x_otl = state[2*((j-(row_size+1))-1)]
                #x_obr = state[2*((j+(row_size))-1)]
                #x_obl = state[2*((j+(row_size-1))-1)]
                
                #For even rows
                #x_etr = state[2*((j-(row_size-1))-1)]
                #x_etl = state[2*((j-(row_size))-1)]
                #x_ebr = state[2*((j+(row_size+1))-1)]
                #x_ebl = state[2*((j+(row_size))-1)]
                
                
                
                #1st row
                #Note: First row is always odd
                if (j - 1) // row_size == 0:
                    #1st column
                    if (j - 1) % row_size == 0:
                        x_r = state[2*(j-1) + 2]
                        x_obr = state[2*((j+(row_size))-1)]
                        a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_r + g[gi]*x_obr - 2*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                    
                    
                    #middle columns
                    elif (j - 1) % row_size != 0 and j % row_size != 0:
                        x_r = state[2*(j-1) + 2]
                        x_l = state[2*(j-1) - 2]
                        x_obr = state[2*((j+(row_size))-1)]
                        x_obl = state[2*((j+(row_size-1))-1)]
                        a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_r + g[gi]*x_obl + g[gi]*x_obr - 4*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                    
                    
                    #last column
                    elif j % row_size == 0:
                        x_l = state[2*(j-1) - 2]
                        x_obr = state[2*((j+(row_size))-1)]
                        x_obl = state[2*((j+(row_size-1))-1)]
                        a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_obl + g[gi]*x_obr - 3*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                    
                    
                    
                #middle rows
                elif (j - 1) // row_size != 0 and (j - 1) // row_size != (row_size - 1):
                    
                    #even rows
                    if ((j - 1) // row_size) % 2 != 0:
                        
                        #1st column
                        if (j - 1) % row_size == 0:
                            x_r = state[2*(j-1) + 2]
                            x_etr = state[2*((j-(row_size-1))-1)]
                            x_etl = state[2*((j-(row_size))-1)]
                            x_ebr = state[2*((j+(row_size+1))-1)]
                            x_ebl = state[2*((j+(row_size))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_r + g[gi]*x_etl + g[gi]*x_etr + g[gi]*x_ebl + g[gi]*x_ebr - 5*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt)    
                            
                            
                        #middle columns
                        elif (j - 1) % row_size != 0 and j % row_size != 0:
                            x_r = state[2*(j-1) + 2]
                            x_l = state[2*(j-1) - 2]
                            x_etr = state[2*((j-(row_size-1))-1)]
                            x_etl = state[2*((j-(row_size))-1)]
                            x_ebr = state[2*((j+(row_size+1))-1)]
                            x_ebl = state[2*((j+(row_size))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_r + g[gi]*x_etl + g[gi]*x_etr + g[gi]*x_ebl + g[gi]*x_ebr - 6*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt)
                            
                    
                        #last column
                        elif j % row_size == 0:
                            x_l = state[2*(j-1) - 2]
                            x_etl = state[2*((j-(row_size))-1)]
                            x_ebl = state[2*((j+(row_size))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_etl + g[gi]*x_ebl - 3*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt)
                            
                            
                            
        
                    #odd rows
                    else: #((j - 1) // row_size) % 2 == 0
                    
                        #1st column
                        if (j - 1) % row_size == 0:
                            x_r = state[2*(j-1) + 2]
                            x_otr = state[2*((j-(row_size))-1)]
                            x_obr = state[2*((j+(row_size))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_r + g[gi]*x_otr + g[gi]*x_obr - 3*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt)
                            
                            
                            
                    
                        #middle columns
                        elif (j - 1) % row_size != 0 and j % row_size != 0:
                            x_r = state[2*(j-1) + 2]
                            x_l = state[2*(j-1) - 2]
                            x_otr = state[2*((j-(row_size))-1)]
                            x_otl = state[2*((j-(row_size+1))-1)]
                            x_obr = state[2*((j+(row_size))-1)]
                            x_obl = state[2*((j+(row_size-1))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_r + g[gi]*x_otl + g[gi]*x_otr + g[gi]*x_obl + g[gi]*x_obr - 6*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt)
                            
                            
                            
                    
                        #last column
                        elif j % row_size == 0:
                            x_l = state[2*(j-1) - 2]
                            x_otr = state[2*((j-(row_size))-1)]
                            x_otl = state[2*((j-(row_size+1))-1)]
                            x_obr = state[2*((j+(row_size))-1)]
                            x_obl = state[2*((j+(row_size-1))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_otl + g[gi]*x_otr + g[gi]*x_obl + g[gi]*x_obr - 5*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt)
                            
                
                
                
                #last row
                elif (j - 1) // row_size == (row_size - 1):
                    
                    #even rows
                    if ((j - 1) // row_size) % 2 != 0:
                        
                        #1st column
                        if (j - 1) % row_size == 0:
                            x_r = state[2*(j-1) + 2]
                            x_etr = state[2*((j-(row_size-1))-1)]
                            x_etl = state[2*((j-(row_size))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_r + g[gi]*x_etl + g[gi]*x_etr - 3*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt)
                            
        
                    
                        #middle columns
                        elif (j - 1) % row_size != 0 and j % row_size != 0:
                            x_r = state[2*(j-1) + 2]
                            x_l = state[2*(j-1) - 2]
                            x_etr = state[2*((j-(row_size-1))-1)]
                            x_etl = state[2*((j-(row_size))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_r + g[gi]*x_etl + g[gi]*x_etr - 4*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt)
                            
                            
                    
                        #last column
                        elif j % row_size == 0:
                            x_l = state[2*(j-1) - 2]
                            x_etl = state[2*((j-(row_size))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_etl - 2*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt)   
                            
                            
                            
                        
                    #odd rows
                    else: #((j - 1) // row_size) % 2 == 0
                    
                        #1st column
                        if (j - 1) % row_size == 0:
                            x_r = state[2*(j-1) + 2]
                            x_otr = state[2*((j-(row_size))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_r + g[gi]*x_otr - 2*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt) 
                            
                    
                        #middle columns
                        elif (j - 1) % row_size != 0 and j % row_size != 0:
                            x_r = state[2*(j-1) + 2]
                            x_l = state[2*(j-1) - 2]
                            x_otr = state[2*((j-(row_size))-1)]
                            x_otl = state[2*((j-(row_size+1))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_r + g[gi]*x_otl + g[gi]*x_otr - 4*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt) 
                            
                            
                            
                    
                        #last column
                        elif j % row_size == 0:
                            x_l = state[2*(j-1) - 2]
                            x_otr = state[2*((j-(row_size))-1)]
                            x_otl = state[2*((j-(row_size+1))-1)]
                            a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                            dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g[gi]*x_l + g[gi]*x_otl + g[gi]*x_otr - 3*g[gi]*x
                            ODEs = np.append(ODEs, dx1dt)
                            dy1dt = d*x - e*y
                            ODEs = np.append(ODEs, dy1dt) 
                
                
                
                
                j += 1
                
            
                
            return ODEs

        #Solving ODEs
        IC = IC(2*number_of_cells)
        X = odeint(MM, IC, t, args = (number_of_cells,k)).T
        
        
        #Finding peak times and recording them into a list
        
        #Recording time series of of X molecules into a list called x_list
        jj = 0
        #x_list = np.array([X[0], X[2], X[4], X[6], X[8], X[10], X[12], X[14], X[16]])
        x_list = np.array([t]) #t is just a place holder so more can be appended
        while jj < 2*number_of_cells:
            x_list = np.append(x_list, [X[jj]], axis = 0)
            jj += 2
        jj = 0
        #print(np.size(x_list[0]))
        
        #Recording the time of each peak (past threshold) in a time series of x into a list called tmax_list.
        #If a cell peaks past threshold, an 1 is places in a number_of_peaks
        #If the time series does not pass threshold a 0 is used as a place holder
        jj = 1
        tmax_list = np.array([])
        number_of_peaks = np.array([])
        while jj < number_of_cells + 1:
            xmax = np.max(x_list[jj])
            #print(xmax)
            if xmax > 250:
                number_of_peaks = np.append(number_of_peaks, 1)
                tmax = np.where(x_list[jj] == xmax)
                #print(tmax)
                tmax_list = np.append(tmax_list, t[tmax])   
            else:
                number_of_peaks = np.append(number_of_peaks, 0)
                tmax_list = np.append(tmax_list, 0)
            jj += 1
        
        #The number of peaks for each trial is added to number_of_peaks_list
        number_of_peaks_list = np.append(number_of_peaks_list, np.sum(number_of_peaks))
        
        #Counting the number of edges between nearest-neighbors
        #j is the cell's position in the lattice. 0 starts at top left
        #j is also the cell's index in tmax_list
        #Cells are only compared to their neighbors down and to the right
        number_of_edges = np.array([])
        j = 0
        lower_time = .5
        upper_time = 2.5
        while j < np.size(tmax_list) - 1: #there is a (-1) because nothing is done for the last cell
        

            #Not last row
            if j // row_size != (row_size - 1):
                #Odd rows
                if (j // row_size) % 2 == 0:
                    
                    #1st column
                    if j % row_size == 0:
                        #right
                        if tmax_list[j] != 0 and tmax_list[j+1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+1]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+1]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                            
                        #bottom right
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                   
                    #Middle columns
                    elif j % row_size != 0 and (j + 1) % row_size != 0:
                        #right
                        if tmax_list[j] != 0 and tmax_list[j+1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+1]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+1]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                        #bottom left
                        if tmax_list[j] != 0 and tmax_list[j+row_size-1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size-1]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size-1]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                        #bottom right
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
    
                    #Last column
                    elif (j + 1) % row_size == 0:
                        #bottom left
                        if tmax_list[j] != 0 and tmax_list[j+row_size-1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size-1]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size-1]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                          
                        #bottom right
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)

                    
                    
                #Even rows 
                else: #((j // row_size) % 2 != 0
                    #Not last column
                    if (j + 1) % row_size != 0:
                        #right
                        if tmax_list[j] != 0 and tmax_list[j+1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+1]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+1]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                        #bottom left
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                        #bottom right
                        if tmax_list[j] != 0 and tmax_list[j+row_size+1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size+1]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size+1]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                    
                    #Last column
                    else: #(j + 1) % row_size == 0
                        #bottom left
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                
                
        
            #Last row
            else: #j // row_size == (row_size - 1)
            #The last row only needs to check to the right regardless of the column
                #right
                if tmax_list[j] != 0 and tmax_list[j+1] != 0:
                    if np.absolute(tmax_list[j] - tmax_list[j+1]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+1]) <= upper_time:
                        number_of_edges = np.append(number_of_edges, 1)
                    else:
                        number_of_edges = np.append(number_of_edges, 0)
                        
                else:
                    number_of_edges = np.append(number_of_edges, 0)
        
            j += 1
            
        #The number of edges for each trial is added to number_of_edges_list
        number_of_edges_list = np.append(number_of_edges_list, np.sum(number_of_edges))    
        
        

        
        i += 1
        
    
    #Averaging the number of peaks for all trials for a given g value
    #print(number_of_peaks_list)
    ave_number_of_peaks = np.average(number_of_peaks_list)
    #print(ave_number_of_peaks)
    std_number_of_peaks = np.sqrt(np.var(number_of_peaks_list, ddof = 1))
    #print(std_number_of_peaks)
    
    ave_number_of_peaks_list = np.append(ave_number_of_peaks_list, ave_number_of_peaks)
    std_number_of_peaks_list = np.append(std_number_of_peaks_list, std_number_of_peaks)

    
    #Averaging the number of edges for all trials for a given g value
    #print(number_of_edges_list)
    ave_number_of_edges = np.average(number_of_edges_list)
    #print(ave_number_of_edges)
    std_number_of_edges = np.sqrt(np.var(number_of_edges_list, ddof = 1))
    #print(std_number_of_edges)
    
    ave_number_of_edges_list = np.append(ave_number_of_edges_list, ave_number_of_edges)
    std_number_of_edges_list = np.append(std_number_of_edges_list, std_number_of_edges)
    
    

    
    #print(gi)
    gi += 1
    


np.savez("625_Cells_Pedge_FiringFrac_Results", g = g, ave_number_of_edges_list = ave_number_of_edges_list, std_number_of_edges_list = std_number_of_edges_list, ave_number_of_peaks_list = ave_number_of_peaks_list, std_number_of_peaks_list = std_number_of_peaks_list)
