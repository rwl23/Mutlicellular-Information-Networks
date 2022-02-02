import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

##########################################################################
#Arbitrary number of cells in a triangular lattice stimulated by a square 
#wave. Each cells has a randomly choosen hmax value

##########################################################################
#Number of trials
trials = 25

#Total number of cells in lattice
row_size = 18 
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

g = .025 


I_off = 40
I_on = 40
I_total = I_off + I_on

number_of_edges_list1 = np.array([])
number_of_edges_list2 = np.array([]) 
number_of_edges_list3 = np.array([]) 
number_of_edges_list4 = np.array([]) 
number_of_edges_list5 = np.array([]) 
number_of_edges_list6 = np.array([]) 
number_of_edges_list7 = np.array([]) 
number_of_edges_list8 = np.array([])  

number_of_edges_lists = [number_of_edges_list1, number_of_edges_list2, number_of_edges_list3, number_of_edges_list4, number_of_edges_list5, number_of_edges_list6, number_of_edges_list7, number_of_edges_list8]

number_of_peaks_list1 = np.array([])
number_of_peaks_list2 = np.array([])
number_of_peaks_list3 = np.array([])
number_of_peaks_list4 = np.array([])
number_of_peaks_list5 = np.array([])
number_of_peaks_list6 = np.array([])
number_of_peaks_list7 = np.array([])
number_of_peaks_list8 = np.array([])

number_of_peaks_lists = [number_of_peaks_list1, number_of_peaks_list2, number_of_peaks_list3, number_of_peaks_list4, number_of_peaks_list5, number_of_peaks_list6, number_of_peaks_list7, number_of_peaks_list8]


t = np.linspace(0, 8*(I_total) + I_off, 10*(8*(I_total) + (I_off)) + 1)


i = 0
while i < trials:

    #Each cell has a random normally distributed hmax
    #The boundary between mono stable low and excitable it about -0.128
    #The boundary for the excitable regime is -0.068
    mu_hmax = -.127 #This ensures ~52% of cells are placed in excitable regime and ~48% stay in mono low
    sigma_hmax = .0198
    

    hmax = np.random.normal(mu_hmax, sigma_hmax, number_of_cells)

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

        if t > ((t//I_total)*(I_total)) and t < ((t//I_total)*(I_total) + I_off):
            h = hmin     
        else:
            h = hmax 
        
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
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_obr - 2*g*x
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
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_obl + g*x_obr - 4*g*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)
                    
                    
                
                
                #last column
                elif j % row_size == 0:
                    x_l = state[2*(j-1) - 2]
                    x_obr = state[2*((j+(row_size))-1)]
                    x_obl = state[2*((j+(row_size-1))-1)]
                    a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_obl + g*x_obr - 3*g*x
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
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_etl + g*x_etr + g*x_ebl + g*x_ebr - 5*g*x
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
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_etl + g*x_etr + g*x_ebl + g*x_ebr - 6*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                
                    #last column
                    elif j % row_size == 0:
                        x_l = state[2*(j-1) - 2]
                        x_etl = state[2*((j-(row_size))-1)]
                        x_ebl = state[2*((j+(row_size))-1)]
                        a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_etl + g*x_ebl - 3*g*x
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
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_otr + g*x_obr - 3*g*x
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
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_otl + g*x_otr + g*x_obl + g*x_obr - 6*g*x
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
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_otl + g*x_otr + g*x_obl + g*x_obr - 5*g*x
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
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_etl + g*x_etr - 3*g*x
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
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_etl + g*x_etr - 4*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                
                    #last column
                    elif j % row_size == 0:
                        x_l = state[2*(j-1) - 2]
                        x_etl = state[2*((j-(row_size))-1)]
                        a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_etl - 2*g*x
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
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_otr - 2*g*x
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
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_otl + g*x_otr - 4*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt) 
                        
                        
                        
                
                    #last column
                    elif j % row_size == 0:
                        x_l = state[2*(j-1) - 2]
                        x_otr = state[2*((j-(row_size))-1)]
                        x_otl = state[2*((j-(row_size+1))-1)]
                        a = (N*(h[j-1] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_otl + g*x_otr - 3*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt) 
            
            
            
            
            j += 1
            
        
            
        return ODEs

    #Solving ODEs
    IC = IC(2*number_of_cells)
    X = odeint(MM, IC, t, args = (number_of_cells,k)).T
    
        
    
    
    #Getting numeric time series of the X molecules for each cell
    jj = 0
    x_list = np.array([t]) #t is just a placeholder so an array can be appended
    while jj < 2*number_of_cells:
        x_list = np.append(x_list, [X[jj]], axis = 0)
        jj += 2
    jj = 0
    #print(np.size(x_list[0]))
    
    
    #Finding peak times and average number of peaks for each cycles
    II = 10*(I_off) # sets interval where to search for peaks for a cycle
    list_i = 0 # index corresponding to the cycle
    while II < 10*(8*(I_total) + (I_off)):
        jj = 1 #starts at 1 because 0 element is t (t is just a placeholder)
        tmax_list = np.array([])
        number_of_peaks = np.array([])
###################################################################################
        while jj < number_of_cells + 1:
            xmax = np.max(x_list[jj][(II):(II + (10*I_total) + 1)])
            #print(xmax, x_list[jj][II], x_list[jj][II + (10*2*I)])
            #print(xmax)
            if xmax != x_list[jj][II] and xmax != x_list[jj][II + (10*I_total)]:
                #number_of_peaks = np.append(number_of_peaks, 1)
                tmax = np.where(x_list[jj][(II):(II + (10*I_total)) + 1] == xmax)
                #print(tmax)
                tmax_list = np.append(tmax_list, t[tmax])
            else:
                #number_of_peaks = np.append(number_of_peaks, 0)
                tmax_list = np.append(tmax_list, 0)
            if xmax > 250:
                number_of_peaks = np.append(number_of_peaks, 1)
            else:
                number_of_peaks = np.append(number_of_peaks, 0)
            jj += 1


        #Adding number of peaks past threshold in cycle (list_i)
        number_of_peaks_lists[list_i] = np.append(number_of_peaks_lists[list_i], np.sum(number_of_peaks))
        

            
        #Counting the number of edges between nearest-neighbors
        #j is the cell's position in the lattice. 0 starts at top left
        #j is also the cell's index in tmax_list
        #Cells are only compared to their neighbors down and to the right
        number_of_edges = np.array([])
        j = 0
        while j < np.size(tmax_list) - 1: #there is a (-1) because nothing is done for the last cell
        

            #Not last row
            if j // row_size != (row_size - 1):
                #Odd rows
                if (j // row_size) % 2 == 0:
                    
                    #1st column
                    if j % row_size == 0:
                        #right
                        if tmax_list[j] != 0 and tmax_list[j+1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+1]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+1]) <= 2.5:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                            
                        #bottom right
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= 2.5:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                   
                    #Middle columns
                    elif j % row_size != 0 and (j + 1) % row_size != 0:
                        #right
                        if tmax_list[j] != 0 and tmax_list[j+1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+1]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+1]) <= 2.5:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                        #bottom left
                        if tmax_list[j] != 0 and tmax_list[j+row_size-1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size-1]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+row_size-1]) <= 2.5:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                        #bottom right
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= 2.5:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
    
                    #Last column
                    elif (j + 1) % row_size == 0:
                        #bottom left
                        if tmax_list[j] != 0 and tmax_list[j+row_size-1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size-1]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+row_size-1]) <= 2.5:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                          
                        #bottom right
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= 2.5:
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
                            if np.absolute(tmax_list[j] - tmax_list[j+1]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+1]) <= 2.5:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                        #bottom left
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= 2.5:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                        #bottom right
                        if tmax_list[j] != 0 and tmax_list[j+row_size+1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size+1]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+row_size+1]) <= 2.5:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
                        
                    
                    #Last column
                    else: #(j + 1) % row_size == 0
                        #bottom left
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= 2.5:
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
                    if np.absolute(tmax_list[j] - tmax_list[j+1]) >= .5 and np.absolute(tmax_list[j] - tmax_list[j+1]) <= 2.5:
                        number_of_edges = np.append(number_of_edges, 1)
                    else:
                        number_of_edges = np.append(number_of_edges, 0)
                        
                else:
                    number_of_edges = np.append(number_of_edges, 0)
        
            j += 1
            

            







######################################################################################
        
        #The number of edges for each trial is added to number_of_edges_list
        number_of_edges_lists[list_i] = np.append(number_of_edges_lists[list_i], np.sum(number_of_edges))    
        #print(np.sum(number_of_edges))

        list_i += 1
        
        II += 10*I_total
        
    print(i)   
    i += 1



# finding the average number of edges for each cycle

edges1 = number_of_edges_lists[0]
edges2 = number_of_edges_lists[1]
edges3 = number_of_edges_lists[2]
edges4 = number_of_edges_lists[3]
edges5 = number_of_edges_lists[4]
edges6 = number_of_edges_lists[5]
edges7 = number_of_edges_lists[6]
edges8 = number_of_edges_lists[7]



#finding the average number of peaks for each cycle
peaks1 = number_of_peaks_lists[0]
peaks2 = number_of_peaks_lists[1]
peaks3 = number_of_peaks_lists[2]
peaks4 = number_of_peaks_lists[3]
peaks5 = number_of_peaks_lists[4]
peaks6 = number_of_peaks_lists[5]
peaks7 = number_of_peaks_lists[6]
peaks8 = number_of_peaks_lists[7]


# plotting results


#plt.figure(1)
#plt.scatter(cycle_list, ave_number_of_edges_list)
#plt.xlabel("Cycle")
#plt.ylabel("Average Number of Edges")
#plt.title(f"Average Number of Edges vs Cycle (Off = {I_off} On = {I_on})")



#plt.figure(2)
#plt.scatter(cycle_list, ave_number_of_peaks_list)
#plt.xlabel("Cycle")
#plt.ylabel("Average Number of Peaks")
#plt.title(f"Average Number of Peaks vs Cycle (Off = {I_off} On = {I_on})")


#plt.figure(3)
#plt.plot(t, X[0], label = "x1")
#plt.plot(t, X[2], label = "x2")
#plt.plot(t, X[4], label = "x3")
#plt.plot(t, X[6], label = "x4")
#plt.plot(t, X[8], label = "x5")
#plt.plot(t, X[10], label = "x6")
#plt.plot(t, X[12], label = "x7")
#plt.plot(t, X[14], label = "x8")
#plt.plot(t, X[16], label = "x9")
#plt.xlabel("time")
#plt.ylabel("Number of Molecules")
#plt.title("Deterministic x(t)")
#plt.legend()   

#plt.figure(4)
#plt.plot(t, X[0], label = "x1")
#plt.plot(t, X[2], label = "x2")
#plt.plot(t, X[4], label = "x3")
#plt.plot(t, X[6], label = "x4")
#plt.plot(t, X[8], label = "x5")
#plt.plot(t, X[10], label = "x6")
#plt.plot(t, X[12], label = "x7")
#plt.plot(t, X[14], label = "x8")
#plt.plot(t, X[16], label = "x9")
#plt.xlim(150,200)
#plt.xlabel("time")
#plt.ylabel("Number of Molecules")
#plt.title("Deterministic x(t)")
#plt.legend()


#plt.figure(5)
#plt.plot(t, X[0], label = "x1")
#plt.plot(t, X[2], label = "x2")
#plt.plot(t, X[4], label = "x3")
#plt.plot(t, X[6], label = "x4")
#plt.plot(t, X[8], label = "x5")
#plt.plot(t, X[10], label = "x6")
#plt.plot(t, X[12], label = "x7")
#plt.plot(t, X[14], label = "x8")
#plt.plot(t, X[16], label = "x9")
#plt.xlim(20,40)
#plt.xlabel("time")
#plt.ylabel("Number of Molecules")
#plt.title("Deterministic x(t)")
#plt.legend()

#Saving data
np.savez("324_Cells_g_025_8_Cycles_40_40_Results", edges1 = edges1, edges2 = edges2, edges3 = edges3, edges4 = edges4, edges5 = edges5, edges6 = edges6, edges7 = edges7, edges8 = edges8, peaks1 = peaks1, peaks2 = peaks2, peaks3 = peaks3, peaks4 = peaks4, peaks5 = peaks5, peaks6 = peaks6, peaks7 = peaks7, peaks8 = peaks8)

 
