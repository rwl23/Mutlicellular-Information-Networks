import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

##########################################################################
#Arbitrary number of cells in a triangular lattice stimulated by a square 
#wave. Each cells has a randomly choosen hmax value

##########################################################################
#Number of trials
trials = 50

#Total number of cells in lattice
row_size = 18 
number_of_cells = row_size**2
number_of_possible_edges = (3*row_size-1)*(row_size-1)

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


period = 10
t_total = 16*period
t = np.linspace(0,t_total, 10*t_total)

Pedge_list = np.array([])
FiringFrac_list = np.array([])

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
    

    tau_list = np.array([0])
    while np.size(tau_list) < 5: # I don't want to deal with trials with few periods
        tau_list = np.array([0])
        tau = 0
        while tau < t_total:
            r1 = np.random.rand()
            r2 = period*np.log(1/r1)
            tau += r2
            tau_list = np.append(tau_list, tau)



    # h_list is for printing to see when h changes
    h_list = np.array([])
    for i_t in t:
        i_tau = 0
        while i_tau < np.size(tau_list):
            if i_t >= tau_list[i_tau]:
                i_tau += 1
            else:
                if i_tau % 2 != 0:
                    h_list = np.append(h_list,0)
                    break
                else:
                    h_list = np.append(h_list,100)
                    break



    
    
    def MM(state,t, number_of_cells, k):
        ODEs = np.array([])


        i_tau = 0
        while i_tau < np.size(tau_list):
            if t >= tau_list[i_tau]:
                i_tau += 1
            else:
                if i_tau % 2 != 0:
                    h = hmin
                    break
                else:
                    h = hmax
                    break

        
        j = 0
        while j < number_of_cells:
            x = state[2*j]
            y = state[2*j + 1] 
            
            
            
            #x_r = state[2*j + 2]
            #x_l = state[2*j - 2]
            
            #For odd rows
            #x_otr = state[2*(j-row_size)]
            #x_otl = state[2*(j-(row_size+1))]
            #x_obr = state[2*(j+(row_size))]
            #x_obl = state[2*(j+(row_size-1))]
            
            #For even rows
            #x_etr = state[2*(j-(row_size-1))]
            #x_etl = state[2*(j-(row_size))]
            #x_ebr = state[2*(j+(row_size+1))]
            #x_ebl = state[2*(j+(row_size))]
            
            
            
            #1st row
            #Note: First row is always odd
            if j // row_size == 0:
                #1st column
                if j % row_size == 0:
                    x_r = state[2*j + 2]
                    x_obr = state[2*(j+(row_size))]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_obr - 2*g*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)
                    
                    
                
                
                #middle columns
                elif j % row_size != 0 and j % row_size != (row_size-1):
                    x_r = state[2*j + 2]
                    x_l = state[2*j - 2]
                    x_obr = state[2*(j+(row_size))]
                    x_obl = state[2*(j+(row_size-1))]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_obl + g*x_obr - 4*g*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)
                    
                    
                
                
                #last column
                elif j % row_size == (row_size-1):
                    x_l = state[2*j - 2]
                    x_obr = state[2*(j+(row_size))]
                    x_obl = state[2*(j+(row_size-1))]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_obl + g*x_obr - 3*g*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)
                    
                
                
                
            #middle rows
            elif j // row_size != 0 and j // row_size != (row_size - 1):
                
                #even rows
                if (j // row_size) % 2 != 0:
                    
                    #1st column
                    if j % row_size == 0:
                        x_r = state[2*j + 2]
                        x_etr = state[2*(j-(row_size-1))]
                        x_etl = state[2*(j-(row_size))]
                        x_ebr = state[2*(j+(row_size+1))]
                        x_ebl = state[2*(j+(row_size))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_etl + g*x_etr + g*x_ebl + g*x_ebr - 5*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)    
                        
                        
                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        x_r = state[2*j + 2]
                        x_l = state[2*j - 2]
                        x_etr = state[2*(j-(row_size-1))]
                        x_etl = state[2*(j-(row_size))]
                        x_ebr = state[2*(j+(row_size+1))]
                        x_ebl = state[2*(j+(row_size))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_etl + g*x_etr + g*x_ebl + g*x_ebr - 6*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        x_etl = state[2*(j-(row_size))]
                        x_ebl = state[2*(j+(row_size))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_etl + g*x_ebl - 3*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                        
    
                #odd rows
                else: #(j  // row_size) % 2 == 0
                
                    #1st column
                    if j % row_size == 0:
                        x_r = state[2*j + 2]
                        x_otr = state[2*(j-(row_size))]
                        x_obr = state[2*(j+(row_size))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_otr + g*x_obr - 3*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                        
                
                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        x_r = state[2*j + 2]
                        x_l = state[2*j - 2]
                        x_otr = state[2*(j-(row_size))]
                        x_otl = state[2*(j-(row_size+1))]
                        x_obr = state[2*(j+(row_size))]
                        x_obl = state[2*(j+(row_size-1))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_otl + g*x_otr + g*x_obl + g*x_obr - 6*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                        
                
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        x_otr = state[2*(j-(row_size))]
                        x_otl = state[2*(j-(row_size+1))]
                        x_obr = state[2*(j+(row_size))]
                        x_obl = state[2*(j+(row_size-1))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_otl + g*x_otr + g*x_obl + g*x_obr - 5*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
            
            
            
            #last row
            elif j // row_size == (row_size - 1):
                
                #even rows
                if (j // row_size) % 2 != 0:
                    
                    #1st column
                    if j % row_size == 0:
                        x_r = state[2*j + 2]
                        x_etr = state[2*(j-(row_size-1))]
                        x_etl = state[2*(j-(row_size))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_etl + g*x_etr - 3*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
    
                
                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-1):
                        x_r = state[2*j + 2]
                        x_l = state[2*j - 2]
                        x_etr = state[2*(j-(row_size-1))]
                        x_etl = state[2*(j-(row_size))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_etl + g*x_etr - 4*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        x_etl = state[2*(j-(row_size))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_etl - 2*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)   
                        
                        
                        
                    
                #odd rows
                else: #(j // row_size) % 2 == 0
                
                    #1st column
                    if j % row_size == 0:
                        x_r = state[2*j + 2]
                        x_otr = state[2*(j-(row_size))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_r + g*x_otr - 2*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt) 
                        
                
                    #middle columns
                    elif j  % row_size != 0 and j % row_size != (row_size-1):
                        x_r = state[2*j + 2]
                        x_l = state[2*j - 2]
                        x_otr = state[2*(j-(row_size))]
                        x_otl = state[2*(j-(row_size+1))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + g*x_l + g*x_r + g*x_otl + g*x_otr - 4*g*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt) 
                        
                        
                        
                
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        x_otr = state[2*(j-(row_size))]
                        x_otl = state[2*(j-(row_size+1))]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
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
    
    #Array to hold Pedge of all periods in a trial
    Pedge_per_period_list = np.array([])

    #Array to hold FiringFrac of all periods in a trail
    FiringFrac_per_period_list = np.array([])


    ii_tau = 1 #the first element in tau_list is when hmax first comes on
    if np.size(tau_list) % 2 == 0:
        ii_tau_max = np.size(tau_list) - 3
    else:
        ii_tau_max = np.size(tau_list) - 2

    
    while ii_tau + 2 <= ii_tau_max:
        I = int(10*tau_list[ii_tau])
        II = int(10*tau_list[ii_tau+2])
        jj = 1 #starts at 1 because 0 element is t (t is just a placeholder)
        tmax_list = np.array([])
        number_of_peaks = 0
        while jj < number_of_cells + 1:
            xmax = np.max(x_list[jj][I:II])
            if xmax != x_list[jj][I] and xmax != x_list[jj][II-1]:
                #number_of_peaks = np.append(number_of_peaks, 1)
                tmax = np.where(x_list[jj][I:II] == xmax)
                #print(tmax)
                tmax_list = np.append(tmax_list, t[tmax])
            else:
                #number_of_peaks = np.append(number_of_peaks, 0)
                tmax_list = np.append(tmax_list, 0)
            if xmax > 250:
                number_of_peaks += 1
            else:
                number_of_peaks += 0
            jj += 1

        #Adding periods firing fraction to list
        FiringFrac_per_period_list = np.append(FiringFrac_per_period_list, number_of_peaks/number_of_cells)



        #Counting the number of edges between nearest-neighbors
        #j is the cell's position in the lattice. 0 starts at top left
        #j is also the cell's index in tmax_list
        #Cells are only compared to their neighbors down and to the right
        number_of_edges = np.array([])
        j = 0
        while j < number_of_cells - 1: #there is a (-1) because nothing is done for the last cell
        

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

        Pedge_per_period_list = np.append(Pedge_per_period_list, sum(number_of_edges)/number_of_possible_edges)



        ii_tau += 2


    #print(Pedge_per_period_list)
    #print(FiringFrac_per_period_list)

    Pedge_list = np.append(Pedge_list, np.average(Pedge_per_period_list))
    FiringFrac_list = np.append(FiringFrac_list, np.average(FiringFrac_per_period_list))


    # plotting results
    #plt.figure()
    #plt.plot(t,h_list,c = 'black')
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
    #plt.show()

    i += 1

#print(Pedge_list)
#print(FiringFrac_list)

#Saving data
np.savez("324_Cells_Exponential_Distribution_10_10_Results", Pedge_list = Pedge_list, FiringFrac_list = FiringFrac_list)

 