import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

##########################################################################
#Arbitrary number of cells in a disordered lattice stimulated by a square 
#wave. Each cells has a randomly choosen hmax value

##########################################################################
#Number of trials
trials = 50

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

g = np.array([.025]) 
gi = 0

p_couple_nn = .8
p_couple_nnn = .15

I_off = 30
I_on = 30
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


    #Creating arrays to see if links are there or not
    #note that extra links that are never used are added onto 
    # certain nodes. Doing this simplifies the indexing
    #left-right links: 
    num_lr_links = row_size*row_size
    r_lr = np.random.rand(num_lr_links)
    for ri in range(np.size(r_lr)):
        if r_lr[ri] < p_couple_nn:
            r_lr[ri] = 1
        else:
            r_lr[ri] = 0

    #bottom-right/top-left links:
    num_brtl_links = (row_size - 1)*row_size
    r_brtl = np.random.rand(num_brtl_links)
    for ri in range(np.size(r_brtl)):
        if r_brtl[ri] < p_couple_nn:
            r_brtl[ri] = 1
        else:
            r_brtl[ri] = 0

    #bottom-left/top-right links:
    num_bltr_links = (row_size - 1)*row_size
    r_bltr = np.random.rand(num_bltr_links)
    for ri in range(np.size(r_bltr)):
        if r_bltr[ri] < p_couple_nn:
            r_bltr[ri] = 1
        else:
            r_bltr[ri] = 0
    

    #next nearest neighbor down/up links:
    num_nnn_du_links = (row_size - 2)*row_size
    r_nnn_du = np.random.rand(num_nnn_du_links)
    for ri in range(np.size(r_nnn_du)):
        if r_nnn_du[ri] < p_couple_nnn:
            r_nnn_du[ri] = 1
        else:
            r_nnn_du[ri] = 0
        
    
    #next nearest neighbor bottom-right/top-left links:
    num_nnn_brtl_links = (row_size - 1)*row_size
    r_nnn_brtl = np.random.rand(num_nnn_brtl_links)
    for ri in range(np.size(r_nnn_brtl)):
        if r_nnn_brtl[ri] < p_couple_nnn:
            r_nnn_brtl[ri] = 1
        else:
            r_nnn_brtl[ri] = 0

    #next nearest neighbor bottom-left/top-right links:
    num_nnn_bltr_links = (row_size - 1)*row_size
    r_nnn_bltr = np.random.rand(num_nnn_bltr_links)
    for ri in range(np.size(r_nnn_bltr)):
        if r_nnn_bltr[ri] < p_couple_nnn:
            r_nnn_bltr[ri] = 1
        else:
            r_nnn_bltr[ri] = 0





    
    
    #Creating array of initial conditions. Note: they are all 0
    def IC(i):
        IC_list = np.array([])
        j = 0
        while j < i:
            IC_list = np.append(IC_list, 0)
            j += 1
        return IC_list
    
    
    
    def MM(state,t, num_of_cells, k):
        ODEs = np.array([])

        if t > ((t//I_total)*(I_total)) and t < ((t//I_total)*(I_total) + I_off):
            h = hmin     
        else:
            h = hmax 
        
        j = 0
        while j < num_of_cells:
            x = state[2*j]
            y = state[2*j + 1] 
            
            #x_l = state[2*j - 2]
            #l_l = r_lr[j-1]
            #x_r = state[2*j + 2]
            #l_r = r_lr[j]
            
            #x_nnnd = state[2*(j+2*row_size)]
            #l_nnnd = r_nnn_du[j]
            #x_nnnu = state[2*(j-2*row_size)]
            #l_nnnu = r_nnn_du[j-2*row_size]
            
            #For odd rows
            #x_otl = state[2*(j-(row_size+1))]
            #l_otl = r_brtl[j-(row_size+1)]
            #x_otr = state[2*(j-row_size)]
            #l_otr = r_bltr[j-row_size]
            #x_obl = state[2*(j+(row_size-1))]
            #l_obl = r_bltr[j]
            #x_obr = state[2*(j+row_size)]
            #l_obr = r_brtl[j]
            #x_onnntl = state[2*(j-(row_size+2))]
            #l_onnntl = r_nnn_brtl[j-(row_size+2)]
            #x_onnntr = state[2*(j-(row_size-1))]
            #l_onnntr = r_nnn_bltr[j-(row_size-1)]
            #x_onnnbl = state[2*(j+(row_size-2))]
            #l_onnnbl = r_nnn_bltr[j]
            #x_onnnbr = state[2*(j+(row_size+1))]
            #l_onnnbr = r_nnn_brtl[j]
            
            #For even rows
            #x_etl = state[2*(j-row_size)]
            #l_etl = r_brtl[j-row_size]
            #x_etr = state[2*(j-(row_size-1))]
            #l_etr = r_bltr[j-(row_size-1)]
            #x_ebl = state[2*(j+row_size)]
            #l_ebl = r_bltr[j]
            #x_ebr = state[2*(j+(row_size+1))]
            #l_ebr = r_brtl[j]
            #x_ennntl = state[2*(j-(row_size+1))]
            #l_ennntl = r_nnn_brtl[j-(row_size+1)]
            #x_ennntr = state[2*(j-(row_size-2))]
            #l_ennntr = r_nnn_bltr[j-(row_size-2)]
            #x_ennnbl = state[2*(j+(row_size-1))]
            #l_ennnbl = r_nnn_bltr[j]
            #x_ennnbr = state[2*(j+(row_size+2))]
            #l_ennnbr = r_nnn_brtl[j]
            
            
            
            #1st row
            #Note: First row is always odd
            if j  // row_size == 0:
                #1st column
                if j % row_size == 0:
                    x_r = state[2*j + 2]
                    l_r = r_lr[j]
                    x_obr = state[2*(j+row_size)]
                    l_obr = r_brtl[j]
                    x_nnnd = state[2*(j+2*row_size)]
                    l_nnnd = r_nnn_du[j]
                    x_onnnbr = state[2*(j+(row_size+1))]
                    l_onnnbr = r_nnn_brtl[j]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_r*g[gi]*x_r + l_obr*g[gi]*x_obr + l_nnnd*g[gi]*x_nnnd + l_onnnbr*g[gi]*x_onnnbr - (l_r+l_obr+l_nnnd+l_onnnbr)*g[gi]*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)
                

                #2nd column
                elif j % row_size == 1:
                    x_l = state[2*j - 2]
                    l_l = r_lr[j-1]
                    x_r = state[2*j + 2]
                    l_r = r_lr[j]
                    x_obl = state[2*(j+(row_size-1))]
                    l_obl = r_bltr[j]
                    x_obr = state[2*(j+row_size)]
                    l_obr = r_brtl[j]
                    x_nnnd = state[2*(j+2*row_size)]
                    l_nnnd = r_nnn_du[j]
                    x_onnnbr = state[2*(j+(row_size+1))]
                    l_onnnbr = r_nnn_brtl[j]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_obl*g[gi]*x_obl + l_obr*g[gi]*x_obr + l_nnnd*g[gi]*x_nnnd + l_onnnbr*g[gi]*x_onnnbr - (l_l+l_r+l_obl+l_obr+l_nnnd+l_onnnbr)*g[gi]*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)
                    
                
                
                #middle columns
                elif j % row_size != 0 and j % row_size != 1 and j % row_size != (row_size-1):
                    x_l = state[2*j - 2]
                    l_l = r_lr[j-1]
                    x_r = state[2*j + 2]
                    l_r = r_lr[j]
                    x_obl = state[2*(j+(row_size-1))]
                    l_obl = r_bltr[j]
                    x_obr = state[2*(j+row_size)]
                    l_obr = r_brtl[j]
                    x_nnnd = state[2*(j+2*row_size)]
                    l_nnnd = r_nnn_du[j]
                    x_onnnbl = state[2*(j+(row_size-2))]
                    l_onnnbl = r_nnn_bltr[j]
                    x_onnnbr = state[2*(j+(row_size+1))]
                    l_onnnbr = r_nnn_brtl[j]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_obl*g[gi]*x_obl + l_obr*g[gi]*x_obr + l_nnnd*g[gi]*x_nnnd + l_onnnbl*g[gi]*x_onnnbl + l_onnnbr*g[gi]*x_onnnbr - (l_l+l_r+l_obl+l_obr+l_nnnd+l_onnnbl+l_onnnbr)*g[gi]*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)

                
                #last column
                elif j % row_size == (row_size-1):
                    x_l = state[2*j - 2]
                    l_l = r_lr[j-1]
                    x_obl = state[2*(j+(row_size-1))]
                    l_obl = r_bltr[j]
                    x_obr = state[2*(j+row_size)]
                    l_obr = r_brtl[j]
                    x_nnnd = state[2*(j+2*row_size)]
                    l_nnnd = r_nnn_du[j]
                    x_onnnbl = state[2*(j+(row_size-2))]
                    l_onnnbl = r_nnn_bltr[j]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_obl*g[gi]*x_obl + l_obr*g[gi]*x_obr + l_nnnd*g[gi]*x_nnnd + l_onnnbl*g[gi]*x_onnnbl - (l_l+l_obl+l_obr+l_nnnd+l_onnnbl)*g[gi]*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)
                    
                
            

            #2nd row
            #Note: Second row is always even
            elif j // row_size == 1:
                #1st column
                if j % row_size == 0:
                    x_r = state[2*j + 2]
                    l_r = r_lr[j]
                    x_etl = state[2*(j-row_size)]
                    l_etl = r_brtl[j-row_size]
                    x_etr = state[2*(j-(row_size-1))]
                    l_etr = r_bltr[j-(row_size-1)]
                    x_ebl = state[2*(j+row_size)]
                    l_ebl = r_bltr[j]
                    x_ebr = state[2*(j+(row_size+1))]
                    l_ebr = r_brtl[j]
                    x_nnnd = state[2*(j+2*row_size)]
                    l_nnnd = r_nnn_du[j]
                    x_ennntr = state[2*(j-(row_size-2))]
                    l_ennntr = r_nnn_bltr[j-(row_size-2)]
                    x_ennnbr = state[2*(j+(row_size+2))]
                    l_ennnbr = r_nnn_brtl[j]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_ebl*g[gi]*x_ebl + l_ebr*g[gi]*x_ebr + l_nnnd*g[gi]*x_nnnd + l_ennntr*g[gi]*x_ennntr + l_ennnbr*g[gi]*x_ennnbr - (l_r+l_etl+l_etr+l_ebl+l_ebr+l_nnnd+l_ennntr+l_ennnbr)*g[gi]*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)

                #middle columns
                elif j % row_size != 0 and j % row_size != (row_size-2) and j % row_size != (row_size-1):
                    x_l = state[2*j - 2]
                    l_l = r_lr[j-1]
                    x_r = state[2*j + 2]
                    l_r = r_lr[j]
                    x_etl = state[2*(j-row_size)]
                    l_etl = r_brtl[j-row_size]
                    x_etr = state[2*(j-(row_size-1))]
                    l_etr = r_bltr[j-(row_size-1)]
                    x_ebl = state[2*(j+row_size)]
                    l_ebl = r_bltr[j]
                    x_ebr = state[2*(j+(row_size+1))]
                    l_ebr = r_brtl[j]
                    x_nnnd = state[2*(j+2*row_size)]
                    l_nnnd = r_nnn_du[j]
                    x_ennntl = state[2*(j-(row_size+1))]
                    l_ennntl = r_nnn_brtl[j-(row_size+1)]
                    x_ennntr = state[2*(j-(row_size-2))]
                    l_ennntr = r_nnn_bltr[j-(row_size-2)]
                    x_ennnbl = state[2*(j+(row_size-1))]
                    l_ennnbl = r_nnn_bltr[j]
                    x_ennnbr = state[2*(j+(row_size+2))]
                    l_ennnbr = r_nnn_brtl[j]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_ebl*g[gi]*x_ebl + l_ebr*g[gi]*x_ebr + l_nnnd*g[gi]*x_nnnd + l_ennntl*g[gi]*x_ennntl + l_ennntr*g[gi]*x_ennntr + l_ennnbl*g[gi]*x_ennnbl + l_ennnbr*g[gi]*x_ennnbr - (l_l+l_r+l_etl+l_etr+l_ebl+l_ebr+l_nnnd+l_ennntl+l_ennntr+l_ennnbl+l_ennnbr)*g[gi]*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)

                #2nd to last column
                elif j % row_size == (row_size-2):
                    x_l = state[2*j - 2]
                    l_l = r_lr[j-1]
                    x_r = state[2*j + 2]
                    l_r = r_lr[j]
                    x_etl = state[2*(j-row_size)]
                    l_etl = r_brtl[j-row_size]
                    x_etr = state[2*(j-(row_size-1))]
                    l_etr = r_bltr[j-(row_size-1)]
                    x_ebl = state[2*(j+row_size)]
                    l_ebl = r_bltr[j]
                    x_ebr = state[2*(j+(row_size+1))]
                    l_ebr = r_brtl[j]
                    x_nnnd = state[2*(j+2*row_size)]
                    l_nnnd = r_nnn_du[j]
                    x_ennntl = state[2*(j-(row_size+1))]
                    l_ennntl = r_nnn_brtl[j-(row_size+1)]
                    x_ennnbl = state[2*(j+(row_size-1))]
                    l_ennnbl = r_nnn_bltr[j]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_ebl*g[gi]*x_ebl + l_ebr*g[gi]*x_ebr + l_nnnd*g[gi]*x_nnnd + l_ennntl*g[gi]*x_ennntl + l_ennnbl*g[gi]*x_ennnbl - (l_l+l_r+l_etl+l_etr+l_ebl+l_ebr+l_nnnd+l_ennntl+l_ennnbl)*g[gi]*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)

                #last column
                elif j % row_size == (row_size-1):
                    x_l = state[2*j - 2]
                    l_l = r_lr[j-1]
                    x_etl = state[2*(j-row_size)]
                    l_etl = r_brtl[j-row_size]
                    x_ebl = state[2*(j+row_size)]
                    l_ebl = r_bltr[j]
                    x_nnnd = state[2*(j+2*row_size)]
                    l_nnnd = r_nnn_du[j]
                    x_ennntl = state[2*(j-(row_size+1))]
                    l_ennntl = r_nnn_brtl[j-(row_size+1)]
                    x_ennnbl = state[2*(j+(row_size-1))]
                    l_ennnbl = r_nnn_bltr[j]
                    a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                    dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_etl*g[gi]*x_etl + l_ebl*g[gi]*x_ebl + l_nnnd*g[gi]*x_nnnd + l_ennntl*g[gi]*x_ennntl + l_ennnbl*g[gi]*x_ennnbl- (l_l+l_etl+l_ebl+l_nnnd+l_ennntl+l_ennnbl)*g[gi]*x
                    ODEs = np.append(ODEs, dx1dt)
                    dy1dt = d*x - e*y
                    ODEs = np.append(ODEs, dy1dt)

                
            #middle rows
            elif j // row_size != 0 and j // row_size != 1 and j // row_size != (row_size - 2) and j // row_size != (row_size - 1):
                
                #even rows
                if (j  // row_size) % 2 != 0:
                    
                    #1st column
                    if j % row_size == 0:
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_etr = state[2*(j-(row_size-1))]
                        l_etr = r_bltr[j-(row_size-1)]
                        x_ebl = state[2*(j+row_size)]
                        l_ebl = r_bltr[j]
                        x_ebr = state[2*(j+(row_size+1))]
                        l_ebr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_nnnd = state[2*(j+2*row_size)]
                        l_nnnd = r_nnn_du[j]
                        x_ennntr = state[2*(j-(row_size-2))]
                        l_ennntr = r_nnn_bltr[j-(row_size-2)]
                        x_ennnbr = state[2*(j+(row_size+2))]
                        l_ennnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_ebl*g[gi]*x_ebl + l_ebr*g[gi]*x_ebr + l_nnnu*g[gi]*x_nnnu + l_nnnd*g[gi]*x_nnnd + l_ennntr*g[gi]*x_ennntr + l_ennnbr*g[gi]*x_ennnbr - (l_r+l_etl+l_etr+l_ebl+l_ebr+l_nnnu+l_nnnd+l_ennntr+l_ennnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)   
                        
                        
                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-2) and j % row_size != (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_etr = state[2*(j-(row_size-1))]
                        l_etr = r_bltr[j-(row_size-1)]
                        x_ebl = state[2*(j+row_size)]
                        l_ebl = r_bltr[j]
                        x_ebr = state[2*(j+(row_size+1))]
                        l_ebr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_nnnd = state[2*(j+2*row_size)]
                        l_nnnd = r_nnn_du[j]
                        x_ennntl = state[2*(j-(row_size+1))]
                        l_ennntl = r_nnn_brtl[j-(row_size+1)]
                        x_ennntr = state[2*(j-(row_size-2))]
                        l_ennntr = r_nnn_bltr[j-(row_size-2)]
                        x_ennnbl = state[2*(j+(row_size-1))]
                        l_ennnbl = r_nnn_bltr[j]
                        x_ennnbr = state[2*(j+(row_size+2))]
                        l_ennnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_ebl*g[gi]*x_ebl + l_ebr*g[gi]*x_ebr + l_nnnu*g[gi]*x_nnnu + l_nnnd*g[gi]*x_nnnd + l_ennntl*g[gi]*x_ennntl + l_ennntr*g[gi]*x_ennntr + l_ennnbl*g[gi]*x_ennnbl + l_ennnbr*g[gi]*x_ennnbr - (l_l+l_r+l_etl+l_etr+l_ebl+l_ebr+l_nnnu+l_nnnd+l_ennntl+l_ennntr+l_ennnbl+l_ennnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)

                    
                    #2nd to last column
                    elif j % row_size == (row_size-2):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_etr = state[2*(j-(row_size-1))]
                        l_etr = r_bltr[j-(row_size-1)]
                        x_ebl = state[2*(j+row_size)]
                        l_ebl = r_bltr[j]
                        x_ebr = state[2*(j+(row_size+1))]
                        l_ebr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_nnnd = state[2*(j+2*row_size)]
                        l_nnnd = r_nnn_du[j]
                        x_ennntl = state[2*(j-(row_size+1))]
                        l_ennntl = r_nnn_brtl[j-(row_size+1)]
                        x_ennnbl = state[2*(j+(row_size-1))]
                        l_ennnbl = r_nnn_bltr[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_ebl*g[gi]*x_ebl + l_ebr*g[gi]*x_ebr + l_nnnu*g[gi]*x_nnnu + l_nnnd*g[gi]*x_nnnd + l_ennntl*g[gi]*x_ennntl + l_ennnbl*g[gi]*x_ennnbl - (l_l+l_r+l_etl+l_etr+l_ebl+l_ebr+l_nnnu+l_nnnd+l_ennntl+l_ennnbl)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_ebl = state[2*(j+row_size)]
                        l_ebl = r_bltr[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_nnnd = state[2*(j+2*row_size)]
                        l_nnnd = r_nnn_du[j]
                        x_ennntl = state[2*(j-(row_size+1))]
                        l_ennntl = r_nnn_brtl[j-(row_size+1)]
                        x_ennnbl = state[2*(j+(row_size-1))]
                        l_ennnbl = r_nnn_bltr[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_etl*g[gi]*x_etl + l_ebl*g[gi]*x_ebl + l_nnnu*g[gi]*x_nnnu + l_nnnd*g[gi]*x_nnnd + l_ennntl*g[gi]*x_ennntl + l_ennnbl*g[gi]*x_ennnbl- (l_l+l_etl+l_ebl+l_nnnu+l_nnnd+l_ennntl+l_ennnbl)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                        
    
                #odd rows
                else: #(j // row_size) % 2 == 0
                
                    #1st column
                    if j % row_size == 0:
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_obr = state[2*(j+row_size)]
                        l_obr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_nnnd = state[2*(j+2*row_size)]
                        l_nnnd = r_nnn_du[j]
                        x_onnntr = state[2*(j-(row_size-1))]
                        l_onnntr = r_nnn_bltr[j-(row_size-1)]
                        x_onnnbr = state[2*(j+(row_size+1))]
                        l_onnnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_r*g[gi]*x_r + l_otr*g[gi]*x_otr + l_obr*g[gi]*x_obr + l_nnnu*g[gi]*x_nnnu + l_nnnd*g[gi]*x_nnnd + l_onnntr*g[gi]*x_onnntr + l_onnnbr*g[gi]*x_onnnbr - (l_r+l_otr+l_obr+l_nnnu+l_nnnd+l_onnntr+l_onnnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)

                        
                
                    #2nd column
                    elif j % row_size == 1:
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_otl = state[2*(j-(row_size+1))]
                        l_otl = r_brtl[j-(row_size+1)]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_obl = state[2*(j+(row_size-1))]
                        l_obl = r_bltr[j]
                        x_obr = state[2*(j+row_size)]
                        l_obr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_nnnd = state[2*(j+2*row_size)]
                        l_nnnd = r_nnn_du[j]
                        x_onnntr = state[2*(j-(row_size-1))]
                        l_onnntr = r_nnn_bltr[j-(row_size-1)]
                        x_onnnbr = state[2*(j+(row_size+1))]
                        l_onnnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_otl*g[gi]*x_otl + l_otr*g[gi]*x_otr + l_obl+g[gi]*x_obl +l_obr*g[gi]*x_obr + l_nnnu*g[gi]*x_nnnu + l_nnnd*g[gi]*x_nnnd + l_onnntr*g[gi]*x_onnntr + l_onnnbr*g[gi]*x_onnnbr - (l_l+l_r+l_otl+l_otr+l_obl+l_obr+l_nnnu+l_nnnd+l_onnntr+l_onnnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)


                    #middle columns
                    elif j % row_size != 0 and j % row_size != 1 and j % row_size != (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_otl = state[2*(j-(row_size+1))]
                        l_otl = r_brtl[j-(row_size+1)]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_obl = state[2*(j+(row_size-1))]
                        l_obl = r_bltr[j]
                        x_obr = state[2*(j+row_size)]
                        l_obr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_nnnd = state[2*(j+2*row_size)]
                        l_nnnd = r_nnn_du[j]
                        x_onnntl = state[2*(j-(row_size+2))]
                        l_onnntl = r_nnn_brtl[j-(row_size+2)]
                        x_onnntr = state[2*(j-(row_size-1))]
                        l_onnntr = r_nnn_bltr[j-(row_size-1)]
                        x_onnnbl = state[2*(j+(row_size-2))]
                        l_onnnbl = r_nnn_bltr[j]
                        x_onnnbr = state[2*(j+(row_size+1))]
                        l_onnnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_otl*g[gi]*x_otl + l_otr*g[gi]*x_otr + l_obl+g[gi]*x_obl +l_obr*g[gi]*x_obr + l_nnnu*g[gi]*x_nnnu + l_nnnd*g[gi]*x_nnnd + l_onnntl*g[gi]*x_onnntl + l_onnntr*g[gi]*x_onnntr + l_onnnbl*g[gi]*x_onnnbl + l_onnnbr*g[gi]*x_onnnbr - (l_l+l_r+l_otl+l_otr+l_obl+l_obr+l_nnnu+l_nnnd+l_onnntl+l_onnntr+l_onnnbl+l_onnnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_otl = state[2*(j-(row_size+1))]
                        l_otl = r_brtl[j-(row_size+1)]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_obl = state[2*(j+(row_size-1))]
                        l_obl = r_bltr[j]
                        x_obr = state[2*(j+row_size)]
                        l_obr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_nnnd = state[2*(j+2*row_size)]
                        l_nnnd = r_nnn_du[j]
                        x_onnntl = state[2*(j-(row_size+2))]
                        l_onnntl = r_nnn_brtl[j-(row_size+2)]
                        x_onnnbl = state[2*(j+(row_size-2))]
                        l_onnnbl = r_nnn_bltr[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_otl*g[gi]*x_otl + l_otr*g[gi]*x_otr + l_obl+g[gi]*x_obl +l_obr*g[gi]*x_obr + l_nnnu*g[gi]*x_nnnu + l_nnnd*g[gi]*x_nnnd + l_onnntl*g[gi]*x_onnntl + l_onnnbl*g[gi]*x_onnnbl - (l_l+l_otl+l_otr+l_obl+l_obr+l_nnnu+l_nnnd+l_onnntl+l_onnnbl)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
            #2nd to last row
            elif j // row_size == (row_size - 2):

                #even rows
                if (j // row_size) % 2 != 0:
                    
                    #1st column
                    if j % row_size == 0:
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_etr = state[2*(j-(row_size-1))]
                        l_etr = r_bltr[j-(row_size-1)]
                        x_ebl = state[2*(j+row_size)]
                        l_ebl = r_bltr[j]
                        x_ebr = state[2*(j+(row_size+1))]
                        l_ebr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_ennntr = state[2*(j-(row_size-2))]
                        l_ennntr = r_nnn_bltr[j-(row_size-2)]
                        x_ennnbr = state[2*(j+(row_size+2))]
                        l_ennnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_ebl*g[gi]*x_ebl + l_ebr*g[gi]*x_ebr + l_nnnu*g[gi]*x_nnnu + l_ennntr*g[gi]*x_ennntr + l_ennnbr*g[gi]*x_ennnbr - (l_r+l_etl+l_etr+l_ebl+l_ebr+l_nnnu+l_ennntr+l_ennnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)

                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-2) and j % row_size != (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_etr = state[2*(j-(row_size-1))]
                        l_etr = r_bltr[j-(row_size-1)]
                        x_ebl = state[2*(j+row_size)]
                        l_ebl = r_bltr[j]
                        x_ebr = state[2*(j+(row_size+1))]
                        l_ebr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_ennntl = state[2*(j-(row_size+1))]
                        l_ennntl = r_nnn_brtl[j-(row_size+1)]
                        x_ennntr = state[2*(j-(row_size-2))]
                        l_ennntr = r_nnn_bltr[j-(row_size-2)]
                        x_ennnbl = state[2*(j+(row_size-1))]
                        l_ennnbl = r_nnn_bltr[j]
                        x_ennnbr = state[2*(j+(row_size+2))]
                        l_ennnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_ebl*g[gi]*x_ebl + l_ebr*g[gi]*x_ebr + l_nnnu*g[gi]*x_nnnu + l_ennntl*g[gi]*x_ennntl + l_ennntr*g[gi]*x_ennntr + l_ennnbl*g[gi]*x_ennnbl + l_ennnbr*g[gi]*x_ennnbr - (l_l+l_r+l_etl+l_etr+l_ebl+l_ebr+l_nnnu+l_ennntl+l_ennntr+l_ennnbl+l_ennnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)

                    #2nd to last column
                    elif j % row_size == (row_size-2):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_etr = state[2*(j-(row_size-1))]
                        l_etr = r_bltr[j-(row_size-1)]
                        x_ebl = state[2*(j+row_size)]
                        l_ebl = r_bltr[j]
                        x_ebr = state[2*(j+(row_size+1))]
                        l_ebr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_ennntl = state[2*(j-(row_size+1))]
                        l_ennntl = r_nnn_brtl[j-(row_size+1)]
                        x_ennnbl = state[2*(j+(row_size-1))]
                        l_ennnbl = r_nnn_bltr[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_ebl*g[gi]*x_ebl + l_ebr*g[gi]*x_ebr + l_nnnu*g[gi]*x_nnnu + l_ennntl*g[gi]*x_ennntl + l_ennnbl*g[gi]*x_ennnbl - (l_l+l_r+l_etl+l_etr+l_ebl+l_ebr+l_nnnu+l_ennntl+l_ennnbl)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_ebl = state[2*(j+row_size)]
                        l_ebl = r_bltr[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_ennntl = state[2*(j-(row_size+1))]
                        l_ennntl = r_nnn_brtl[j-(row_size+1)]
                        x_ennnbl = state[2*(j+(row_size-1))]
                        l_ennnbl = r_nnn_bltr[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_etl*g[gi]*x_etl + l_ebl*g[gi]*x_ebl + l_nnnu*g[gi]*x_nnnu + l_ennntl*g[gi]*x_ennntl + l_ennnbl*g[gi]*x_ennnbl- (l_l+l_etl+l_ebl+l_nnnu+l_ennntl+l_ennnbl)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)


                #odd rows
                else: #(j // row_size) % 2 == 0
                
                    #1st column
                    if j % row_size == 0:
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_obr = state[2*(j+row_size)]
                        l_obr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_onnntr = state[2*(j-(row_size-1))]
                        l_onnntr = r_nnn_bltr[j-(row_size-1)]
                        x_onnnbr = state[2*(j+(row_size+1))]
                        l_onnnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_r*g[gi]*x_r + l_otr*g[gi]*x_otr + l_obr*g[gi]*x_obr + l_nnnu*g[gi]*x_nnnu + l_onnntr*g[gi]*x_onnntr + l_onnnbr*g[gi]*x_onnnbr - (l_r+l_otr+l_obr+l_nnnu+l_onnntr+l_onnnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)

                        
                
                    #2nd column
                    elif j % row_size == 1:
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_otl = state[2*(j-(row_size+1))]
                        l_otl = r_brtl[j-(row_size+1)]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_obl = state[2*(j+(row_size-1))]
                        l_obl = r_bltr[j]
                        x_obr = state[2*(j+row_size)]
                        l_obr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_onnntr = state[2*(j-(row_size-1))]
                        l_onnntr = r_nnn_bltr[j-(row_size-1)]
                        x_onnnbr = state[2*(j+(row_size+1))]
                        l_onnnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_otl*g[gi]*x_otl + l_otr*g[gi]*x_otr + l_obl+g[gi]*x_obl +l_obr*g[gi]*x_obr + l_nnnu*g[gi]*x_nnnu + l_onnntr*g[gi]*x_onnntr + l_onnnbr*g[gi]*x_onnnbr - (l_l+l_r+l_otl+l_otr+l_obl+l_obr+l_nnnu+l_onnntr+l_onnnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)


                    #middle columns
                    elif j % row_size != 0 and j % row_size != 1 and j % row_size != (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_otl = state[2*(j-(row_size+1))]
                        l_otl = r_brtl[j-(row_size+1)]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_obl = state[2*(j+(row_size-1))]
                        l_obl = r_bltr[j]
                        x_obr = state[2*(j+row_size)]
                        l_obr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_onnntl = state[2*(j-(row_size+2))]
                        l_onnntl = r_nnn_brtl[j-(row_size+2)]
                        x_onnntr = state[2*(j-(row_size-1))]
                        l_onnntr = r_nnn_bltr[j-(row_size-1)]
                        x_onnnbl = state[2*(j+(row_size-2))]
                        l_onnnbl = r_nnn_bltr[j]
                        x_onnnbr = state[2*(j+(row_size+1))]
                        l_onnnbr = r_nnn_brtl[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_otl*g[gi]*x_otl + l_otr*g[gi]*x_otr + l_obl+g[gi]*x_obl +l_obr*g[gi]*x_obr + l_nnnu*g[gi]*x_nnnu + l_onnntl*g[gi]*x_onnntl + l_onnntr*g[gi]*x_onnntr + l_onnnbl*g[gi]*x_onnnbl + l_onnnbr*g[gi]*x_onnnbr - (l_l+l_r+l_otl+l_otr+l_obl+l_obr+l_nnnu+l_onnntl+l_onnntr+l_onnnbl+l_onnnbr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_otl = state[2*(j-(row_size+1))]
                        l_otl = r_brtl[j-(row_size+1)]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_obl = state[2*(j+(row_size-1))]
                        l_obl = r_bltr[j]
                        x_obr = state[2*(j+row_size)]
                        l_obr = r_brtl[j]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_onnntl = state[2*(j-(row_size+2))]
                        l_onnntl = r_nnn_brtl[j-(row_size+2)]
                        x_onnnbl = state[2*(j+(row_size-2))]
                        l_onnnbl = r_nnn_bltr[j]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_otl*g[gi]*x_otl + l_otr*g[gi]*x_otr + l_obl+g[gi]*x_obl +l_obr*g[gi]*x_obr + l_nnnu*g[gi]*x_nnnu + l_onnntl*g[gi]*x_onnntl + l_onnnbl*g[gi]*x_onnnbl - (l_l+l_otl+l_otr+l_obl+l_obr+l_nnnu+l_onnntl+l_onnnbl)*g[gi]*x
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
                        l_r = r_lr[j]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_etr = state[2*(j-(row_size-1))]
                        l_etr = r_bltr[j-(row_size-1)]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_ennntr = state[2*(j-(row_size-2))]
                        l_ennntr = r_nnn_bltr[j-(row_size-2)]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_nnnu*g[gi]*x_nnnu + l_ennntr*g[gi]*x_ennntr - (l_r+l_etl+l_etr+l_nnnu+l_ennntr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)

                    #middle columns
                    elif j % row_size != 0 and j % row_size != (row_size-2) and j % row_size != (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_etr = state[2*(j-(row_size-1))]
                        l_etr = r_bltr[j-(row_size-1)]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_ennntl = state[2*(j-(row_size+1))]
                        l_ennntl = r_nnn_brtl[j-(row_size+1)]
                        x_ennntr = state[2*(j-(row_size-2))]
                        l_ennntr = r_nnn_bltr[j-(row_size-2)]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_nnnu*g[gi]*x_nnnu + l_ennntl*g[gi]*x_ennntl + l_ennntr*g[gi]*x_ennntr - (l_l+l_r+l_etl+l_etr+l_nnnu+l_ennntl+l_ennntr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)

                    #2nd to last column
                    elif j % row_size == (row_size-2):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_etr = state[2*(j-(row_size-1))]
                        l_etr = r_bltr[j-(row_size-1)]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_ennntl = state[2*(j-(row_size+1))]
                        l_ennntl = r_nnn_brtl[j-(row_size+1)]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_etl*g[gi]*x_etl + l_etr*g[gi]*x_etr + l_nnnu*g[gi]*x_nnnu + l_ennntl*g[gi]*x_ennntl - (l_l+l_r+l_etl+l_etr+l_nnnu+l_ennntl)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                        
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_etl = state[2*(j-row_size)]
                        l_etl = r_brtl[j-row_size]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_ennntl = state[2*(j-(row_size+1))]
                        l_ennntl = r_nnn_brtl[j-(row_size+1)]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_etl*g[gi]*x_etl + l_nnnu*g[gi]*x_nnnu + l_ennntl*g[gi]*x_ennntl - (l_l+l_etl+l_nnnu+l_ennntl)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)


                    
                #odd rows
                else: #(j // row_size) % 2 == 0
                
                    #1st column
                    if j % row_size == 0:
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_onnntr = state[2*(j-(row_size-1))]
                        l_onnntr = r_nnn_bltr[j-(row_size-1)]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_r*g[gi]*x_r + l_otr*g[gi]*x_otr + l_nnnu*g[gi]*x_nnnu + l_onnntr*g[gi]*x_onnntr - (l_r+l_otr+l_nnnu+l_onnntr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)

                        
                
                    #2nd column
                    elif j % row_size == 1:
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_otl = state[2*(j-(row_size+1))]
                        l_otl = r_brtl[j-(row_size+1)]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_onnntr = state[2*(j-(row_size-1))]
                        l_onnntr = r_nnn_bltr[j-(row_size-1)]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_otl*g[gi]*x_otl + l_otr*g[gi]*x_otr + l_nnnu*g[gi]*x_nnnu + l_onnntr*g[gi]*x_onnntr - (l_l+l_r+l_otl+l_otr+l_nnnu+l_onnntr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)


                    #middle columns
                    elif j % row_size != 0 and j % row_size != 1 and j % row_size != (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_r = state[2*j + 2]
                        l_r = r_lr[j]
                        x_otl = state[2*(j-(row_size+1))]
                        l_otl = r_brtl[j-(row_size+1)]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_onnntl = state[2*(j-(row_size+2))]
                        l_onnntl = r_nnn_brtl[j-(row_size+2)]
                        x_onnntr = state[2*(j-(row_size-1))]
                        l_onnntr = r_nnn_bltr[j-(row_size-1)]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_r*g[gi]*x_r + l_otl*g[gi]*x_otl + l_otr*g[gi]*x_otr + l_nnnu*g[gi]*x_nnnu + l_onnntl*g[gi]*x_onnntl + l_onnntr*g[gi]*x_onnntr - (l_l+l_r+l_otl+l_otr+l_nnnu+l_onnntl+l_onnntr)*g[gi]*x
                        ODEs = np.append(ODEs, dx1dt)
                        dy1dt = d*x - e*y
                        ODEs = np.append(ODEs, dy1dt)
                        
                
                    #last column
                    elif j % row_size == (row_size-1):
                        x_l = state[2*j - 2]
                        l_l = r_lr[j-1]
                        x_otl = state[2*(j-(row_size+1))]
                        l_otl = r_brtl[j-(row_size+1)]
                        x_otr = state[2*(j-row_size)]
                        l_otr = r_bltr[j-row_size]
                        x_nnnu = state[2*(j-2*row_size)]
                        l_nnnu = r_nnn_du[j-2*row_size]
                        x_onnntl = state[2*(j-(row_size+2))]
                        l_onnntl = r_nnn_brtl[j-(row_size+2)]
                        a = (N*(h[j] + theta + (1/3)))/(theta + 1)
                        dx1dt = a - x + b*x**2 - c*x**3 - f*x*y + l_l*g[gi]*x_l + l_otl*g[gi]*x_otl + l_otr*g[gi]*x_otr + l_nnnu*g[gi]*x_nnnu + l_onnntl*g[gi]*x_onnntl - (l_l+l_otl+l_otr+l_nnnu+l_onnntl)*g[gi]*x
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
        number_of_peaks_lists[list_i] = np.append(number_of_peaks_lists[list_i], np.sum(number_of_peaks)/number_of_cells)


        total_neighbor_pairs = sum(r_nnn_bltr) + sum(r_nnn_brtl) + sum(r_nnn_du) + sum(r_bltr) + sum(r_brtl) + sum(r_lr)

            
        #Counting the number of causal edges between nearest-neighbors
        #j is the link index in the link matrices, (the total number of 1's represents the total neighbor pairs)
        #j is also the cell's index in tmax_list
        number_of_edges = np.array([])

        lower_time = .5
        upper_time = 2.5

        #r_lr matrix
        j = 0
        while j < np.size(r_lr):
            if j % row_size != (row_size-1):
                if r_lr[j] == 1:
                    if tmax_list[j] != 0 and tmax_list[j+1] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+1]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+1]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                    else:
                        number_of_edges = np.append(number_of_edges, 0)
            j += 1
        
        #r_nnn_du
        j = 0
        while j < np.size(r_nnn_du):
            if j % row_size != (row_size-1):
                if r_nnn_du[j] == 1:
                    if tmax_list[j] != 0 and tmax_list[j+(2*row_size)] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+(2*row_size)]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+(2*row_size)]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                    else:
                        number_of_edges = np.append(number_of_edges, 0)
            j += 1

        #r_bltr
        j = 0
        while j < np.size(r_bltr):
            if j % (2*row_size) != 0:
                if r_bltr[j] == 1:
                
                    #even rows
                    if (j // row_size) % 2 != 0:
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)

                    else: #(j // row_size) % 2 == 0
                        if tmax_list[j] != 0 and tmax_list[j+(row_size-1)] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+(row_size-1)]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+(row_size-1)]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
            
            j += 1
        
        #r_brtl
        j = 0
        while j < np.size(r_brtl):
            if j % (2*row_size) != (2*row_size - 1):
                if r_brtl[j] == 1:
                
                    #even rows
                    if (j // row_size) % 2 != 0:
                        if tmax_list[j] != 0 and tmax_list[j+(row_size+1)] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+(row_size+1)]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+(row_size+1)]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)

                    else: #(j // row_size) % 2 == 0
                        if tmax_list[j] != 0 and tmax_list[j+row_size] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+row_size]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+row_size]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
            
            j += 1

        #r_nnn_brtl
        j = 0
        while j < np.size(r_nnn_brtl):
            if j % row_size != (row_size - 1) and j % (2*row_size) != (2*row_size - 2):
                if r_nnn_brtl[j] == 1:
                
                    #even rows
                    if (j // row_size) % 2 != 0:
                        if tmax_list[j] != 0 and tmax_list[j+(row_size+2)] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+(row_size+2)]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+(row_size+2)]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)

                    else: #(j // row_size) % 2 == 0
                        if tmax_list[j] != 0 and tmax_list[j+(row_size+1)] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+(row_size+1)]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+(row_size+1)]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
            
            j += 1
        

        #r_nnn_bltr
        j = 0
        while j < np.size(r_nnn_bltr):
            if j % row_size != (row_size - 1) and j % (2*row_size) != (2*row_size - 2):
                if r_nnn_bltr[j] == 1:
                
                    #even rows
                    if (j // row_size) % 2 != 0:
                        if tmax_list[j] != 0 and tmax_list[j+(row_size-1)] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+(row_size-1)]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+(row_size-1)]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)

                    else: #(j // row_size) % 2 == 0
                        if tmax_list[j] != 0 and tmax_list[j+(row_size-2)] != 0:
                            if np.absolute(tmax_list[j] - tmax_list[j+(row_size-2)]) >= lower_time and np.absolute(tmax_list[j] - tmax_list[j+(row_size-2)]) <= upper_time:
                                number_of_edges = np.append(number_of_edges, 1)
                            else:
                                number_of_edges = np.append(number_of_edges, 0)
                            
                        else:
                            number_of_edges = np.append(number_of_edges, 0)
            
            j += 1


######################################################################################
        
        #The number of edges for each trial is added to number_of_edges_list
        number_of_edges_lists[list_i] = np.append(number_of_edges_lists[list_i], np.sum(number_of_edges)/total_neighbor_pairs)    
        #print(np.sum(number_of_edges))

        list_i += 1
        
        II += 10*I_total
        
    #print(i)   
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



#Saving data
np.savez("324_cells_Pedge_Period_30_30_disordered_Results", edges1 = edges1, edges2 = edges2, edges3 = edges3, edges4 = edges4, edges5 = edges5, edges6 = edges6, edges7 = edges7, edges8 = edges8, peaks1 = peaks1, peaks2 = peaks2, peaks3 = peaks3, peaks4 = peaks4, peaks5 = peaks5, peaks6 = peaks6, peaks7 = peaks7, peaks8 = peaks8)




