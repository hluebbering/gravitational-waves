
"""
Parts of the procedure for simulating the inspiral portions of gravitational
waves, collected into modular functions.
"""
import numpy as np


def get_M_and_eta(**kwargs):
    """Obtain total mass and symmetric mass ratio
    Parameters:
        kwargs (floats): all of the inputs are floats
    Returns:
        M (float): Total mass
        eta (float): Symmetric mass ratio
    """    
    # making sure that all of the inputs are floats
    for each_variable in kwargs.values():
        assert type(each_variable) == float, 'All inputs should be floats.'

    #making sure that the user doesn't put arguments from both methods at once
    if ('m1' in kwargs or 'm2' in kwargs) and ('logMc' in kwargs or 'q' in kwargs):
        raise TypeError('Please don\'t mix the two input methods.')

    #making sure that both of the arguments for a method are specified
    if (('m1' in kwargs) != ('m2' in kwargs)) or (('logMc' in kwargs) != ('q' in kwargs)):
        raise TypeError('Please specify both arguments required by a method.')

    if 'q' in kwargs and kwargs['q'] > 1:  # making sure that q is not above 1
        raise ValueError('The mass ratio needs to be less than or equal to 1, '
                         'i.e. the smaller mass divided by the larger.')

    if 'm1' in kwargs:  # first method
        M = kwargs['m1'] + kwargs['m2']
        eta = (kwargs['m1']*kwargs['m2'])/M**2
    elif 'q' in kwargs:  # second method
        eta = kwargs['q'] / ((1+kwargs['q'])**2)
        M = 10**kwargs['logMc'] / (eta**(3/5))

    return (M, eta)


def startx(M, flow):
    """Calculate one of two parameters used by subsequent integration
    Parameters:
        M (float): Total Mass
        flow (float): _description_
    Returns:
        value (float): Value based on Buskirk eq. 22
    """
    #input type checking
    assert type(M) == float, 'M should be a float.'
    assert type(flow) == float, 'flow should be a float.'

    #basic constants
    pi = np.pi
    c = 2.99792458e8
    #geometric unit conversion
    Msunkg = 1.9891e30

    value = (6.673e-11 * M*Msunkg*pi*flow*c**-3)**(2/3)  # Buskirk eq. 22

    return value


def endx(eta, merger_type):
    """Calculate one of two parameters used by subsequent integration
    Parameters:
        eta (float): Symmetric mass ratio
        merger_type (string): Binary/black hole neutron star mergers
    Returns:
        value (float): Value based on Buskirk eq. 23
    """
    #input type checking
    assert type(eta) == float, 'eta should be a float.'

    if merger_type == 'BH':
        value = (1/3)*(1 + (7/18)*eta) # Buskirk eq. 23
    elif merger_type == 'NS':
        value = (1/6)*(1 + (7/18)*eta)  # Buskirk eq. 23
    else:
        raise ValueError('merger_type must be either \'BH\' for BH-BH or '
                         '\'NS\' for BH-NS or NS-NS.')
    return value


def PNderiv(x, M, eta):
    """Calculate the post-Newtonian Derivative
    Parameters:
        x (float): Post-Newtonian parameter
        M (float): Total mass
        eta (float): Symmetric mass ratio
    Returns:
        Mdxdt (float): Post-Newtonian parameter derivative
    """
    pi = np.pi
    # coefficients from Huerta article
    a0 = 12.8*eta
    a2 = eta*(-16*(924*eta + 743))/420
    a3 = 51.2*pi*eta
    a4 = (1903104*eta**2 + 3934368*eta + 1091296)*(eta/45360)
    # Buskirk also includes a8 to a12 (6PN)
    a8 = 170.799 - 742.551*eta + 370.173*eta**2 - 43.4703*eta**3 - \
        0.0249486*eta**4 + (14.143 - 150.692*eta)*np.log(x)
    a9 = 1047.25 - 2280.56*eta + 923.756 * \
        eta**2 + 22.7462*eta**3 - 102.446*np.log(x)
    a10 = 714.739 - 1936.48*eta + 3058.95*eta**2 - 514.288*eta**3 + \
        29.5523*eta**4 - 0.185941*eta**5 + \
        (-3.00846 + 1019.71*eta + 1146.13*eta**2)*np.log(x)
    a11 = 3622.99 - 11498.7*eta + 12973.5*eta**2 - 1623*eta**3 + \
        25.5499*eta**4 + (83.1435 - 1893.65*eta)*np.log(x)
    a12 = 11583.1 - 45878.3*eta + 33371.8*eta**2 - 7650.04*eta**3 + \
        648.748*eta**4 - 14.5589*eta**5 - 0.0925075*eta**6 + \
        (-1155.61 + 7001.79*eta - 2135.6*eta**2 -
         2411.92*eta**3)*np.log(x) + 33.2307*np.log(x**2)
    Mdxdt = a0*x**5 + a2*x**6 + a3*x**6.5 + a4*x**7 + ((1/M)*(64/5)*eta*x**5
                                                       * (1 + a8*x**4 + a9*x**4.5 + a10*x**5 + a11*x**5.5 + a12*x**6))
    return Mdxdt


def PN_parameter_integration(start, end, M, eta):
    """Integrate post-Newtonian parameter
    Parameters:
        start (float): Value based on Buskirk eq. 22
        end (float): Value based on Buskirk eq. 23
        M (float): Total mass
        eta (float): Symmetric mass ratio
    Returns:
        x (array): Post-Newtonian parameter
        xtimes (array): Timesteps in geometric units
        dt (array): Holds the timestep
    """
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == float, 'All inputs should be floats.'

    dt = []  # holds the timestep
    xtimes = [0]  # times when x is evaluated
    i = 0  # loop counter
    # Euler integration with variable timestep
    x = [start]  # post-Newtonian parameter
    while x[i] <= end:
        dt.append((10**-6*x[i])/PNderiv(x[i], M, eta)) # timestep with 10^-6 threshold
        x.append(x[i] + 10**-6*x[i])  # Euler step
        xtimes.append(xtimes[i] + dt[i])  # increment time-storer
        i += 1  # increment counter

    return [x, xtimes, dt]


def inspiral_time_conversion(xtimes, M):
    """Convert timesteps (xtimes) from geometric units into seconds
    Parameters:
        xtimes (list): Timesteps in geometric units
        M (float): Total mass
    Returns:
        realtimes (list): Timesteps in seconds
    """
    #input type checking
    assert type(xtimes) == list, 'xtimes should be a list.'
    assert type(M) == float, 'M should be a float.'

    #geometric unit conversion
    Msuns = 4.923e-6

    realtimes = np.empty((len(xtimes)))  # initialisation of list
    for i in range(len(xtimes)):
        realtimes[i] = xtimes[i]*M*Msuns

    #output type conversion
    realtimes = list(realtimes)

    return realtimes


def inspiral_phase_freq_integration(x, dt, M):
    """Calculate phase, angular, and linear frequency of inspiralling binary
    Parameters:
        x (list): Post-Newtonian parameter
        dt (list): Holds the timestep
        M (float): Total mass
    Returns:
        i_phase (list): Orbital phase
        omega (list): Angular frequency
        freq (list): Frequency
    """
    #input type checking
    assert type(x) == list, 'x should be a list.'
    assert type(dt) == list, 'dt should be a list.'
    assert type(M) == float, 'M should be a float.'

    pi = np.pi
    # geometric unit conversion
    Msuns = 4.923e-6

    i_phase = np.empty((len(x)))  # orbital phase
    omega = np.empty((len(x)))  # angular frequency
    freq = np.empty((len(x)))  # frequency

    # initial values
    i_phase[0] = 0
    omega[0] = x[0]**1.5
    freq[0] = omega[0]/(M*Msuns)

    # for the phase integration, we are tied to the times at which we evaluated x
    for i in range(len(x)-1):  # -1 because always i+1 terms
        omega[i+1] = x[i+1]**1.5
        i_phase[i+1] = i_phase[i] + omega[i+1]*dt[i]
        freq[i+1] = omega[i+1] / (M*Msuns*pi)

    # output type conversion
    i_phase = list(i_phase)
    omega = list(omega)
    freq = list(freq)

    return [i_phase, omega, freq]


def radius_calculation(x, M, eta):
    """Calculate radius and its time-derivative
    Parameters:
        x (list): Post-Newtonian parameter
        M (float): Total mass
        eta (float): Symmetric mass ratio
    Returns:
        r (list): Orbital radius
        rdot (list): Derivative of radius
    """
    # input type checking
    assert type(x) == list, 'x should be a list.'
    assert type(M) == float, 'M should be a float.'
    assert type(eta) == float, 'eta should be a float.'

    # post-Newtonian correction coefficients from Buskirk
    r0pn = 1
    r1pn = -1 + 0.333333*eta
    r2pn = 4.75*eta + 0.111111*eta**2
    r3pn = -7.51822*eta - 3.08333*eta**2 + 0.0246914*eta**3

    r = np.empty((len(x)))  # orbital radius (geometric u.)
    rdot = np.empty((len(x)))  # derivative of radius

    for i in range(len(x)):  # Buskirk radius equations
        r[i] = r0pn*(1/x[i]) + r1pn + r2pn*x[i] + r3pn*x[i]**2
        rdot[i] = PNderiv(x[i], M, eta) * (-2*r0pn*x[i]
                                           ** -2 + r2pn + 2*r3pn*x[i])

    #output type conversion
    r = list(r)
    rdot = list(rdot)

    return [r, rdot]


def a1_a2_calculation(r, rdot, omega, D, M, eta):
    """Calculate values over time of two polarisations waveforms
    Parameters:
        r (list): Orbital radius
        rdot (list): Derivative of radius
        omega (list): Angular frequency
        D (float): Mpc
        M (float): Total mass
        eta (float): Symmetric mass ratio
    Returns:
        A1 (list): Values of first polarisation wave
        A2 (list): Values of second polarisation wave
    """
    #input type checking
    assert type(r) == list, 'r should be a list.'
    assert type(rdot) == list, 'rdot should be a list.'
    assert type(omega) == list, 'omega should be a list.'
    assert type(D) == float, 'D should be a float.'
    assert type(M) == float, 'M should be a float.'
    assert type(eta) == float, 'eta should be a float.'

    Dkm = D * 3.086e19  # conversion from Mpc to km
    A1 = np.empty((len(r)))
    A2 = np.empty((len(r)))

    for i in range(len(r)):  # based on Buskirk eq. 9
        A1[i] = (-2*M*eta*(1/Dkm))*(rdot[i]**2 + (r[i]*omega[i])**2 + 1/r[i])
        A2[i] = (-2*M*eta*(1/Dkm))*(2*r[i]*rdot[i]*omega[i])

    # output type conversion
    A1 = list(A1)
    A2 = list(A2)

    return [A1, A2]


def inspiral_strain_polarisations(A1, A2, i_phase):
    """Calculate strain polarisations
    Parameters:
        A1 (list): Values of first polarisation wave
        A2 (list): Values of second polarisation wave
        i_phase (list): Orbital phase
    Returns:
        Aorth (list): Orthogonal/plus polarisation
        Adiag (list): Diagonal/cross polarisation
    """
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == list, 'All inputs should be lists.'

    Aorth = np.empty((len(i_phase)))  # orthogonal/plus polarisation
    Adiag = np.empty((len(i_phase)))  # diagonal/cross polarisation

    for i in range(len(i_phase)):
        Aorth[i] = A1[i]*np.cos(2*i_phase[i]) + A2[i]*np.sin(2*i_phase[i])
        Adiag[i] = A1[i]*np.sin(2*i_phase[i]) - A2[i]*np.cos(2*i_phase[i])

    #output type conversion
    Aorth = list(Aorth)
    Adiag = list(Adiag)

    return [Aorth, Adiag]


def inspiral_strain_amplitude(Aorth, Adiag):
    """Compute overall strain amplitude of signal over duration of inspiral
    Parameters:
        Aorth (list): Orthogonal/plus polarisation
        Adiag (list): Diagonal/cross polarisation
    Returns:
        i_amp (list): Strain amplitude of signal over inspiral portion
    """
    #input type checking
    assert type(Aorth) == list, 'Aorth should be a list.'
    assert type(Adiag) == list, 'Adiag should be a list.'

    i_amp = np.empty((len(Aorth)))
    for i in range(len(Aorth)):
        i_amp[i] = np.sqrt(Aorth[i]**2 + Adiag[i]**2)

    #output type conversion
    i_amp = list(i_amp)

    return i_amp


def list_size_reducer(reduction_factor, your_list):
    """_summary_
    Parameters:
        reduction_factor (int): _description_
        your_list (list): _description_
    Returns:
        reduced_list (_type_): _description_
    """
    # input type checking
    assert type(reduction_factor) == int, 'reduction_factor should be an int.'
    assert type(
        your_list) == list, ('The thing to be reduced needs to be a  list.')

    # create new list with every nth point of your_list
    reduced_list = [your_list[0]]
    for i in range(reduction_factor, len(your_list), reduction_factor):
        reduced_list.append(your_list[i])

    return reduced_list






"""
The pre-matching parts of the procedure for simulating the merger/ringdown
portions of gravitational waves from binary black holes, collected into modular
functions.
"""

def quasi_normal_modes(eta):
    #input type checking
    assert type(eta) == float, 'eta should be a float.'
    
    assert eta <= 0.25 and eta > 0, ('eta should be positive and no larger '
                                     'than 0.25.')
    
    sfin = 2*np.sqrt(3)*eta - (390/79)*eta**2 + (2379/287)*eta**3 - \
        (4621/276)*eta**4               #final spin; based on Buskirk eq. 20
    sfin = float(sfin)                  #becomes numpy.float64 otherwise
    wqnm = 1 - 0.63*(1 - sfin)**0.3     #based on Buskirk eq. 19
    
    return (sfin,wqnm)

def gIRS_coefficients(eta,sfin):
    #input type checking
    assert type(eta) == float, 'eta should be a float.'
    assert type(sfin) == float, 'sfin should be a float.'
    
    #constants from Buskirk Appendix C
    Q = 2/((1 - sfin)**0.45)
    alpha = Q**-2 * (16313/562 + (21345/124)*eta)
    b = 16014/979 - (29132/1343)*eta**2
    C = 206/903 + (180/1141)*np.sqrt(eta) + (424/1205)*eta**2*(1/np.log(eta))
    kappa = 713/1056 - (23/193)*eta
    
    #output type conversion
    C = float(C)
    
    return (alpha,b,C,kappa)

def merger_freq_calculation(wqnm,b,C,kappa):
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == float, 'All inputs should be floats.'
    
    #time in geometric units    
    time = np.empty((201))    #can't set it directly as range because of typing
    for i in range(201):
        time[i] = range(-100,101)[i]
        
    fhat = np.empty((201))
    m_omega = np.empty((201))
    for i in range(201):
        fhat[i] = (C/2) * (1 + 1/kappa)**(1 + kappa) * (1 - (1 + \
            (1/kappa)*np.exp(-2*time[i]*(1/b)))**-kappa)
        #based on Buskirk eq. 17
        m_omega[i] = 0.5 * wqnm * (1 - fhat[i])
        #based on Buskirk eq. 18, with 0.5 to standardise weird angle
        #definitions in Buskirk paper
        
    #output type conversion
    fhat = list(fhat)
    m_omega = list(m_omega)
        
    return [fhat,m_omega]

def fhat_differentiation(fhat):
    #input type checking
    assert type(fhat) == list, 'fhat should be a list.'
    
    #making sure fhat has the correct length, to avoid IndexErrors
    assert len(fhat) == 201, 'fhat should have length 201.'
    
    fhatdot = np.empty((201))
    #unusual averages at ends
    fhatdot[0] = fhat[1] - fhat[0]
    fhatdot[200] = fhat[200] - fhat[199]
    for i in range(1,200):              #derivative approximated as differences
        fhatdot[i] = 0.5*(fhat[i+1] - fhat[i-1])
        
    #output type conversion
    fhatdot = list(fhatdot)
    
    return fhatdot

def merger_time_conversion(M):
    #input type checking
    assert type(M) == float, 'M should be a float.'
    
    #geometric unit conversion
    Msuns=4.923e-6
    
    time = range(-100,101)                          #time in geometric units
    m_time = np.empty((201))                        #time in s (initialisation)
    for i in range(201):
        m_time[i] = time[i]*M*Msuns
        
    #output type conversion
    m_time = list(m_time)
    
    return m_time






"""
Parts of the procedure for matching together the inspiral and merger/ringdown
portions of the gravitational waveforms of binary black holes, collected into
modular functions.
"""
def MQdiff(i,i_time,i_omega,m_time,m_omega):
    #input type checking
    assert type(i) == int, 'i should be an int.'
    assert type(i_time) == list, 'i_time should be a list.'
    assert type(i_omega) == list, 'i_omega should be a list.'
    assert type(m_time) == list, 'm_time should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    
    #making sure i is in the range of indices used for the merger/ringdown
    assert 0 <= i <= 200, 'i should be in range(201).'
    
    #with this MQ method, we iterate for each frequency value in the merger
    #part, look for the closest frequency in the inspiral part, measure the df
    #difference between these points in the two waveforms, and then select the
    #frequency with the minimum df difference as the switching point, adjusting
    #the waveform times appropriately so the combined waveform is continuous in
    #f and t
    
    try:
        closest_index = np.searchsorted(i_omega, m_omega[i], side='right')
        #with this method, we use searchsorted in the frequency domain instead
        #of the time domain
        #it does assume frequency increases monotonously, which it should
        df_diff = abs((m_omega[i] - m_omega[i-1])/(m_time[i] - m_time[i-1]) - \
                      (i_omega[closest_index] - i_omega[closest_index - 1])/ \
                      (i_time[closest_index] - i_time[closest_index - 1]))
    except IndexError:
        df_diff = np.nan
        #to get rid of known searchsorted errors for some irrelevant
        #frequencies
        
    return df_diff

def min_switch_ind_finder(i_time,i_omega,m_time,m_omega):
    #input type checking
    for each_variable in locals().values():
        assert type(each_variable) == list, 'All inputs should be lists.'
        
    minMQ = 1000                                    #deliberately overly high
    MQ = []                                         #holds MQ values
    min_switch_ind = []
    
    #(method description is in MQdiff() code)
    
    for i in range(1,len(m_time)):
        #MQ calculation at each point in merger waveform
        #EXCEPT 0, because then the [i-1] in MQdiff rolls over to the last
        #value in the array (and the first index shouldn't be the switch point
        #anyway)
        MQ = MQdiff(i,i_time,i_omega,m_time,m_omega)
        if MQ < minMQ:
            minMQ = MQ                              #update minimum
            min_switch_ind = i                      #time of new minimum
            
    return min_switch_ind

def final_i_index_finder(min_switch_ind,i_omega,m_omega):
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(i_omega) == list, 'i_omega should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    
    final_i_index = np.searchsorted(i_omega, m_omega[min_switch_ind], \
                                    side='right')
    
    #output type conversion
    final_i_index = int(final_i_index)    
    
    return final_i_index

def time_offset_finder(min_switch_ind,final_i_index,i_time,m_time):
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_time) == list, 'i_time should be a list.'
    assert type(m_time) == list, 'm_time should be a list.'
    
    time_offset = i_time[final_i_index] - m_time[min_switch_ind]
    
    #output type conversion
    time_offset = float(time_offset)
    
    return time_offset

def time_frequency_stitching(min_switch_ind,final_i_index,time_offset,i_time,\
                             i_omega,m_time,m_omega):
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(time_offset) == float, 'time_offset should be a float.'
    assert type(i_time) == list, 'i_time should be a list.'
    assert type(i_omega) == list, 'i_omega should be a list.'
    assert type(m_time) == list, 'm_time should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    
    min_offset_m_time = np.empty((len(m_time)))
    for i in range(len(m_time)):                    #offsetting to match i_time
        min_offset_m_time[i] = m_time[i] + time_offset
        
    #now we stitch the inspiral and merger frequency waveforms together    
    i_m_omega = []                                  #combined omega
    i_m_time = []                                   #combined time
    for i in range(final_i_index):                  #inspiral segment
        i_m_omega.append(i_omega[i])
        i_m_time.append(i_time[i])
        
    for i in range(min_switch_ind,len(m_time)):     #merger segment
        i_m_omega.append(m_omega[i])
        i_m_time.append(min_offset_m_time[i])
        
    return [i_m_time,i_m_omega]

def frequency_SI_units(i_m_omega,M):
    #input type checking
    assert type(i_m_omega) == list, 'i_m_omega should be a list.'
    assert type(M) == float, 'M should be a float.'
    
    pi=np.pi
    Msuns=4.923e-6                                  #geometric unit conversion
    
    i_m_freq = np.empty((len(i_m_omega)))
    for i in range(len(i_m_omega)):
        i_m_freq[i] = i_m_omega[i] / (M*Msuns*pi)
        
    #output type conversion
    i_m_freq = list(i_m_freq)
        
    return i_m_freq




"""
The post-matching parts of the procedure for simulating the merger/ringdown
portions of gravitational waves from binary black holes, collected into modular
functions.
"""
def merger_phase_calculation(min_switch_ind,final_i_index,i_phase,m_omega):
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_phase) == list, 'i_phase should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    
    m_phase = np.empty((201 - min_switch_ind))
    #starts at min_switch_ind instead of 0 because the part of the
    #merger/ringdown data before that gets cut off by the matching
    m_phase[0] = i_phase[final_i_index]
    #matching phase at start of merger/ringdown to phase at end of inspiral
    for i in range(min_switch_ind + 1,201):
        m_phase[i - min_switch_ind] = m_phase[i - min_switch_ind - 1] + \
            m_omega[i]
                                        #Euler integration Buskirk of eq. 21
    
    #output type conversion
    m_phase = list(m_phase)
    
    return m_phase

def phase_stitching(final_i_index,i_phase,m_phase):
    #input type checking
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_phase) == list, 'i_phase should be a list.'
    assert type(m_phase) == list, 'm_phase should be a list.'
    
    i_m_phase = np.concatenate((i_phase[:final_i_index],m_phase))
    
    #output type conversion
    i_m_phase = list(i_m_phase)
    
    return i_m_phase

def merger_strain_amplitude(min_switch_ind,final_i_index,alpha,i_amp,m_omega,\
                            fhat,fhatdot):
    #input type checking
    assert type(min_switch_ind) == int, 'min_switch_ind should be an int.'
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(alpha) == float, 'alpha should be a float.'
    assert type(i_amp) == list, 'i_amp should be a list.'
    assert type(m_omega) == list, 'm_omega should be a list.'
    assert type(fhat) == list, 'fhat should be a list.'
    assert type(fhatdot) == list, 'fhatdot should be a list.'
    
    m_amp = np.empty((201 - min_switch_ind))
    for i in range(min_switch_ind,201):         #initial unscaled calculation
        m_amp[i - min_switch_ind] = (1/(2e0*m_omega[i])) * ((abs(fhatdot[i]))/\
             (1 + alpha*(fhat[i]**2 - fhat[i]**4)))**(1/2)
                                                #Buskirk eq. 16 (note 2e0)           
        
    matching_amplitude = i_amp[final_i_index]   #amplitude at end of inspiral
    scaling_ratio = matching_amplitude / m_amp[0]
    for i in range(len(m_amp)):
        m_amp[i] = scaling_ratio * m_amp[i]
        #rescaling for continuity between inspiral and merger/ringdown
        
    #output type conversion
    m_amp = list(m_amp)
    
    return m_amp

def amplitude_stitching(final_i_index,i_amp,m_amp):
    #input type checking
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_amp) == list, 'i_amp should be a list.'
    assert type(m_amp) == list, 'm_amp should be a list.'
    
    i_m_amp = np.concatenate((i_amp[:final_i_index],m_amp))
    
    #output type conversion
    i_m_amp = list(i_m_amp)
    
    return i_m_amp

def merger_polarisations(final_i_index,m_amp,m_phase,i_Aorth):
    #input type checking
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(m_amp) == list, 'm_amp should be a list.'
    assert type(m_phase) == list, 'm_phase should be a list.'
    assert type(i_Aorth) == list, 'i_Aorth should be a list.' 
    
    m_Aorth = np.empty((len(m_amp)))            #orthogonal/plus polarisation
    m_Adiag = np.empty((len(m_amp)))            #diagonal/cross polarisation
    
    #adjusting phase to ensure continuity between inspiral and merger
    #polarisations at the switching point
    comparison_phase_1 = (1/2) * np.arccos(i_Aorth[final_i_index] / m_amp[0])
    comparison_phase_2 = np.pi - comparison_phase_1
    sign_i_Aorth = np.sign(i_Aorth[final_i_index] - i_Aorth[final_i_index-1])
    sign_comparison = np.sign(np.cos(2*(m_phase[1] - m_phase[0] + \
        comparison_phase_1)) - np.cos(2*comparison_phase_1))
    if sign_i_Aorth == sign_comparison: #signs match, phase_1 is correct
        comparison_phase = comparison_phase_1
    else:                               #signs don't match, phase_2 is correct
        comparison_phase = comparison_phase_2
    
    phase_difference = m_phase[0] - comparison_phase
    adjusted_m_phase = np.empty((len(m_phase)))
    for i in range(len(m_amp)):    #adjusting phase for polarisation continuity
        adjusted_m_phase[i] = m_phase[i] - phase_difference
    
    for i in range(len(m_amp)):
        m_Aorth[i] = m_amp[i] * np.cos(2*adjusted_m_phase[i])
        m_Adiag[i] = m_amp[i] * np.sin(2*adjusted_m_phase[i])
    
    #output type conversion
    m_Aorth = list(m_Aorth)
    m_Adiag = list(m_Adiag)
    
    return [m_Aorth,m_Adiag]

def polarisation_stitching(final_i_index,i_Aorth,i_Adiag,m_Aorth,m_Adiag):
    #input type checking
    assert type(final_i_index) == int, 'final_i_index should be an int.'
    assert type(i_Aorth) == list, 'i_Aorth should be a list.'
    assert type(i_Adiag) == list, 'i_Adiag should be a list.'
    assert type(m_Aorth) == list, 'm_Aorth should be a list.'
    assert type(m_Adiag) == list, 'm_Adiag should be a list.'
    
    i_m_Aorth = np.concatenate((i_Aorth[:final_i_index],m_Aorth))
    i_m_Adiag = np.concatenate((i_Adiag[:final_i_index],m_Adiag))
    
    #output type conversion
    i_m_Aorth = list(i_m_Aorth)
    i_m_Adiag = list(i_m_Adiag)
    
    return [i_m_Aorth,i_m_Adiag]