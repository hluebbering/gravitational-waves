# Generate a synthetic gravitational wave signal for a binary black hole merger
# 1. compute waveform for inspiral portion of merger
# 2. compute merger/ringdown portion 
# 3. match them together

# -----------------------------------------


# Import libraries 
import numpy as np
import riroriro.inspiralfuns as ins


# -----------------------------------------


# Inspiral portion (inspiralfuns) 

# GW150914-like
logMc = 1.4 # log(Mc/Msun) between 0.0 and 2.0 
q = 0.8 # q between 0.1 and 1.0

# defaults
flow = 10.0 # (Hz)
merger_type = 'BH' # Binary black hole merger 
                   # Use 'NS' for binary/black hole neutron star mergers
D = 100.0 # (Mpc)

# Obtain total mass (M) and symmetric mass ratio (eta)
M, eta = ins.get_M_and_eta(logMc = logMc, q = q)

# Calculate two parameters used by subsequent integration
start_x = ins.startx(M, flow)
end_x = ins.endx(eta, merger_type)

# Integrate post-Newtonian parameter, 
x, xtimes, dt = ins.PN_parameter_integration(start_x, end_x, M,eta)

# Convert timesteps (xtimes) from geometric units into seconds
realtimes = ins.inspiral_time_conversion(xtimes, M)

# Calculate phase, angular, and linear frequency of inspiralling binary
i_phase, omega, freq = ins.inspiral_phase_freq_integration(x, dt, M)

# Calculate radius and its time-derivative
r, rdot = ins.radius_calculation(x, M, eta)

# Calculate values over time of two polarisations waveforms
# (orthogonal/plus and diagonal/cross) 
A1, A2 = ins.a1_a2_calculation(r, rdot, omega, D, M, eta)
i_Aorth, i_Adiag = ins.inspiral_strain_polarisations(A1, A2, i_phase) # strain polarisations

# Compute overall strain amplitude of signal over duration of inspiral
i_amp = ins.inspiral_strain_amplitude(i_Aorth, i_Adiag)

# Reduce size of lists to reduce file size (optional) 
i_time = ins.list_size_reducer(100, realtimes)
i_omega = ins.list_size_reducer(100, omega)
i_phase = ins.list_size_reducer(100, i_phase)
i_amp = ins.list_size_reducer(100, i_amp)
i_Aorth = ins.list_size_reducer(100, i_Aorth)
i_Adiag = ins.list_size_reducer(100, i_Adiag)


# -----------------------------------------


# Merger/ringdown portion (mergerfirstfuns)
import riroriro.mergerfirstfuns as me1

# Constant coefficients, depending on symmetric mass ratio
sfin, wqnm = me1.quasi_normal_modes(eta)
alpha, b, C, kappa = me1.gIRS_coefficients(eta, sfin)

# Calculate frequencies for merger/ringdown portion
fhat, m_omega = me1.merger_freq_calculation(wqnm, b, C, kappa)
fhatdot = me1.fhat_differentiation(fhat)

# Convert timesteps from geometric units to seconds
m_time = me1.merger_time_conversion(M)


# -----------------------------------------


# Matching inspiral and merger/ringdown (matchingfuns) into single waveform
import riroriro.matchingfuns as mat

# Calculate timestep indices for which switch from inspiral to merger/ringdown should occur 
# to give smoothest transition in frequency and frequency gradient 
min_switch_ind = mat.min_switch_ind_finder(i_time,i_omega,m_time,m_omega)
final_i_index = mat.final_i_index_finder(min_switch_ind,i_omega,m_omega)

# Calculate offset applied to time for merger/ringdown portion 
# so that two portions are aligned in time at transition
time_offset = mat.time_offset_finder(min_switch_ind,final_i_index,i_time,m_time)

# Create combined lists for time and angular frequency values for entire binary merger
i_m_time, i_m_omega = mat.time_frequency_stitching(min_switch_ind,final_i_index,time_offset,i_time,i_omega,m_time,m_omega)

# Convert angular frequency of merging binary to a linear frequency 
# which is what would be observed by gravitational wave detectors
i_m_freq = mat.frequency_SI_units(i_m_omega, M)


# -----------------------------------------

# Rest of merger/ringdown functions (mergersecondfuns)
import riroriro.mergersecondfuns as me2

# Calculate values of phase for merger/ringdown portion 
m_phase = me2.merger_phase_calculation(min_switch_ind,final_i_index,i_phase,m_omega)

# Combine inspiral and merger/ringdown phase lists
# phase values are smooth across transition - computed after matching so we know when transition is
i_m_phase = me2.phase_stitching(final_i_index, i_phase, m_phase) 

# Calculate strain amplitude for merger/ringdown portion 
m_amp = me2.merger_strain_amplitude(min_switch_ind,final_i_index,alpha,i_amp,m_omega,fhat,fhatdot)

# Combine amplitude lists for two portions of waveform
i_m_amp = me2.amplitude_stitching(final_i_index,i_amp,m_amp)



# -----------------------------------------


# Graphs show values of strain (envelope) amplitude & frequency over time for whole waveform
import matplotlib.pyplot as plt

plt.figure(1)
plt.plot(i_m_time,i_m_amp)
plt.xlabel('Time (s)')
plt.ylabel('Strain amplitude')

plt.figure(2)
plt.plot(i_m_time,i_m_freq)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')

# Calculate strain polarisations for merger/ringdown portion 
m_Aorth, m_Adiag = me2.merger_polarisations(final_i_index,m_amp,m_phase,i_Aorth)

# Combine polarisation lists for two portions of waveform
i_m_Aorth, i_m_Adiag = me2.polarisation_stitching(final_i_index,i_Aorth,i_Adiag,m_Aorth,m_Adiag)

# Plot evolution of strain polarisations at end of inspiral and full merger/ringdown
plt.figure(1)
plt.plot(i_m_time,i_m_Aorth)
plt.plot(i_m_time,i_m_Adiag)
plt.axis([6.1,6.25,-6.2e-21,6.2e-21])
plt.xlabel('Time (s)')
plt.ylabel('Strain amplitude')


