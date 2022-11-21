using DelimitedFiles

include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")

# path to save files
global path = "$(@__DIR__)/plots/test_01/"
mkpath(path)

global out_path = "$(@__DIR__)/data/test_01/"


#= If any of these files don't exist, or you want more information, check lines 153-158 in src/main.jl and uncomment.
Also remember to uncomment the corresponding `end` statements for each of these in lines 418-427 in src/main.jl.
=#

# Global variables
yr2sec = 365*24*60*60

# Read data
event_time = readdlm(string(out_path, "event_time.out"), header=false)
tStart = event_time[:,1]        # List of start times for all earthquakes (in seconds)
tEnd = event_time[:,2]          # List of end times for all earthquakes (in seconds)
hypo = event_time[:,3]          # List of hypocenters for all earthquakes (in meters depth)

event_stress = readdlm(string(out_path, "event_stress.out"), header=false)
indx = Int(length(event_stress[1,:])/2)
taubefore = event_stress[:,1:indx]      # List of shear stress along depth before each earthquake
tauafter = event_stress[:,indx+1:end]   # List of shear stress along depth after each earthquake

delfafter = readdlm(string(out_path, "coseismic_slip.out"), header=false)       # Coseismic slip = (slip_after - slip_before) each earthquake along depth


# These are slip and sliprate stored at every certain timestep (e.g. every 10 or 100 timesteps)
# These are very big files, so I usually only store them when necessary
#  slip = readdlm(string(out_path, "slip.out"), header=false)                   
#  sliprate = readdlm(string(out_path, "sliprate.out"), header=false)

# Order of storage: Seff, tauo, FltX, cca, ccb, xLf
params = readdlm(string(out_path, "params.out"), header=false)

Seff = params[1,:]      # Effective normal stress
tauo = params[2,:]      # Initial shear stress
FltX = params[3,:]      # Fault location along depth
cca = params[4,:]       # Rate-state-friction parameter 'a'
ccb = params[5,:]       # Rate-state-friction parameter 'b'
Lc = params[6,:]        # Rate-state-friction parameter 'Lc'

# Index (location) of fault from 0 to 18 km
flt18k = findall(FltX .<= 18)[1]

time_vel = readdlm(string(out_path, "time_velocity.out"), header=false)
t = time_vel[:,1]               # timesteps
Vfmax = time_vel[:,2]           # Max. slip rate
Vsurface = time_vel[:,3]        # Surface slip rate


rho1 = 2670
vs1 = 3464
rho2 = 2500
vs2 = 0.6*vs1
mu = rho2*vs2^2

delfsec = readdlm(string(out_path, "delfsec.out"))              # slip along depth of fault plotted every 0.1 second (edit line 84 in src/main.jl to change)
delfyr = readdlm(string(out_path, "delfyr.out"))                # slip along depth plotted every 2 years (edit line 81 to change)
#stress = readdlm(string(out_path, "stress.out"), header=false)  # this is shear stress along depth and time: again, very large file so I typically don't save this 

#start_index = get_index(stress', taubefore')
stressdrops = taubefore .- tauafter

# Mw = moment magnitude
# del_sigma = stress drop
# fault_slip = average slip for one rupture
# rupture_len = length of the rupture

Mw, del_sigma, fault_slip, rupture_len =
        moment_magnitude_new(mu, FltX, delfafter', stressdrops', t);

