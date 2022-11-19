using DelimitedFiles

include("$(@__DIR__)/post/event_details.jl")
include("$(@__DIR__)/post/plotting_script.jl")

# path to save files
global path = "$(@__DIR__)/plots/test_01/"
mkpath(path)

global out_path = "$(@__DIR__)/data/test_01/"

# Global variables
yr2sec = 365*24*60*60

# Read data
# event_time = readdlm(string(out_path, "event_time.out"), header=false)
# tStart = event_time[:,1]
# tEnd = event_time[:,2]
# hypo = event_time[:,3]

# event_stress = readdlm(string(out_path, "event_stress.out"), header=false)
# indx = Int(length(event_stress[1,:])/2)
# taubefore = event_stress[:,1:indx]
# tauafter = event_stress[:,indx+1:end]

# delfafter = readdlm(string(out_path, "coseismic_slip.out"), header=false)
#  slip = readdlm(string(out_path, "slip.out"), header=false)
# sliprate = readdlm(string(out_path, "sliprate.out"), header=false)

# Read SCEC data
fltst000 = readdlm(string(out_path, "fltst_dp000.out"), header=false)
time_ = fltst000[:,1]
slip = fltst000[:,2]
sliprate = fltst000[:,3]
stress = fltst000[:,4]

# Comparison data
comp_scec = readdlm(string(out_path, "../scec-data/abdelmeguid_fltst_dp_000.txt"), header=false, skipstart=13)
comp_time = comp_scec[:,1]
comp_sliprate = comp_scec[:,3]


# Order of storage: Seff, tauo, FltX, cca, ccb, xLf
params = readdlm(string(out_path, "params.out"), header=false)

Seff = params[1,:]
tauo = params[2,:]
FltX = params[3,:]
cca = params[4,:]
ccb = params[5,:]
Lc = params[6,:];

# Index of fault from 0 to 18 km
# flt18k = findall(FltX .<= 18)[1]

# time_vel = readdlm(string(out_path, "time_velocity.out"), header=false)
# t = time_vel[:,1]
# Vfmax = time_vel[:,2]
# Vsurface = time_vel[:,3]
# alphaa = time_vel[:,4]

# delfsec = readdlm(string(out_path, "delfsec.out"))
# delfyr = readdlm(string(out_path, "delfyr.out"))
# stress = readdlm(string(out_path, "stress.out"), header=false)

#  start_index = get_index(seismic_stress, taubefore)
#  stressdrops = taubefore .- tauafter

#  Mw, del_sigma, fault_slip, rupture_len =
        #  moment_magnitude_new(mu, P1, P3.FltX, delfafter, stressdrops, time_);
