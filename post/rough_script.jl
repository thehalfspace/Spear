#################################
# SOME ROUGH SCRIPTS
# FOR TRYING OUT STUFF
#################################

using StatsBase
using PyPlot

# Get index for each event
function event_indx(tStart, tEnd, time_)
    start_indx = zeros(size(tStart))
    end_indx = zeros(size(tEnd))

    for i = 1:length(tStart)
        
        start_indx[i] = findall(time_ .== tStart[i])[1]
        end_indx[i] = findall(time_ .== tEnd[i])[1]
    end

    return start_indx, end_indx
end

# Plot sliprates for each event with depth
function test1(S, O, evno)
    start_indx = zeros(size(O.tStart))
    end_indx = zeros(size(O.tEnd))

    for i = 1:length(O.tStart)
        
        start_indx[i] = findall(O.time_ .== O.tStart[i])[1]
        end_indx[i] = findall(O.time_ .== O.tEnd[i])[1]
    end
    
    start_indx = Int.(start_indx)[evno]
    end_indx = Int.(end_indx)[evno]
    sv = zeros(size(O.seismic_slipvel))
    
    inc = 0.1       # time interval = 0.1 sec for plotting
    to = O.time_[start_indx]    # start time
    j = 1
    for i=start_indx:end_indx
        if O.time_[i] >= to
            sv[:,j] = O.seismic_slipvel[:,i]
            to = to + inc
            j = j+1
        end
    end
    
    sv = sv[:,1:j]

    fig = PyPlot.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    
    ax.plot(sv, S.FltX/1e3, ".--", label="a", lw = 1)
    ax.set_xlabel("Slip rate (m/s)")
    ax.set_ylabel("Depth (km)")
    ax.set_title("Slip rate for one event")
    #  ax.set_xlim([0, 0.02])
    ax.set_ylim([-24, 0])
    show()
    
    figname = string(path, "slipvel.pdf")
    fig.savefig(figname, dpi = 300)
end

