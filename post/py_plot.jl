#########################################################
##
##  Plotting templates using Matplotlib and PyPlot
##
#########################################################

using PyPlot
using StatsBase
using LaTeXStrings
using PyCall
mpl = pyimport("matplotlib")


# Matplotlib global tweaks
#  mpl.rcParams["axes.linewidth"] = 2.0

# Plot Vfmax
function VfmaxPlot(Vfmax, time_, yr2sec)

    #  Vfmax = maximum(SlipVel, dims = 1)[:]
    
    
    fig = PyPlot.figure(figsize=(12,5))
    ax = fig.add_subplot(111)
    
    plt.rc("font",size=24)
    plt.rc("axes", linewidth=2.0)

    ax.plot(time_./yr2sec, Vfmax, lw = 2.0)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    ax.get_yaxis().set_tick_params(which="both", direction="in")
    ax.get_xaxis().set_tick_params(which="both", direction="in")
    #  plt.rc("grid", linestyle="--", color="black", alpha=0.5)
    #  plt.grid("True")

    ax.set_xlim([230,400])
    
    plt.tight_layout()
    show()
    
    figname = string(path, "Vfmax01.png")
    fig.savefig(figname, dpi = 300)
end

# Plot friction parameters
function icsPlot(a_b, Seff, tauo, FltX)
    fig = PyPlot.figure(figsize=(8,10))
    ax = fig.add_subplot(111)
    
    plt.rc("font",size=16)
    plt.rc("axes", linewidth=2.0)

    ax.plot(Seff, -FltX./1e3, lw = 2.0, label="Normal Stress")
    ax.plot(tauo, -FltX./1e3, lw = 2.0, label="Shear Stress")
    ax.set_xlabel("Stresses (Mpa)")
    ax.set_ylabel("Depth (km)")
    
    ax2 = ax.twiny()
    ax2.plot(a_b, -FltX./1e3, lw = 2.0, label="(a-b)")
    ax2.set_xlabel("Rate-state friction value")
    
    ax.set_ylim([0,16])
    ax.invert_yaxis()
    ax.get_yaxis().set_tick_params(which="both", direction="in")
    ax.get_xaxis().set_tick_params(which="both", direction="in")
    #  plt.rc("grid", linestyle="--", color="black", alpha=0.5)
    #  plt.grid("True")

    
    #  plt.tight_layout()
    show()
    
    figname = string(path, "ics_02.png")
    fig.savefig(figname, dpi = 300)
end

# Plot slip vs event number
function slipPlot(delfafter, rupture_len, FltX, Mw, tStart)

    fig = PyPlot.figure(figsize=(12,8))
    ax1 = fig.add_subplot(121)

    xaxis = tStart[Mw .>2.8]   #collect(1:length(delfafter[1,:]))

    Mw2 = Mw[Mw .> 2.8]

    plt.rc("font",size=14)
    
    # Normalize colorbar
    norm = matplotlib.colors.Normalize(vmin = minimum(Mw2), vmax=maximum(Mw2)) 
    colors = matplotlib.cm.inferno_r(norm(Mw2))

    ax1.barh(xaxis, delfafter[end-1,:], height=10, 
             label="At trench: 60 m depth", color=colors, align="center"); 
    ax1.set_xlabel("Coseismic Slip (m) at 60 m depth")
    #ax.set_xlabel("Coseismic slip (m)")
    ax1.set_ylabel("Time (yr)")
    ax1.get_yaxis().set_tick_params(which="both", direction="in")
    ax1.get_xaxis().set_tick_params(which="both", direction="in")
    #  ax.set_yticks(collect(0:50:maximum(xaxis)))
    #  ax.set_xlim([0, 2])

    trench_depth = findall(abs.(FltX) .< 6.0e3)[1]

    ax2 = fig.add_subplot(122)
    
    ax2.barh(xaxis, delfafter[trench_depth,:], height=10, 
             label="At trench: 6 km depth", color=colors, align="center"); 
    
    #  ax2.set_xlim([0, 2])
    ax2.set_xlabel("Coseismic Slip (m) at 6 km depth")
    
    #  ax.set_yticks(round.(xaxis, digits = 0))
    #  plt.legend()
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap="inferno_r")
    sm.set_array([])
    fig.colorbar(sm, shrink=0.9, label="Mw")

    show()
    
    figname = string(path, "mw2_time.png")
    fig.savefig(figname, dpi = 300)
end

# Plot Alphaa
function alphaaPlot(alphaa, time_, yr2sec)
    
    fig = PyPlot.figure(figsize=(12,4))
    ax = fig.add_subplot(111)
    
    plt.rc("font",size=16)

    ax.plot(time_./yr2sec, alphaa, lw = 2)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Shear Modulus Contrast (%)")
    ax.get_yaxis().set_tick_params(which="both", direction="in")
    ax.get_xaxis().set_tick_params(which="both", direction="in")
    #  plt.rc("grid", linestyle="--", color="black", alpha=0.5)
    #  plt.grid("True")
    
    plt.tight_layout()
    show()


    figname = string(path, "alphaa_01.png")
    fig.savefig(figname, dpi = 300)
end

# Plot recurrence
function recurrencePlot(tStart, Mw, yr2sec)

    #  Vfmax = maximum(SlipVel, dims = 1)[:]
    
    rec2 = diff(tStart)./yr2sec
    rec = rec2[rec2 .> 0.1]
    xa = collect(1:length(rec))
    Mw = Mw[2:end]

    Mw2 = Mw[rec2 .> 0.1]

    fig = PyPlot.figure(figsize=(8,6), dpi = 300)
    ax = fig.add_subplot(111)
    
    plt.rc("font",size=12)

    ax.scatter(xa, rec, s=Mw.^2)
    ax.set_xlabel("Interevent Number")
    ax.set_ylabel("Reccurence Interval (yr)")
    #  ax.set_title("Max. slip rate on fault")
    #  ax.get_yaxis().set_tick_params(which="both", direction="in")
    #  ax.get_xaxis().set_tick_params(which="both", direction="in")
    #  plt.rc("grid", linestyle="--", color="black", alpha=0.5)
    #  plt.grid("True")
    
    show()


    figname = string(path, "recurrence.png")
    fig.savefig(figname, dpi = 300)
end

# Plot cumulative slip
function cumSlipPlot(delfsec, delf5yr, FltX)
    indx = findall(abs.(FltX) .<= 18e3)[1]

    delfsec2 = delfsec[indx:end, :]

    fig = PyPlot.figure(figsize=(10,7), dpi=300)
    ax = fig.add_subplot(111)
    plt.rc("font",size=12)

    ax.plot(delf5yr, -FltX/1e3, color="royalblue", lw=0.5, alpha=1.0)
    ax.plot(delfsec2, -FltX[indx:end]/1e3, "-", color="chocolate", lw=0.5, alpha=1.0)
    ax.set_xlabel("Accumulated Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,24])
    #  ax.set_xlim([1,20])
    
    ax.invert_yaxis()
    #  ax.get_yaxis().set_tick_params(which="both", direction="in")
    #  ax.get_xaxis().set_tick_params(which="both", direction="in")
    #  plt.rc("grid", linestyle="--", color="black", alpha=0.5)
    
    show()
    
    figname = string(path, "cumslip_02.png")
    fig.savefig(figname, dpi = 300)

end

function csPlot2(vfsec, delf5yr, start_index, FltX)
    indx = findall(abs.(FltX) .<= 18e3)[1]
    vfsec = log10.(vfsec[indx:end, :])

    #  s_end = size(delfsec)[2]
    #  seis_indx = push!(start_index,s_end)
    
    y = -FltX./1e3

    Xsec = LinRange(0, maximum(vfsec), length(vfsec[1,:]))

    xgrid_sec,ygrid_sec = Matlab.meshgrid(Xsec, y)

    fig = PyPlot.figure(figsize=(10,7), dpi=300)
    ax = fig.add_subplot(111)
    plt.rc("font",size=12)

    cbar = ax.contourf(xgrid_sec, ygrid_sec, vfsec)
    #  ax.contourf(xgrid_yr, ygrid_yr, delf5yr)

    ax.set_xlabel("Timesteps (Variable)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,24])
    #  ax.set_xlim([1,20])
    
    ax.invert_yaxis()
    #  ax.get_yaxis().set_tick_params(which="both", direction="in")
    #  ax.get_xaxis().set_tick_params(which="both", direction="in")
    #  plt.rc("grid", linestyle="--", color="black", alpha=0.5)
    plt.colorbar(cbar)
    show()
    
    figname = string(path, "cs2.png")
    fig.savefig(figname, dpi = 300)

end

function MwPlot(Mw)

    hist = fit(Histogram, Mw, nbins = 15)

    # Cumulative
    cum = cumsum(hist.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(8,7))
    ax = fig.add_subplot(111)
    plt.rc("font",size=12)

    #  ax.plot](hist.edges[1][1:end-1], hist.weights, ".", label="Non-cumulative")
    ax.plot(hist.edges[1][1:end-1], cum, ".", markersize=20) #, label="Cumulative")
    ax.set_xlabel("Moment Magnitude (Mw)")
    ax.set_ylabel("Number of Earthquakes")
    ax.set_yscale("log")
    show()

    figname = string(path, "mfd.png")
    fig.savefig(figname, dpi = 300)
end


## Plot paleoseismic record (Depth vs event no)
function paleoRecordPlot(FltX, Mw, delfafter, tStart, yr2sec)
    """
        Parameters should be scaled with Mw
    """

    depth = -FltX./1e3
    coseismic_slip = delfafter  #[:,Mw .> 4.6] 
    evno = length(Mw[Mw .> 4.6])
    
    x_axis = tStart./yr2sec

    sampling_depth = [0, 2, 4, 6, 8, 10]
    sampling_index = Int.(zeros(length(sampling_depth)))

    for i=1:length(sampling_depth)
        sampling_index[i] = findall(depth .> sampling_depth[i])[end]
    end

    fig = PyPlot.figure(figsize=(12,2), dpi = 300)
    ax = fig.add_subplot(111)
    
    plt.rc("font",size=4)

    for i = 1:length(sampling_depth)
        _size = coseismic_slip[sampling_index[i],:]
        ax.scatter(x_axis, sampling_depth[i]*ones(evno), color="k", s = 20*_size.^2)
    end
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Depth (km)")
    ax.set_xticks(round.(x_axis, digits = 0))
    ax.invert_yaxis()
    ax.get_yaxis().set_tick_params(which="both", direction="in")
    ax.get_xaxis().set_tick_params(which="both", direction="in")
    plt.rc("grid", linestyle="--", color="black", alpha=0.3)
    plt.grid("True", axis="x")
    
    show()
    figname = string(path, "paleo_record.png")
    fig.savefig(figname, dpi = 300)
end

# spatiotemporal imshow
function sptempPlot(seismic_slipvel2, FltX)
    
    indx = findall(abs.(FltX) .<= 18e3)[1]
    value = seismic_slipvel2[indx:end,:]
    
    depth = -FltX[indx:end]./1e3


    fig = PyPlot.figure(figsize=(7,6), dpi=100)
    plt.rc("font",size=12)
    ax = fig.add_subplot(111)

    c = ax.imshow(value, cmap="viridis", aspect="auto", 
                  norm=matplotlib.colors.LogNorm(vmin=1e-2, vmax=1e0), 
                         interpolation="bicubic",
                         extent=[0,length(seismic_slipvel2[1,:])/10, 0,18])

    #   ax.set_yticks(ax.get_yticks()[1:2:end])
    #   ax.set_xticks(ax.get_yticks()[1:2:end])

    #   ax.get_yaxis().set_tick_params(which="both", direction="in")
    #   ax.get_xaxis().set_tick_params(which="both", direction="in")
    
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Depth (km)")

    ax.invert_yaxis()
    cbar = fig.colorbar(c)
    #   cbar.set_ticks(cbar.get_ticks()[1:2:end])

    
    show()
    figname = string(path, "spatiotemporal.png")
    fig.savefig(figname, dpi = 300)
end


# Spatiotemporal Slip Rate evolution
function spatiotemporalPlot(seismic_sliprate, FltX)
    
    indx = findall(abs.(FltX) .<= 18e3)[1]
    value = seismic_slipvel[indx:end,:]
    
    depth = -FltX[indx:end]./1e3

    x = range(0, length(value[1,:])/1, length=length(value[1,:]))

    fig = PyPlot.figure(figsize=(7,6), dpi=100)
    ax = fig.add_subplot(111)

    c = ax.pcolormesh(x, depth[end:-1:1], value, vmin=0, vmax=1, 
                      cmap="inferno")
    ax.set_yticks(ax.get_yticks()[1:2:end])
    ax.set_xticks(ax.get_yticks()[1:2:end])

    ax.get_yaxis().set_tick_params(which="both", direction="in")
    ax.get_xaxis().set_tick_params(which="both", direction="in")

    ax.invert_yaxis()
    cbar = fig.colorbar(c)
    cbar.set_ticks(cbar.get_ticks()[1:2:end])
    
    show()
    figname = string(path, "spatiotemporal.png")
    fig.savefig(figname, dpi = 300)


end
