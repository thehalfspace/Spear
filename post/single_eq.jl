########################################
#   ISOLATE EACH EARTHQUAKE STUFF
########################################

function stress_evol(seismic_stress, FltX)

    #  x_axis = seismic_stress[1:5:end,:]

    fig = PyPlot.figure(figsize=(4,4))
    ax = fig.add_subplot(111)

    for i in 1:length(seismic_stress[1,:])
        ax.plot(seismic_stress[:,i].*1e-1 .+ 1e-1*i, -FltX./1e3, "k-")
    end
    
    #  ax.plot(seismic_stress, -FltX./1e3, "k.-")
    
    ax.xaxis.set_ticklabels([])
    ax.set_xlabel("Sliprate (m/s)")
    ax.set_ylabel("Depth (km)")
    #  ax.set_title("Shear Stress evolut")
    ax.set_ylim([0, 15])
    ax.invert_yaxis()
    fig.tight_layout()
    show()

    figname = string(path, "ss_evol.png")
    fig.savefig(figname, dpi = 300)

end

# process zone
function process(stress1, slip1, stress2, slip2)

    #  x_axis = seismic_stress[1:5:end,:]

    fig = PyPlot.figure(figsize=(4,4))
    ax = fig.add_subplot(111)

    ax.plot(slip1, stress1, "b-", label="Model 2")
    ax.plot(slip2, stress2, "r-", label="Model 3")
    ax.set_xlabel("Slip (m)")
    ax.set_ylabel("Shear Stress (MPa)")
    #  ax.set_title("Shear Stress evolut")
    #  ax.set_ylim([0, 15])
    #  ax.invert_yaxis()
    ax.legend(loc="upper right")
    show()

    figname = string(path, "process.png")
    fig.savefig(figname, dpi = 300)

end

# Proposal figure stuff
function fig1(delfsec, delf5yr, FltX, Mw, hypo)
    indx = findall(abs.(FltX) .<= 18e3)[1]

    delfsec2 = delfsec[indx:end, :]

    fig = PyPlot.figure(figsize=(12,5))

    ax1 = plt.subplot2grid((1,6), (0,0), colspan=3)
    #  ax2 = plt.subplot2grid((1,6), (0,3), colspan=1)
    ax3 = plt.subplot2grid((1,6), (0,4), colspan=2)

    #  ax1 = fig.add_subplot(131)
    #  ax2 = fig.add_subplot(132)
    #  ax3 = fig.add_subplot(133)
    
    # Shade the fault zone region
    x_shade = LinRange(5,25,25)
    y1 = repeat([0],25)
    y3 = repeat([24],25)

    ax1.plot(delf5yr, -FltX/1e3, color="royalblue", lw=1, alpha=1.0)
    ax1.plot(delfsec2, -FltX[indx:end]/1e3, "-", color="chocolate", lw=1, alpha=1.0)
    ax1.fill_between(x_shade, y3, y1, color="chocolate", alpha=0.2)
    #  ax1.set_xlabel("Accumulated Slip (m)")
    ax1.set_ylabel("Depth (km)")
    ax1.set_title("Cumulative Slip History")
    ax1.set_ylim([0,24])
    ax1.set_xlim([5,20])
    ax1.invert_yaxis()

    # hypocenter
    #  hist = fit(Histogram, -hypo./1e3, closed=:right, nbins=10)

    #  y1 = repeat([0],10)
    #  y3 = repeat([24],10)
    #  ax2.barh(hist.edges[1][1:end-1], hist.weights)
    #  #  ax2.fill_between(LinRange(0,10,10),y3,y1, color="tab:blue", alpha=0.3)
    #  ax2.set_xlabel("Number of Earthquakes")
    #  #  ax2.set_ylabel("Depth (km)")
    #  ax2.set_title("Hypocenter Location")
    #  ax2.set_ylim([0, 24])
    #  ax2.invert_yaxis()
    #  fig.tight_layout()

    # mfd
    hist2 = fit(Histogram, Mw, nbins = 9)

    # Cumulative
    cum = cumsum(hist2.weights[end:-1:1])[end:-1:1]

    ax3.plot(hist2.edges[1][1:end-1], cum, "k.", markersize=10) 
    ax3.set_xlabel("Moment Magnitude (Mw)")
    ax3.set_ylabel("Number of Earthquakes")
    ax3.set_yscale("log")
    ax3.set_title("Magnitude-frequency distribution")
    #  ax[:set_xlim]([2, 7])
    ax3.set_ylim([1, 30])


    show()
    
    figname = string(path, "proposal_fig.png")
    fig.savefig(figname, dpi = 300)

end
