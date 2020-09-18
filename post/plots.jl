#################################
#   PLOTS FOR EARTHQUAKE CYCLES
#################################

using PyPlot

# Customize plot
#  PyPlot.matplotlib[:rc]("font", size = 24)
#  PyPlot.matplotlib[:rc]("axes", labelsize = 21)
#  PyPlot.matplotlib[:rc]("axes", titlesize = 21)
#PyPlot.matplotlib[:rc]("xtick", labelsize = 12)
#PyPlot.matplotlib[:rc]("ytick", labelsize = 12)
#  PyPlot.matplotlib[:rc]("legend", fontsize = 24)
#  PyPlot.matplotlib[:rc]("figure", titlesize = 24)

# Plot cumulative slip
function cumSlipPlot(delfsec, delf5yr, FltX)
    
    FZ = -8

    indx = findall(abs.(FltX) .<= 18e3)[1]

    delfsec2 = delfsec[indx:end, :]

    fig = PyPlot.figure(figsize=(4,4))
    ax = fig.add_subplot(111)
    
    # Shade the fault zone region
    x_shade = LinRange(5,25,25)
    y1 = repeat([0],25)
    y2 = repeat([FZ],25)
    y3 = repeat([-24],25)

    ax.plot(delf5yr, -FltX/1e3, color="royalblue", lw=1, alpha=1.0)
    ax.plot(delfsec2, -FltX[indx:end]/1e3, "-", color="chocolate", lw=1, alpha=1.0)
    #  ax[:fill_between](x_shade, y2, y1, color="chocolate", alpha=0.3)
    #  ax[:fill_between](x_shade, y3, y1, color="chocolate", alpha=0.3)
    ax.set_xlabel("Accumulated Slip (m)")
    ax.set_ylabel("Depth (km)")
    #  ax.set_title("Cumulative Slip History")
    ax.set_ylim([0,24])
    #  ax.set_xlim([0.5,5.5])
    ax.invert_yaxis()
    show()
    
    figname = string(path, "cumslip_devol.png")
    fig.savefig(figname, dpi = 300)

end

function plotHypo(hypo)  #S, Slip, SlipVel, Stress, time_)

    # Plot hypocenter
    hist = fit(Histogram, -hypo./1e3, closed=:right, nbins=10)

    fig = PyPlot.figure(figsize=(5,7))
    ax = fig.add_subplot(111)

    ax.barh(hist.edges[1][1:end-1], hist.weights)
    #  ax[:plot](collect(1:80), -8*ones(80), "k--", label="Fault Zone Depth")
    ax.set_xlabel("Number of Earthquakes")
    ax.set_ylabel("Depth (km)")
    ax.set_title("Hypocenter Location")
    ax.set_ylim([0, 24])
    ax.invert_yaxis()
    fig.tight_layout()
    show()

    figname = string(path, "hypo.png")
    fig.savefig(figname, dpi = 300)

end


#...........
# Plot MFD
#...........
function MwPlot(Mw)

    hist = fit(Histogram, Mw, nbins = 9)

    # Cumulative
    cum = cumsum(hist.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(8,7))
    ax = fig.add_subplot(111)

    #  ax.plot](hist.edges[1][1:end-1], hist.weights, ".", label="Non-cumulative")
    ax.plot(hist.edges[1][1:end-1], cum, "k.--", markersize=20) #, label="Cumulative")
    ax.set_xlabel("Moment Magnitude (Mw)")
    ax.set_ylabel("Number of Earthquakes")
    ax.set_yscale("log")
    #  ax.set_title("Magnitude-frequency distribution")
    #  ax.set_xlim([2, 7])
    ax.set_ylim([1, 200])
    #  ax.legend(loc="upper right")
    show()

    figname = string(path, "mfd.png")
    fig.savefig(figname, dpi = 300)
end

function MwPlot2(Mw1, Mw2, Mw3, Mw4, Mw5, Mw6)

    hist1 = fit(Histogram, Mw1, nbins = 6)
    hist2 = fit(Histogram, Mw2, nbins = 6)
    hist3 = fit(Histogram, Mw3, nbins = 6)
    hist4 = fit(Histogram, Mw4, nbins = 6)
    hist5 = fit(Histogram, Mw5, nbins = 6)
    hist6 = fit(Histogram, Mw5, nbins = 6)

    # Cumulative
    cum1 = cumsum(hist1.weights[end:-1:1])[end:-1:1]
    cum2 = cumsum(hist2.weights[end:-1:1])[end:-1:1]
    cum3 = cumsum(hist3.weights[end:-1:1])[end:-1:1]
    cum4 = cumsum(hist4.weights[end:-1:1])[end:-1:1]
    cum5 = cumsum(hist5.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(12,9))
    ax = fig.add_subplot(111)

    ax.plot(hist1.edges[1][1:end-1], cum1, ".-", markersize=10, label="400 m")
    ax.plot(hist2.edges[1][1:end-1], cum2, ".-", markersize=10, label="320 m")
    ax.plot(hist3.edges[1][1:end-1], cum3, ".-", markersize=10, label="250 m")
    ax.plot(hist4.edges[1][1:end-1], cum4, ".-", markersize=10, label="160 m")
    ax.plot(hist5.edges[1][1:end-1], cum5, ".-", markersize=10, label="100 m")
    ax.set_xlabel("Moment Magnitude (Mw)")
    ax.set_ylabel("Number of Earthquakes")
    ax.set_yscale("log")
    ax.set_title("Magnitude-frequency distribution")
    #  ax.set_xlim([1.7, 7])
    ax.set_ylim([1e-1, 1e3])
    ax.legend(loc="upper right")
    show()

    figname = string(path, "mfdcomp.png")
    fig.savefig(figname, dpi = 300)
end

#.................................
# Plot earthquake catalog
# (Earthquake magnitude with time)
#.................................
function eq_catalog(Mw, t_catalog, yr2sec)

    fig = PyPlot.figure(figsize=(12,9))
    ax = fig.add_subplot(111)

    ax.scatter(t_catalog./yr2sec, Mw, s= 30, marker=".")
    ax.set_xlabel("Time (yrs)")
    ax.set_ylabel("Moment Magnitude (Mw)")
    ax.set_title("Earthquake Catalogue")
    show()

    figname = string(path, "catalogue.pdf")
    fig.savefig(figname, dpi = 300)
end



# Plot friction parameters
function fricPlot(cca, ccb, FltX)
    
    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig.add_subplot(111)

    ax.plot(cca, FltX/1e3, "ko--", label="a", lw = 1)
    ax.plot(ccb, FltX/1e3, "ko--", label="b", lw = 1)
    ax.plot(cca-ccb, FltX/1e3, label="a-b", lw = 1)
    ax.set_xlabel("Value")
    ax.set_ylabel("Depth (km)")
    ax.set_title("Rate and State friction parameters")
    ax.legend(loc="upper right")
    ax.set_ylim([-24, 0])
    show()

    #  figname = string(dir, "/plots", name, "/fric.pdf")
    figname = string(path, "friction.pdf")
    fig.savefig(figname, dpi = 300)

end

# Plot shear stress at location as a function of time
function stressPlot(Stress, time_, FltX, yr2sec, loc1 = 8e3)
    
    FltID = findall(abs.(FltX) .<= loc1)[1]
    
    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig.add_subplot(111)

    ax.plot(time_/yr2sec, Stress[FltID, :], lw = 1)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Shear stress (MPa)")
    ax.set_title(string("Shear stress at ", loc1/1e3, "km depth"))   
    show()

    #  figname = string(dir, "/plots", name, "/shear.pdf")
    #  fig[:savefig](figname, dpi = 300)

end


# Plot slip velocity at location as a function of time (same location)
function slipvelPlot(SlipVel, time_, FltX, yr2sec, loc1 = 8e3)
    
    FltID = findall(abs.(FltX) .<= loc1)[1]

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig.add_subplot(111)

    ax.plot(time_/yr2sec, SlipVel[FltID, :], lw = 1)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Slip rate (m/s)")
    ax.set_title(string("Slip rate at ", loc1/1e3, "km depth"))
    ax.set_yscale("log")
    show()

    #  figname = string(dir, "/plots", name, "/sliprate.pdf")
    #  fig[:savefig](figname, dpi = 300)

end

function Vfmaxsingle(Vfmax1, t1, yr2sec)

    #  Vfmax = maximum(SlipVel, dims = 1)[:]
    
    fig = PyPlot.figure(figsize=(6,4))
    ax = fig.add_subplot(111)

    ax.plot(t1./yr2sec, Vfmax1, lw = 1)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    ax.legend(loc="upper right")
    show()

    figname = string(path, "Vfmax.png")
    fig.savefig(figname, dpi = 300)
end
# Plot Vfmax
function VfmaxPlot(Vfmax1, Vfmax2, Vfmax3, t1, t2, t3, yr2sec)

    #  Vfmax = maximum(SlipVel, dims = 1)[:]
    
    fig = PyPlot.figure(figsize=(6,4))
    ax = fig.add_subplot(111)

    ax.plot(t1./yr2sec, Vfmax1, lw = 1, label="Model 1")
    ax.plot(t2./yr2sec, Vfmax2, lw = 1, label="Model 2")
    ax.plot(t3./yr2sec, Vfmax3, lw = 1, label="Model 3")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_title("Max. slip rate on fault")
    ax.set_yscale("log")
    ax.legend(loc="upper right")
    show()

    figname = string(path, "Vfmaxcomp.png")
    fig.savefig(figname, dpi = 300)
end

function Vfmaxcompare(Vfmax1, Vfmax2, Vfmax3, t1, t2, t3, yr2sec)

    #  Vfmax = maximum(SlipVel, dims = 1)[:]
    
    fig = PyPlot.figure(figsize=(4,8))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    ax1.plot(t1, Vfmax1, lw = 1, label="Model 1")
    ax2.plot(t2, Vfmax2, lw = 1, label="Model 2")
    ax3.plot(t3, Vfmax3, lw = 1, label="Model 3")
    ax3.set_xlabel("Time (seconds from the start of the event)")
    #  ax1.set_ylabel("Max. Slip rate (m/s)")
    ax2.set_ylabel("Max. Slip rate (m/s)")
    ax1.set_xlim([-20,100])
    ax2.set_xlim([-20,100])
    ax3.set_xlim([-20,100])
    ax1.set_yscale("log")
    ax2.set_yscale("log")
    ax3.set_yscale("log")
    
    plt.suptitle("Max. slip rate on fault")
    show()

    figname = string(path, "Vfcomp.png")
    fig.savefig(figname, dpi = 300)
end


# Compare cumulative slip
function cumpare(dfsec1, dfsec2, df5yr1, df5yr2, FltX1, FltX2)
    indx = findall(abs.(FltX) .<= 18e3)[1]

    delfsec1 = dfsec1[indx:end, :]
    delfsec2 = dfsec2[indx:end, :]

    fig = PyPlot.figure(figsize=(8,7))
    ax = fig.add_subplot(111)
    

    ax.plot(df5yr2, -FltX2/1e3, color="royalblue", lw=1, alpha=1.0)
    ax.plot(delfsec2, -FltX2[indx:end]/1e3, color="chocolate", lw=1, alpha=1.0)
    ax.plot(df5yr1, -FltX1/1e3, color="royalblue", lw=1, alpha=0.5)
    ax.plot(delfsec1, -FltX1[indx:end]/1e3, color="chocolate", lw=1, alpha=0.5)
    ax.set_xlabel("Accumulated Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_title("Cumulative Slip History")
    ax.set_ylim([0,24])
    #  ax.set_xlim([5,15])
    ax.invert_yaxis()
    show()
    
    figname = string(path, "cumpare.png")
    fig.savefig(figname, dpi = 300)

end


