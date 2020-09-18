##############################
#  PLOTTING SCRIPTS
##############################

using PyPlot
using StatsBase
using LaTeXStrings
using PyCall
mpl = pyimport("matplotlib")

# Default plot params
function plot_params()
  plt.rc("xtick", labelsize=16)
  plt.rc("ytick", labelsize=16)
  plt.rc("xtick", direction="in")
  plt.rc("ytick", direction="in")
  plt.rc("font", size=15)
  plt.rc("figure", autolayout="True")
  plt.rc("axes", titlesize=16)
  plt.rc("axes", labelsize=17)
  plt.rc("xtick.major", width=1.5)
  plt.rc("xtick.major", size=5)
  plt.rc("ytick.major", width=1.5)
  plt.rc("ytick.major", size=5)
  plt.rc("lines", linewidth=2)
  plt.rc("axes", linewidth=1.5)
  plt.rc("legend", fontsize=13)
  plt.rc("mathtext", fontset="stix")
  plt.rc("font", family="STIXGeneral")

  # Default width for Nature is 7.2 inches, 
  # height can be anything
  #plt.rc("figure", figsize=(7.2, 4.5))
end

# Plot slip vs event number
function slipPlot(delfafter2, rupture_len, FltX, Mw, tStart)
    plot_params()
    fig, ax = PyPlot.subplots(nrows=1, ncols=4, sharex="all", sharey="all", figsize=(9.2, 5.00))

    xaxis = tStart[Mw .>2.8]
    delfafter = delfafter2[:,Mw .> 2.8]
    Mw2 = Mw[Mw .> 2.8]

    # Normalize colorbar
    norm = matplotlib.colors.Normalize(vmin = minimum(Mw2), 
                                       vmax=maximum(Mw2)) 
    colors = matplotlib.cm.inferno_r(norm(Mw2))

    ax[1].barh(xaxis, delfafter[end-1,:], height=6, 
              color=colors, align="center"); 
    ax[1].set_ylabel("Time (yr)")
    ax[1].invert_yaxis()
    ax[1].set_title("At 60 m depth")

    trench_depth1 = findall(abs.(FltX) .< 4.0e3)[1]
    trench_depth2 = findall(abs.(FltX) .< 6.0e3)[1]
    trench_depth3 = findall(abs.(FltX) .< 8.0e3)[1]
    
    ax[2].barh(xaxis, delfafter[trench_depth1,:], height=6, 
              color=colors, align="center"); 
    ax[2].invert_yaxis()
    ax[2].set_title("At 4 km depth")
    
    ax[3].barh(xaxis, delfafter[trench_depth2,:], height=6, 
              color=colors, align="center"); 
    ax[3].invert_yaxis()
    ax[3].set_title("At 6 km depth")
    
    ax[4].barh(xaxis, delfafter[trench_depth3,:], height=6, 
              color=colors, align="center"); 
    ax[4].invert_yaxis()
    ax[4].set_title("At 8 km depth")
    
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap="inferno_r")
    sm.set_array([])
    fig.colorbar(sm, shrink=0.9, label="Mw")
    plt.xlabel("Coseismic Slip (m)")
    plt.tight_layout()
    show()
    
    figname = string(path, "coseismic_slip.png")
    fig.savefig(figname, dpi = 300)
end

# Cumulative sliprate plot
function eqCyclePlot(sliprate, FltX)
    indx = findall(abs.(FltX) .<= 16e3)[1]
    value = sliprate[indx:end,10000:end]
    
    depth = -FltX[indx:end]./1e3

    plot_params()
    fig = PyPlot.figure(figsize=(9.2, 4.45))
    ax = fig.add_subplot(111)
    
    c = ax.imshow(value, cmap="inferno", aspect="auto",
                  norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e0),
                  interpolation="bicubic",
                  extent=[0,length(sliprate[1,:]), 0,16])
    
    ax.set_xlabel("Timesteps")
    ax.set_ylabel("Depth (km)")

    ax.invert_yaxis()
    cbar = fig.colorbar(c)
    #   cbar.set_ticks(cbar.get_ticks()[1:2:end])
    
    show()
    figname = string(path, "slr03.png")
    fig.savefig(figname, dpi = 300)
    
end

# spatiotemporal imshow
function sptempPlot(seismic_slipvel, FltX)
    
    indx = findall(abs.(FltX) .<= 16)[1]
    value = transpose(seismic_slipvel[:, indx:end])
    
    depth = -FltX[indx:end]

    plot_params()
    fig = PyPlot.figure(figsize=(12.2, 4.45))
    ax = fig.add_subplot(111)

    # for sliprate
    c = ax.imshow(value, cmap="inferno", aspect="auto",
                  norm=matplotlib.colors.LogNorm(vmin=1e-12, vmax=1e0),
                  interpolation="bicubic",
                  extent=[0,length(seismic_slipvel), 0,16])

    # for stress
    #  c = ax.imshow(value, cmap="inferno", aspect="auto",
                  #  vmin=22.5, vmax=40,
                  #  interpolation="bicubic",
                  #  extent=[0,length(seismic_slipvel[1,:]), 0,16])
    
    ax.set_xlabel("Timesteps")
    ax.set_ylabel("Depth (km)")

    ax.invert_yaxis()
    cbar = fig.colorbar(c)
    #   cbar.set_ticks(cbar.get_ticks()[1:2:end])
    
    show()
    figname = string(path, "slr03.png")
    fig.savefig(figname, dpi = 300)
end

# Plot Vfmax
function VfmaxPlot(Vfmax, time_, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)
    
    ax.plot(time_./yr2sec, Vfmax, lw = 2.0)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    #  ax.set_xlim([230,400])
    show()
    
    figname = string(path, "Vfmax01.png")
    fig.savefig(figname, dpi = 300)
end

# Plot alpha
function alphaaPlot(alphaa, time_, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)

    ax.plot(time_./yr2sec, alphaa, lw = 2)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Shear Modulus Contrast (%)")
    #  ax.set_xlim([230,400])
    show()


    figname = string(path, "alphaa_01.png")
    fig.savefig(figname, dpi = 300)
end

# Plot cumulative slip
function cumSlipPlot(delfsec, delfyr, FltX)
    indx = findall(abs.(FltX) .<= 18)[1]

    delfsec2 = transpose(delfsec[:,indx:end])
    delfyr2 = transpose(delfyr)

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    plt.rc("font",size=12)

    ax.plot(delfyr2, FltX, color="royalblue", lw=1.0)
    ax.plot(delfsec2, FltX[indx:end], color="chocolate", lw=1.0)
    ax.set_xlabel("Accumulated Slip (m)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,24])
    #  ax.set_xlim([1,20])
    
    ax.invert_yaxis()
    
    show()
    
    figname = string(path, "cumslip_02.png")
    fig.savefig(figname, dpi = 300)

end

# Plot friction parameters
function icsPlot(a_b, Seff, tauo, FltX)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.plot(Seff/1e6, -FltX./1e3, "k-", label="Normal Stress")
    ax.plot(tauo/1e6, -FltX./1e3, "k--", label="Shear Stress")
    ax.set_xlabel("Stresses (MPa)")
    ax.set_ylabel("Depth (km)")
    plt.legend(loc="lower right") 
    
    col="tab:blue"
    ax2 = ax.twiny()
    ax2.plot(a_b, -FltX./1e3, label="(a-b)")
    ax2.set_xlabel("Rate-state friction value", color=col)
    ax2.get_xaxis().set_tick_params(color=col)
    ax2.tick_params(axis="x", labelcolor=col)
    
    ax.set_ylim([0,16])
    ax.invert_yaxis()
    show()
    
    figname = string(path, "ics_02.png")
    fig.savefig(figname, dpi = 300)
end
