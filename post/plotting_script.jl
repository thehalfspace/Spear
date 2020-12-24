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


# Plot alpha and Vfmax on the same plot
function healing_analysis(Vf, alphaa, t, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.plot(t./yr2sec, Vf, lw = 2.0, label="Max. Slip rate")
    lab1 = "Max. slip rate"
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    ax.set_xlim([0, 250])
    
    col="tab:red"
    ax2 = ax.twinx()
    
    ax2.plot(t./yr2sec, alphaa.*100, c=col, lw=2.0, label="Shear modulus ratio")
    lab2 = "Shear modulus ratio"
    ax.set_xlabel("Time (years)")
    ax2.set_ylabel("Shear Modulus (% of host rock)")
    ax2.set_ylim([35, 60])
    #  ax2.set_ylim([75, 100])
    ax2.get_xaxis().set_tick_params(color=col)
    ax2.tick_params(axis="x", labelcolor=col)

    #  ax.legend([lab1, lab2], loc=0)
    show()
    
    figname = string(path, "healing_analysis.png")
    fig.savefig(figname, dpi = 300)
end

# Plot the stressdrops after each earthquake
function stressdrop_2(taubefore, tauafter, FltX)
    plot_params()
  
    i = 1;
    #  for i in 1:length(stressdrops[1,:])
      fig = PyPlot.figure(figsize=(7.2, 4.45));
      ax = fig.add_subplot(111);
      ax.plot(taubefore, FltX, lw = 2.0, color="tab:orange", 
              label="Shear stress before the earthquake", alpha=1.0);
      ax.plot(tauafter, FltX, lw = 2.0, color="tab:blue", 
              label="Shear stress after the earthquake", alpha=1.0);
      ax.set_xlabel("Stress drop (MPa)");
      ax.set_ylabel("Depth (km)");
      ax.set_ylim([0,24]);
      ax.set_xlim([15,45]);
      ax.invert_yaxis();
      plt.legend();
      show()
      
      figname = string(path, "shear_stress_im_",i,".png");
      fig.savefig(figname, dpi = 300);
    #  end
end


# Plot shear stress comparison
function shear_stress_comp(shear1b, shear1a, shear2b, shear2a, FltX1, FltX2)
    plot_params()
   
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    ax.plot(shear1b, FltX1, lw = 2.0, color="tab:blue",ls=:dashed, 
            label="Immature Fault Zone: before")
    ax.plot(shear1a, FltX1, lw = 2.0, color="tab:blue", label="Immature Fault Zone: after")
    ax.plot(shear2b, FltX2, lw = 2.0, color="tab:orange", ls=:dashed, 
            label="Mature Fault Zone: before")
    ax.plot(shear2a, FltX2, lw = 2.0, color="tab:orange", label="Mature Fault Zone: after")
    ax.set_xlabel("Shear stress (MPa)")
    ax.set_ylabel("Depth (km)")
    ax.set_ylim([0,24])
    ax.invert_yaxis()
    plt.legend()
    show()
    
    figname = string(path, "Shear_Stress_002.png");
    fig.savefig(figname, dpi = 300);

end

# Plot rupture_length vs event time
function stem_plot(rl1, rl2, rl3, rl4)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.stem(-rl1./1e3, linefmt="C0-", markerfmt="C0o", basefmt=:None, label="Healing time = 10 yr")
    ax.stem(-rl2./1e3, linefmt="C1-", markerfmt="C1o", basefmt=:None, label="Healing time = 12 yr")
    ax.stem(-rl3./1e3, linefmt="C2-", markerfmt="C2o", basefmt=:None, label="Healing time = 15 yr")
    #  ax.stem(-rl4./1e3, linefmt="C3-", markerfmt="C3o", basefmt=:None, label="Healing time = 20 years")
    ax.set_xlabel("Event number")
    #  ax.set_xlabel("Event time (years)")
    ax.set_ylabel("Rupture length (km)")
    #  ax.set_xlim([10,300])
    plt.legend()
    show()
    
    figname = string(path, "stem_plot2.png")
    fig.savefig(figname, dpi = 300)

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
    indx = findall(abs.(FltX) .<= 16)[1]
    value = sliprate[indx:end,:]
    
    depth = FltX[indx:end]

    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    c = ax.imshow(value, cmap="inferno", aspect="auto",
                  norm=matplotlib.colors.LogNorm(vmin=1e-9, vmax=1e0),
                  interpolation="bicubic",
                  extent=[0,length(sliprate[1,:]), 0,16])
    
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
    figname = string(path, "mature_sliprate_3.png")
    fig.savefig(figname, dpi = 300)
    
end

# Plot Vfmax
function VfmaxPlot(Vfmax, t, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)
    
    ax.plot(t./yr2sec, Vfmax, lw = 2.0)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Max. Slip rate (m/s)")
    ax.set_yscale("log")
    #  ax.set_xlim([230,400])
    show()
    
    figname = string(path, "Vfmax01.png")
    fig.savefig(figname, dpi = 300)
end

# Plot Vsurface
function VsurfPlot(Vsurf10, Vsurf12, Vsurf15, t10, t12, t15, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(8.2, 6.00))
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    
    ax.plot(t10./yr2sec, Vsurf10, lw = 2.0, label="Healing time = 10 yr")
    ax2.plot(t12./yr2sec, Vsurf12, lw = 2.0, label="Healing time = 12 yr")
    ax3.plot(t15./yr2sec, Vsurf15, lw = 2.0, label="Healing time = 15 yr")
    ax3.set_xlabel("Time (years)")
    ax2.set_ylabel("Surface. Slip rate (m/s)")
    ax.set_yscale("log")
    ax2.set_yscale("log")
    ax3.set_yscale("log")
    #  ax.set_xlim([230,400])
    #  plt.legend()
    plt.tight_layout()
    show()
    
    figname = string(path, "Vsurf02.png")
    fig.savefig(figname, dpi = 300)
end

# Plot alpha
function alphaaPlot(alphaa, t, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)

    ax.plot(t./yr2sec, alphaa.*100, lw = 2)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Shear Modulus Contrast (%)")
    #  ax.set_xlim([230,400])
    show()


    figname = string(path, "alpha_01.png")
    fig.savefig(figname, dpi = 300)
end

# Compare alpha
function alphaComp(a1, t1, a2, t2, a3, t3, yr2sec)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 3.45))
    ax = fig.add_subplot(111)

    ax.plot(t1./yr2sec, a1.*100, lw = 2, label="10 yr")
    ax.plot(t2./yr2sec, a2.*100, lw = 2, label="12 yr")
    ax.plot(t3./yr2sec, a3.*100, lw = 2, label="15 yr")
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Shear Modulus Ratio (%)")
    #  ax.set_xlim([230,400])
    legend()
    show()


    figname = string(path, "alpha_comp.png")
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
    ax.set_xlim([0,maximum(delfsec2)])
    ax.set_xlim([0,9.0])
    
    ax.invert_yaxis()
    
    show()
    
    figname = string(path, "cumulative_slip.png")
    fig.savefig(figname, dpi = 300)

end

# Plot friction parameters
function icsPlot(a_b, Seff, tauo, FltX)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.plot(Seff, FltX, "k-", label="Normal Stress")
    ax.plot(tauo, FltX, "k--", label="Shear Stress")
    ax.set_xlabel("Stresses (MPa)")
    ax.set_ylabel("Depth (km)")
    plt.legend(loc="lower right") 
    
    col="tab:blue"
    ax2 = ax.twiny()
    ax2.plot(a_b, FltX, label="(a-b)")
    ax2.set_xlabel("Rate-state friction value", color=col)
    ax2.get_xaxis().set_tick_params(color=col)
    ax2.tick_params(axis="x", labelcolor=col)
    
    ax.set_ylim([0,16])
    ax.invert_yaxis()
    show()
    
    figname = string(path, "ics_02.png")
    fig.savefig(figname, dpi = 300)
end


# Plot stressdrop comparison
function sd_comp(ds_im, tS_im, ds_m, tS_m)
    plot_params()
    fig = PyPlot.figure(figsize=(7.2, 4.45))
    ax = fig.add_subplot(111)
    
    ax.scatter(tS_im, ds_im, color="tab:blue", label="Immature")
    ax.scatter(tS_m, ds_m, color="tab:orange", label="Mature")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Stress drops (MPa)")
    plt.legend() 
    show()
    
    figname = string(path, "del_sigma_tStart.png")
    fig.savefig(figname, dpi = 300)
end
