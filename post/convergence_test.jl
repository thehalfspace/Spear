#   Convergence test comparison plots

# Plotnumber of earthquakes in magnitude bins
function resolution_test2(M100, M160, M200, M250, M320, M400)
    
    ref23 = length(M100[2.5 .<= M100 .< 3.5])
    n_160_23 = length(M160[2.5 .<= Mw160 .< 3.5])
    n_200_23 = length(M200[2.5 .<= Mw200 .< 3.5])
    n_250_23 = length(M250[2.5 .<= Mw250 .< 3.5])
    n_320_23 = length(M320[2.5 .<= Mw320 .< 3.5])
    n_400_23 = length(M400[2.5 .<= Mw400 .< 3.5])

    y_axis_23 = [n_160_23-ref23, n_200_23-ref23, n_250_23-ref23, n_320_23-ref23, n_400_23-ref23]
    x_axis = [160, 200, 250, 320, 400]
    
    ref34 = length(M100[3.5 .<= M100 .< 4.5])
    n_160_34 = length(M160[3.5 .<= Mw160 .< 4.5])
    n_200_34 = length(M200[3.5 .<= Mw200 .< 4.5])
    n_250_34 = length(M250[3.5 .<= Mw250 .< 4.5])
    n_320_34 = length(M320[3.5 .<= Mw320 .< 4.5])
    n_400_34 = length(M400[3.5 .<= Mw400 .< 4.5])

    y_axis_34 = [n_160_34-ref34, n_200_34-ref34, n_250_34-ref34, n_320_34-ref34, n_400_34-ref34]
    
    ref45 = length(M100[4.5 .<= M100 .< 5.5])
    n_160_45 = length(M160[4.5 .<= Mw160 .< 5.5])
    n_200_45 = length(M200[4.5 .<= Mw200 .< 5.5])
    n_250_45 = length(M250[4.5 .<= Mw250 .< 5.5])
    n_320_45 = length(M320[4.5 .<= Mw320 .< 5.5])
    n_400_45 = length(M400[4.5 .<= Mw400 .< 5.5])

    y_axis_45 = [n_160_45-ref45, n_200_45-ref45, n_250_45-ref45, n_320_45-ref45, n_400_45-ref45]

    ref5 = length(M100[5.5 .<= M100 .< 9])
    n_160_5 = length(M160[5.5 .<= Mw160 .< 9])
    n_200_5 = length(M200[5.5 .<= Mw200 .< 9])
    n_250_5 = length(M250[5.5 .<= Mw250 .< 9])
    n_320_5 = length(M320[5.5 .<= Mw320 .< 9])
    n_400_5 = length(M400[5.5 .<= Mw400 .< 9])

    y_axis_5 = [n_160_5-ref5, n_200_5-ref5, n_250_5-ref5, n_320_5-ref5, n_400_5-ref5]
    
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    ax.scatter(x_axis, y_axis_23, label="Mw bin 2.5-3.5")
    ax.scatter(x_axis, y_axis_34, label="Mw bin 3.5-4.5")
    ax.scatter(x_axis, y_axis_45, label="Mw bin 4.5-5.5")
    ax.scatter(x_axis, y_axis_5, label="Mw bin >= 5.5")
    ax.set_xlabel("Element size (m)")
    ax.set_ylabel("rms difference")
    ax.set_title("No of earthquakes in different magnitude bins")
    plt.legend()
    show()

    figname = string(path, "resolution_noe.png")
    fig.savefig(figname, dpi = 300)

end


# average coseismic slip rms diff
function cosavg(res, rms)
    fig = PyPlot.figure(figsize=(8,7))
    ax = fig.add_subplot(111)

    ax.plot(res, rms, "k.--", markersize=20)
    
    #  ax.set_xlabel("Average node spacing (m)")
    #  ax.set_ylabel("rms difference")
    #  ax.set_title("rms difference of average coseismic slip (reference = 25 m)")
    #  ax.set_ylim([0,5e-4])
    #  plt.legend()
    show()

    figname = string(path, "diff_slip_rms2.png")
    fig.savefig(figname, dpi = 300)
end

# RMS difference of differential slip
function rms_diff_slip(rms1, res1, rms2, res2, rms3, res3, rms4, res4, rms5, res5, ref)
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    ax.scatter(res1, abs(ref-rms1))#, label="400m")
    ax.scatter(res2, abs(ref-rms2))#, label="320m")
    ax.scatter(res3, abs(ref-rms3))#, label="250m")
    ax.scatter(res4, abs(ref-rms4))#, label="200m")
    ax.scatter(res5, abs(ref-rms5))#, label="160m")
    #  ax.scatter(res6, rms6, label="Element Size = 100m")
    ax.set_xlabel("Average node spacing (m)")
    ax.set_ylabel("rms difference")
    ax.set_title("rms difference of average coseismic slip with respect to 25 m spacing")
    #  ax.set_ylim([0,5e-4])
    #  plt.legend()
    show()

    figname = string(path, "diff_slip_rms2.png")
    fig.savefig(figname, dpi = 300)


end

# Timesteps
function timesteps(t100, t160, t200, t250, t320, t400)

    y_axis = [t400-t100, t320-t100, t250-t100, t200-t100, t160-t100]
    x_axis = [160, 200, 250, 320, 400]
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)
    
    ax.scatter(x_axis, abs.(y_axis))
    ax.set_xlabel("Element size (m)")
    ax.set_ylabel("rms difference")
    ax.set_title("Number of timesteps")
    show()

    figname = string(path, "timesteps1.png")
    fig.savefig(figname, dpi = 300)
end

# Max velocity
function max_vel(v100, v160, v200, v250, v320, v400)

    y_axis = [v400-v100, v320-v100, v250-v100, v200-v100, v160-v100]
    x_axis = [160, 200, 250, 320, 400]
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)
    
    ax.scatter(x_axis, abs.(y_axis))
    ax.set_xlabel("Element size (m)")
    ax.set_ylabel("rms difference")
    ax.set_title("Maximum velocity on fault")
    show()

    figname = string(path, "max_vel.png")
    fig.savefig(figname, dpi = 300)
end

# Vfmax convergence
function Vfcon(Vf100, t100, Vf160, t160,  Vf200, t200, Vf250, t250, Vf320, t320, Vf400, t400, yr2sec)
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    ax.plot(t100./yr2sec, Vf100, lw = 1, label="Element Size = 100m")
    ax.plot(t160./yr2sec, Vf160, lw = 1, label="Element Size = 160m")
    ax.plot(t200./yr2sec, Vf200, lw = 1, label="Element Size = 200m")
    ax.plot(t250./yr2sec, Vf250, lw = 1, label="Element Size = 250m")
    ax.plot(t320./yr2sec, Vf320, lw = 1, label="Element Size = 3200m")
    ax.plot(t400./yr2sec, Vf400, lw = 1, label="Element Size = 400m")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Vfmax")
    ax.set_title("Maximum sliprate on fault")
    plt.legend()
    show()

    figname = string(path, "Vfmax_convg.png")
    fig.savefig(figname, dpi = 300)

end

# 1. Compare the differential slip for the largest earthquake
function diff_slip(delf1, fltx1, delf2, fltx2, delf3, fltx3, delf4, fltx4)
    
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    ax.plot(delf1, -fltx1./1e3, lw = 1, label="Element Size = 250m")
    ax.plot(delf2, -fltx2./1e3, lw = 1, label="Element Size = 200m")
    ax.plot(delf3, -fltx3./1e3, lw = 1, label="Element Size = 160m")
    ax.plot(delf4, -fltx4./1e3, lw = 1, label="Element Size = 100m")
    ax.set_xlabel("Differential Slip")
    ax.set_ylabel("Depth (km)")
    ax.set_title("Coseismic slip (Long Fault Zone)")
    ax.set_ylim([0,24])
    ax.invert_yaxis()
    plt.legend()
    show()

    figname = string(path, "diff_slip_large.png")
    fig.savefig(figname, dpi = 300)
end

# Plot shear stress for events: 1 sec time interval
function res_plot(Mw100, Mw160, Mw200, Mw250, Mw320, Mw400)
    fig = PyPlot.figure(figsize=(12,8))
    ax = fig.add_subplot(111)

    ax.plot(Mw100, "o-", label="100 m")
    ax.plot(Mw160, "o-", label="160 m")
    ax.plot(Mw200, "o-", label="200 m")
    ax.plot(Mw250, "o-", label="250 m")
    ax.plot(Mw320, "o-", label="320 m")
    ax.plot(Mw400, "o-", label="400 m")

    ax.set_xlabel("Event number")
    ax.set_ylabel("Magnitude")
    ax.set_title("Resolution test")
    ax.legend()
    show()

    #  figname = string(dir, "/plots", name, "/fric.png")
    figname = string(path, "pulse5.png")
    fig.savefig(figname, dpi = 300)

end

# Plot shear stress for events: 1 sec time interval
function res_plot2(tausec, FltX)
    fig = PyPlot.figure(figsize=(12,8))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    ax1.plot(tausec[:,1], -FltX./1e3)
    ax2.plot(tausec[:,2], -FltX./1e3)
    ax3.plot(tausec[:,3], -FltX./1e3)

    ax2.set_xlabel("Sliprate (m/s)")
    ax1.set_ylabel("Depth (km)")
    plt.suptitle("Seismic sliprate evolution")
    ax1.set_ylim([0,24])
    ax2.set_ylim([0,24])
    ax3.set_ylim([0,24])
    ax1.set_xlim([0.001,1])
    ax2.set_xlim([0.001,1])
    ax3.set_xlim([0.001,1])
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    ax3.invert_yaxis()
    show()

    #  figname = string(dir, "/plots", name, "/fric.png")
    figname = string(path, "pulse6.png")
    fig.savefig(figname, dpi = 300)

end

# Plot spatiotemporal sliprate
function spatiotemporal(tausec, FltX, rup_time)

    indx = findall(abs.(FltX) .<= 18e3)[1]
    tausec = tausec[indx:end, :]
    
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    #  ims = ax.imshow(tausec, vmin=1e-4, vmax=3, extent=[0,rup_time, 0, 17], aspect="auto", interpolation="bilinear")
    ims = ax.imshow(tausec, vmin=25, vmax=40, extent=[0,rup_time, 0, 17], aspect="auto", interpolation="bilinear")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Depth (km)")
    #  ax.set_title("Seismic shear stress evolution (MPa)")
    #  ax.set_title("Seismic sliprate evolution (m/s)")
    #  ax[.set_ylim]([0,15])
    ax.invert_yaxis()
    plt.colorbar(ims)
    show()

    #  figname = string(dir, "/plots", name, "/fric.png")
    figname = string(path, "stressevol.png")
    fig.savefig(figname, dpi = 300)

end


# Plot timesteps
function timesteps(time250, time200, time160, time100)
    
    fig = PyPlot.figure()
    ax = figadd_subplot(111)

    ax.plot(time250, label="250 m")
    ax.plot(time200, label="200 m")
    ax.plot(time160, label="160 m")
    ax.plot(time100, label="100 m")
    ax.set_xlabel("iterations")
    ax.set_ylabel("Vfmax (m/s)")
    ax.set_title("Convergence test")
    plt.legend()
    show()

    figname = string(path, "Vfmax_convg.png")
    fig.savefig(figname, dpi = 300)
end

# rupture time difference
function rup_time(time250, time200, time160, time100)
    
    fig = PyPlot.figure()
    ax = fig.add_subplot(111)

    ax.plot(time250, label="250 m")
    ax.plot(time200, label="200 m")
    ax.plot(time160, label="160 m")
    ax.plot(time100, label="100 m")
    ax.set_xlabel("iterations")
    ax.set_ylabel("Time (s)")
    ax.set_title("Convergence test")
    plt.legend()
    show()

    figname = string(path, "timestepsC_convg.png")
    fig.savefig(figname, dpi = 300)
end
