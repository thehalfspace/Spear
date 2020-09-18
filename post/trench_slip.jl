#############################
# PLOT SLIP AT VARIOUS DEPTHS 
#############################

function plot_trench_slip(delfsec, delf5yr, FltX)
    """Input parameters:
        delfsec: slip contours every 0.5 seconds during
        the seismic period

        delf5yr: slip contours every 2 years during
        the interseismic period

        FltX: Location of nodes on the fault"""


    # Depth valeus
    d0 = findall(abs.(FltX) .< 1)[1]
    d1 = findall(abs.(FltX) .< 150)[1]
    d2 = findall(abs.(FltX) .< 300)[1]
    d3 = findall(abs.(FltX) .< 500)[1]
    d4 = findall(abs.(FltX) .< 1000)[1]
    d5 = findall(abs.(FltX) .< 5000)[1]
    d6 = findall(abs.(FltX) .< 10000)[1]

    fig = PyPlot.figure(figsize=(8,6), dpi = 300)
    ax = fig.add_subplot(111)
    
    plt.rc("font",size=12)
    ax.plot(delfsec[d0,:], label=:"0 m depth", lw=2)
    ax.plot(delfsec[d1,:], label=:"150 m depth", lw=2)
    #  ax.plot(delfsec[d2,:], label=:"300 m depth", lw=2)
    #  ax.plot(delfsec[d3,:], label=:"500 m depth", lw=2)
    #  ax.plot(delfsec[d4,:], label=:"1 km depth", lw=2)
    ax.plot(delfsec[d5,:], label=:"5 km depth", lw=2)
    ax.plot(delfsec[d6,:], label=:"10 km depth", lw=2)

    ax.set_xlabel("Timestep")
    ax.set_ylabel("Seismic Slip")
    ax.get_yaxis().set_tick_params(which="both", direction="in")
    ax.get_xaxis().set_tick_params(which="both", direction="in")
    #  plt.rc("grid", linestyle="--", color="black", alpha=0.5)
    #  plt.grid("True")
    plt.legend()
    
    show()

    figname = string(path, "trench_slip.png")
    fig.savefig(figname, dpi = 300)
end
