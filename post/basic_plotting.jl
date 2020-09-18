####################################################
#       BASIC PLOTTING FOR EARTHQUAKE CYCLES
####################################################

# Includes plotting cumulative slip contours, peak slip rate on fault,
# hypocenter location, friction parameters and initial stresses,
# earthquake catalogue, slip rate and stresses at a specified location
# on the fault.

using Plots
using LaTeXStrings
using StatsBase
#  pgfplots()  # backed for latex-based plotting:
               # good for publication quality figures

#  gr()    # backend for python based pyplot: good for quick viewing

#   Plots.scalefontsizes(2)
pyplot()

#-------------------------------
# Initial Friction and Stresses
#-------------------------------
function plot_friction(cca, ccb, Seff, tauo, FltX)

    plt = plot(framestyle=[:box :grid],fontfamily=font(28, "Helvetica"), dpi=200)

    plot!(Seff/1e6, -FltX/1e3, axis=:bottom, label=:"Normal Stress", legend=:bottomright, lw = 2)
    plot!(tauo/1e6, -FltX/1e3, axis=:bottom, label=:"Shear Stress", lw = 2)
    plot!((cca-ccb)*1e3 .+ 10 , axis=:top, -FltX/1e3, label=:"Friction Parameters", lw = 2)
    #  xlabel!("Value")
    yaxis!(L"Depth (km)", yflip=true)
    ylims!(0, 24)

    #  figname = string(dir, "/plots", name, "/fric.pdf")
    #  savefig(string(path, "initial_parameters.svg"))
    savefig(string(path, "initial_parameters.pdf"))

    plt
end

#-------------------------
# CUMULATIVE SLIP CONTOURS
#-------------------------
function plot_cumulative_slip(delfsec, delf5yr, FltX)
    """Input parameters:
        delfsec: slip contours every 0.5 seconds during
        the seismic period

        delf5yr: slip contours every 2 years during
        the interseismic period

        FltX: Location of nodes on the fault"""

    indx = findall(abs.(FltX) .<= 18e3)[1]
    delfsec2 = delfsec[indx:end, :]

    plt = plot(framestyle=[:box],grid=false, size=(600,300), dpi=300)

    plot!(delf5yr, -FltX/1e3, lc=:steelblue,label=:"", lw=0.5)
    plot!(delfsec2, -FltX[indx:end]/1e3, lc=:peru, label=:"", lw=0.5)

    xaxis!(L"Accumulated\ Slip\ (m)"); #xlims!(10,30); #xticks!(-1:0.2:1)
    yaxis!(L"Depth\ (km)", yflip=true); ylims!(0,20); #yticks!(0:0.1:1)
    #  savefig(string(path, "cumulative_slip.svg"))
    savefig(string(path, "cumulative_slip.png"))

    plt
end

function plot_Vfmax(Vfmax, time_)
    """Input parameters
        Vfmax: peak slip rate on-fault
        time_: time in years"""

    plt = plot(framestyle=[:box], grid=false,dpi=200)

    plot!(time_, log10.(Vfmax), lc=:steelblue, label=:"",lw=1.5)

    xaxis!(L"Time\ (yr)") 
    yaxis!(L"Peak\ slip\ rate\ on\ fault\ (m/s)")
    #  savefig(string(path, "Vfmax.svg"))
    savefig(string(path, "Vfmax_02.png"))

    plt
end

function plot_Mw(Mw)

    hist = fit(Histogram, Mw, nbins = 9)

    # Cumulative
    cum = cumsum(hist.weights[end:-1:1])[end:-1:1]
    plt = plot(framestyle=[:box :grid],fontfamily=font(28, "Helvetica"), dpi=200)

    scatter!(hist.edges[1][1:end-1], cum, ms=5, label=:"") #, label="Cumulative")

    xaxis!(L"Moment\ Magnitude (Mw)"); #xgrid!(:on, :minorgrid)
    yaxis!(L"Number\ of\ Earthquakes", :log); #ylims!(10^1.45, 10^1.65)
    #gui()
    #  savefig(string(path, "MwPlot.png"))
    savefig(string(path, "MwPlot.pdf"))

    plt
end

function plot_alphaa(alphaa, time_)
    plt = plot(framestyle=[:box], grid=false, size=(600,300),dpi=300)

    plot!(time_[alphaa .< 1.0], alphaa[alphaa .<1.0], lc=:darkblue, label=:"",lw=1.5)

    xaxis!(L"Time\ (yr)") #xminorgrid=true
    #  xticks!(0:100:maximum(time_));
    yaxis!(L"\alpha_D")
    #  savefig(string(path, "Vfmax.svg"))
    savefig(string(path, "alpha_02.png"))

    plt

end

function plot_stress_evolution(seismic_stress, FltX, rupture_time)

    indx = findall(abs.(FltX) .<= 18e3)[1]
    seismic_stress = seismic_stress[indx:end, :]

    plt = plot(framestyle=[:box], size=(600,400), dpi=200)

    heatmap!(seismic_stress)
    savefig(string(path, "heatmap.pdf"))

    plt
end

function plot_stress(seismic_stress, FltX)

    l = length(seismic_stress[1,:])

    plt = plot(framestyle=[:box], grid=false, size=(600,400), dpi=200)

    for i in 1:l
        plot!(seismic_stress[:,i].*1e-0 .+ 2*i, -FltX/1e3, label=:"", lc=:black, lw=0.5)
    end

    xaxis!(L"Shear\ stress\ every\ second"); #xticks!(0:20:l)
    yaxis!(L"Depth\ (km)", yflip=true); ylims!(0,15)

    savefig(string(path, "ssevol3.pdf"))

    plt

end


# Plot stress heterogeneity: norm of differential shear stress
function plot_stress_heterogeneity(stress_dam, stress_homo)

    ld = length(stress_dam[1,:])
    lh = length(stress_homo[1,:])

    nd = zeros(ld)
    nh = zeros(lh)

    plot(framestyle=[:box :grid], dpi=300)

    for i in 1:ld
        nd[i] = norm(stress_dam[:,i])
    end

    for i in 1:lh
        nh[i] = norm(stress_homo[:,i])
    end

    scatter!(range(1,stop=ld), nd, msa=0, color=:darkblue, label=:"")
    plot!(range(1,stop=ld), nd,lc=:darkblue, label=:"Damaged zone")

    scatter!(range(1,stop=lh), nh, msa=0, color=:darkred, label=:"")
    plot!(range(1,stop=lh), nh,lc=:darkred, label=:"Homogeneous medium")

    xaxis!("Timesteps (s)"); #xlims!(150,350)
    yaxis!("Norm of differential shear stress (MPa)")

    gui()
    savefig(string(path, "sshetero.png"))

end
