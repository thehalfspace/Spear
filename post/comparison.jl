#################################################################
#
#   SCRIPT FOR COMPARISON OF RESULTS (MFD, CUMULATIVE SLIP, ETC)
#
#################################################################


#...........
# Plot MFD
#...........
function MwPlot(Mw500, Mw1500, Mw3000, Mw5000, Mw7000)

    hist500 = fit(Histogram, Mw500, nbins = 10)
    hist1500 = fit(Histogram, Mw1500, nbins = 10)
    hist3000 = fit(Histogram, Mw3000, nbins = 10)
    hist5000 = fit(Histogram, Mw5000, nbins = 10)
    hist7000 = fit(Histogram, Mw7000, nbins = 10)

    # Cumulative
    cum500 = cumsum(hist500.weights[end:-1:1])[end:-1:1]
    cum1500 = cumsum(hist1500.weights[end:-1:1])[end:-1:1]
    cum3000 = cumsum(hist3000.weights[end:-1:1])[end:-1:1]
    cum5000 = cumsum(hist5000.weights[end:-1:1])[end:-1:1]
    cum7000 = cumsum(hist7000.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig.add_subplot(111)

    ax.plot(hist500.edges[1][1:end-1], cum500, ".", label="0.5 km width")
    ax.plot(hist1500.edges[1][1:end-1], cum1500, "*", label="1.5 km width")
    ax.plot(hist3000.edges[1][1:end-1], cum3000, ":", label="3.0 km width")
    ax.plot(hist5000.edges[1][1:end-1], cum5000, "^", label="5.0 km width")
    ax.plot(hist7000.edges[1][1:end-1], cum7000, "+", label="7.0 km width")
    ax.set_xlabel("Moment Magnitude (Mw)")
    ax.set_ylabel("Number of Earthquakes")
    ax.set_yscale("log")
    ax.set_title("Magnitude-frequency distribution")
    ax.legend(loc="upper right")
    show()

    figname = string(path, "mfd_comp.png")
    fig.savefig(figname, dpi = 300)
end
