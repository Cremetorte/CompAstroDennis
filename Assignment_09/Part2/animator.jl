using Plots
try
    using DataFrames
    using CSV
catch
    println("Could not import all necessary Packages. Do you wish to install DataFrames and CSV right now?")
    println("type y or n: ")
    usrinput = readline()
    if usrinput == "y"
        using Pkg
        Pkg.add("DataFrames")
        Pkg.add("CSV")
    else
        error("Packages could not be installed")
    end
end
using DataFrames
using CSV

gr() # Use the GR backend for plotting


x(index) = (index-1) * 0.2


filename = "output.csv"


df = CSV.read(filename, DataFrame)
data = Matrix(df)

t = data[:, 1]
ρ = data[:, 2:end]

# println("Shape of data: ", size(data))

nr_frames = size(data)[1]
gif_length = 10 # seconds


# Create a GIF
anim = @animate for i in 1:10:nr_frames
    print("Calculating frame $(i) / $(nr_frames)...   \r")
    plot(x.(eachindex(ρ[i, :])), ρ[i, :], label="t = $(round(t[i], digits=2))", xlabel="x", ylabel="Density", title="Density Evolution")
    xlims!(0, 100)
    ylims!(1, 2.5)
end


println("\nAnimation complete! Saving to file...")
gif(anim, "Animation/density_evolution.gif", fps = round(nr_frames / 10 / gif_length, digits=0))

