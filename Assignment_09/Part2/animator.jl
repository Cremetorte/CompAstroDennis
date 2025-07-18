using Plots
using Base.Threads


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


quantity = "internal_energy"
filename = "Output/output_$quantity.csv"

df = CSV.read(filename, DataFrame)
data = Matrix(df)

t = data[:, 1]
ρ = data[:, 2:end]

# println("Shape of data: ", size(data))

nr_frames = size(data)[1]
gif_length = 10 # seconds


ylim = maximum(ρ) * 1.1
ymin = minimum(ρ) * 0.9


# Create a GIF


frames = Vector{Plots.Plot}(undef, nr_frames)

Threads.@threads for i in 1:nr_frames
    frames[i] = plot(
        x.(eachindex(ρ[i, :])), ρ[i, :],
        label="t = $(round(t[i], digits=2))",
        xlabel="x", ylabel="$quantity",
        title="$quantity evolution",
        legend=:topright,
        ylim=(ymin, ylim), xlims=(0, 100)
    )
end

anim = @animate for i in 1:nr_frames
    print("Rendering frame $(i) / $(nr_frames)...   \r")
    plot!(frames[i])
end


println("\nAnimation complete! Saving to file...")
gif(anim, "Animation/anim_$quantity.gif", fps = round(nr_frames / gif_length, digits=0))

