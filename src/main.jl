using Printf
using Oceananigans
using Oceananigans.Units


Δt=3600seconds
stop_time=48hours

# 1 degree
nlat = 120
nlon = 60
nz = 40
latitude_bounds = range(-60, 60, length=nlat+1) |> collect
longitude_bounds = range(0, 60, length=nlon+1) |> collect
z_bounds = range(-4kilometers, 0, length=nz+1) |> collect

# %%

@printf("Grid generation...\n")
grid = LatitudeLongitudeGrid(CPU(); size=(nlon, nlat, nz), halo=(3, 3, 3), topology=(Bounded, Bounded, Bounded), latitude=latitude_bounds, longitude=longitude_bounds, z=z_bounds)

@printf("Model creation\n")
model = HydrostaticFreeSurfaceModel(grid;
                                    buoyancy = BuoyancyTracer(),
                                    tracers = :b,
                                    momentum_advection = WENO(),
                                    tracer_advection = WENO())

ϵ(x, y, z) = 2rand() - 1
set!(model, u=ϵ, v=ϵ)
simulation = Simulation(model; Δt=Δt, stop_time=stop_time)

@printf("Run model\n")
run!(simulation)

# %%

@printf("Plotting...\n")
using CairoMakie  # paper quality, no X11
#using GLMakie     # with X11
#GLMakie.activate!()

fig = Figure(size = (700, 200))
ax = Axis(fig[1, 1],
          xlabel="x [m]",
          ylabel="z [m]",
          limits=((0, grid.Lx), (-grid.Lz, 0)))

#band!(ax, x, bottom_boundary, top_boundary, color = :mediumblue)
plot!(ax, x, bottom_boundary, color = :mediumblue)
plot!(ax, x, top_boundary, color = :mediumblue)
save("figure.pdf", fig)

