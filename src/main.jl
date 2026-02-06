using Printf
using Oceananigans
using Oceananigans.Units
using CUDA

Δt=3600seconds
stop_time=10days

# 1 degree
nlat = 120
nlon = 60
nz = 40
latitude_bounds = range(-60, 60, length=nlat+1) |> collect
longitude_bounds = range(0, 60, length=nlon+1) |> collect
z_bounds = range(-4kilometers, 0, length=nz+1) |> collect

# %%

@printf("Grid generation...\n")
grid = LatitudeLongitudeGrid(GPU(); size=(nlon, nlat, nz), halo=(3, 3, 3), topology=(Bounded, Bounded, Bounded), latitude=latitude_bounds, longitude=longitude_bounds, z=z_bounds)

@printf("Create boundary conditions\n")

const ρ0 = 1025.0,
const cp = 3999.0,
const τx = 0.1, # N/m^2 
const τy = 0.0,
const Q = 200.0,
const F = 1e-7, # m/s ( 3 m/year )
const S0 = 35.0

Q_T(x, y, t) = Q / (ρ0 * cp)
S_flux(x, y, t) = -F * S0
τu(x, y, t) = τx / ρ0
τv(x, y, t) = τy / ρ0

boundary_conditions = (
    T = FieldBoundaryConditions(
        top = FluxBoundaryCondition(Q_T)
    ),
    S = FieldBoundaryConditions(
        top = FluxBoundaryCondition(S_flux)
    ),
    u = FieldBoundaryConditions(
        top = FluxBoundaryCondition(τu)
    ),
    v = FieldBoundaryConditions(
        top = FluxBoundaryCondition(τv)
    ),
)

@printf("Model creation\n")
model = HydrostaticFreeSurfaceModel(
    grid;
    buoyancy = BuoyancyTracer(),
    tracers = :b,
    momentum_advection = WENO(),
    tracer_advection = WENO(),
    boundary_conditions = boundary_conditions,
)

@printf("Set initial conditions\n")
ϵ(x, y, z) = 2rand() - 1
set!(model, u=ϵ, v=ϵ)


@printf("Construct Simulation object\n")
simulation = Simulation(model; Δt=Δt, stop_time=stop_time)
progress(sim) = @printf("Iter: % 6d, sim time: % 1.3f, wall time: % 10s, Δt: % 1.4f, advective CFL: %.2e, diffusive CFL: %.2e\n",
                        iteration(sim), time(sim), prettytime(sim.run_wall_time),
                        sim.Δt, AdvectiveCFL(sim.Δt)(sim.model), DiffusiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))


@printf("Model information:\n")
display(model)
@printf("Simulation information:\n")
display(simulation)
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

