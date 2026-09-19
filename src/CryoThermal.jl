module CryoThermal

include("materials/Materials.jl")
include("network/Network.jl")
include("solvers/SteadyState.jl")

using .Materials
using .Network
using .SteadyState

# Re-export Material types and functions
export AbstractMaterial,
       CopperOFHC,
       StainlessSteel304,
       BerylliumCopper,
       Beryllium,
       thermal_conductivity,
       specific_heat,
       density,
       youngs_modulus,
       thermal_expansion,
       thermal_conductivity_integral

# Re-export Network types and functions
export ThermalNode,
       AbstractThermalLink,
       ConductionLink,
       ContactLink,
       RadiationLink,
       GenericResistorLink,
       ThermalSystem,
       heat_flow,
       stiffness_axial,
       energy_balance_residuals

# Re-export Solvers
export solve_steady_state, SteadyStateResult

end # module CryoThermal
