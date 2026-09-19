module CryoThermal

include("materials/Materials.jl")
include("network/Network.jl")
include("interfaces/ColdInterfaces.jl")
include("solvers/SteadyState.jl")

using .Materials
using .Network
using .ColdInterfaces
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

# Re-export Cold Interfaces (Indium Bolted Joint & Helium Exchange Gas)
export IndiumBoltedJoint,
       HeliumExchangeGasGap,
       IndiumContactLink,
       ExchangeGasLink,
       clamping_pressure,
       contact_conductance_indium,
       gas_mean_free_path,
       knudsen_number,
       gas_gap_conductance

# Re-export Solvers
export solve_steady_state, SteadyStateResult

end # module CryoThermal
