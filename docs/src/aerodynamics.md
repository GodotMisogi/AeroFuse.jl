# Aerodynamics API

## Doublet-Source Panel API

```@autodocs
Modules = [AeroFuse.DoubletSource]
```

## Coupled Potential-Flow API

`PotentialFlowSystem` solves named components built from populated potential-flow
elements. A component contains one concrete element type; different components
can use horseshoe vortices, vortex rings, source panels, or fuselage lines.
Fuselage-line components require at least two strictly increasing axial stations.

```julia
aircraft = ComponentVector(
    wing = elements(wing_mesh, Horseshoe()),
    body = elements(fuselage, SourcePanel()),
)
system = PotentialFlowSystem(aircraft, freestream, references)
wing_circulation = system.strengths.wing
body_source_density = system.strengths.body
```

The canonical module is `AeroFuse.PotentialFlow`. `AeroFuse.VortexLattice`,
`VortexLatticeSystem`, and `AbstractVortex` remain compatibility aliases.
The old properties `system.vortices` and `system.circulations` remain readable
as aliases for `system.elements` and `system.strengths`. Direct `getfield`
access, serialized type identities, and tools that use field names must migrate.

`AbstractPotentialFlowElement` is available for custom populated element types.
Its docstring describes the solve, transformation, wake, and force interfaces.
Strength units depend on the element type; only vortex strengths are circulation.

`surface_forces`, `surface_moments`, and `surface_dynamics` use the same
element-specific force models. `surface_velocities` samples trailing-induced
flow, freestream, rotation, and propeller slipstream at bound-leg centres
(control points for non-vortex adapters). It excludes bound-vortex induction
and prescribed fuselage thickness-source flow.

Streamlines and body pressure evaluation use solved-element induction but
exclude prescribed fuselage thickness sources and propeller slipstreams.
Matrix-free residuals apply the standard normal-velocity condition; they do
not replay prescribed fields or custom boundary-row overrides. For a coupled
linear solve, check the stored influence matrix and boundary vector with the
normalized strengths. Compressible strength recovery for non-vortex elements
has not been independently validated.

```@docs
AeroFuse.elements
```

```@autodocs
Modules = [AeroFuse.PotentialFlow]
```
