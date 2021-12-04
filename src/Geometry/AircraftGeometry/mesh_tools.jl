abstract type AbstractSpacing end

struct Cosine  <: AbstractSpacing end
struct Sine    <: AbstractSpacing end
struct Uniform <: AbstractSpacing end

spacing(n, :: Uniform, flip = false) = uniform_spacing(0., 1., n + 1)
spacing(n, :: Cosine, flip = false)  = cosine_spacing(0.5, 1., n + 1)
spacing(n, :: Sine, flip = false) = sine_spacing(0., 1., (n + 1) * ifelse(flip, -1, 1))
spacing(n, distribution :: AbstractSpacing = Cosine(), flip = false) = spacing(n, distribution, flip)

Base.reverse(q :: AbstractSpacing) = q
Base.Broadcast.broadcastable(q::AbstractSpacing) = Ref(q)

function chop_sections(set1, set2, n :: Integer, distribution = Cosine(); flip = false)
    space = spacing(n, distribution, flip)
    @views [ weighted_vector.(set1, set2, μ) for μ ∈ space ][1:end-1]
end

chop_coordinates(coords, n, spacing = Cosine(), flip = false) = @views [ reduce(vcat, chop_sections.(coords[1:end-1], coords[2:end], n, spacing; flip = flip)); [ coords[end] ] ]

chop_spans(xyzs, div, spacing = Cosine(), flip = false) = permutedims(combinedimsview(map(xyz -> chop_coordinates(xyz, div, spacing, flip), eachrow(xyzs))))

chop_chords(xyzs, div, spacing = Cosine()) = combinedimsview(map(xyz -> chop_coordinates(xyz, div, spacing), eachcol(xyzs)))

chop_wing(xyzs, span_num, chord_num; span_spacing = Cosine(), chord_spacing = Cosine(), flip = false) = chop_chords(chop_spans(xyzs, span_num, span_spacing, flip), chord_num, chord_spacing)

transform_coordinates(xyz, twist, section) = eachrow(xyz * RotY(-twist)') .+ Ref(section)

function chop_spanwise_sections(scaled_foils, twisties, leading_xyz, span_num, spacings, flip = false)
    # Reverse direction if left side
    if flip
        scaled_foils = reverse(scaled_foils)
        twisties     = reverse(twisties)
    end

    # Rotate and translate coordinates
    foil_coords = combinedimsview(transform_coordinates.(scaled_foils, twisties, leading_xyz))

    # Chop up spanwise sections
    chop_spans(foil_coords, span_num, spacings, flip)
end

function interpolate_grid(mesh, f_interp)
    
     mapslices(mesh, dims = 1:2) do 


        
     end

end