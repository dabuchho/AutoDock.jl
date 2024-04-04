# configuration
export 
    Scale, 
    ConfSize, 
    torsions_set_to_null!, # could rename torsions_ methods because multiple dispatch
    torsions_increment!,
    torsions_randomize!,
    torsions_too_close,
    torsions_generate!,
    RigidChange,
    RigidConf,
    set_to_null!,
    increment!,
    randomize!,
    too_close,
    mutate_position!,
    mutate_orientation!,
    generate!,
    apply!,
    LigandChange,
    LigandConf,
    ResidueChange,
    ResidueConf,
    Change,
    operator,
    Conf



struct Scale{T<:AbstractFloat}
    position::T       
    orientation::T  
    torsion::T
end

# typedef std::size_t sz;           # => size of any object
# typedef std::vector<sz> szv;      # => vector of object sizes ? 
struct ConfSize              
    ligands::Vector{Int64}   # AutoDock typed this szv (number of elements in a vector)
    flex::Vector{Int64}      # szv
    degrees_of_freedom::Int64
    function ConfSize(ligands, flex)
        # 6 * number_of_atoms, because 3 translational and 3 rotational degrees
        # allow en empty ConfSize for now
        if length(ligands) == 0 || length(flex) == 0
            new(ligands, flex, 0)
        else
            new(ligands,
                flex,
                sum(ligands) + sum(flex) + 6*length(ligands)
            ) 
        end
    end
end         # this could be combined with BALL

# used in LigandConf and ResidueConf
function torsions_set_to_null!(torsions::Vector{T}) where T
    for i in eachindex(torsions)
        torsions[i] = 0.0
    end
end

# used in LigandConf and ResidueConf
function torsions_increment!(torsions::Vector{T}, c::Vector{T}, factor::T) where T
    for i in eachindex(torsions)
        torsions[i] += normalize_angle(factor * c[i])
        # AutoDock normalizes torsions[i] again ?
    end
end

# used in LigandConf and ResidueConf
function torsions_randomize!(torsions::Vector{T}) where T
    for i in eachindex(torsions)
        torsions[i] = rand(T) * 2π - π
    end
end     # TODO: add Seeding/ RNG option

# used in Conf
function torsions_too_close(
    torsions1::Vector{T},
    torsions2::Vector{T}, 
    cutoff::T
) where T
    @assert size(torsions1) == size(torsions2)
    for i in eachindex(torsions1)
        if abs(normalize_angle(torsions1[i] - torsions2[i])) > cutoff
            return false
        end
    end
    true
end

# used in Conf
function torsions_generate!(
    torsions::Vector{T}, 
    spread::T, 
    rp::T, 
    rs::Union{Vector{T}, Nothing} = nothing   # rs does not need to exist. TEST THIS!
) where T
    @assert isnothing(rs) || size(rs) == size(torsions) # if rs resent, check size
    for i in eachindex(torsions)
        if !isnothing(rs) && rand(T) < rp
            torsions[i] = rs[i]
        else
            torsions[i] += rand(T) * 2*spread - spread
        end
    end
end     # add Seed/ RNG option

struct RigidChange{T<:AbstractFloat}
    position::Vector{T}
    orientation::Vector{T}
end

RigidChange{T}() where T = RigidChange(T[0,0,0], T[0,0,0])


mutable struct RigidConf{T<:AbstractFloat}
    position::Vector{T}
    orientation::Quaternion{T}
end

RigidConf{T}() where T = RigidConf(
    Vector{T}([0,0,0]), 
    Quaternion{T}(1,0,0,0)     # identity Quaternion
)


function set_to_null!(rconf::RigidConf{T}) where T
    rconf.position = zeros(T, 3)
    rconf.position = zeros(T, 3)
end


function increment!(rconf::RigidConf{T}, rchan::RigidChange{T}, factor::T) where T
    rconf.position += factor * rchan.position
    rotation::Vector{T} = factor * rchan.orientation
    quaternion_increment(rconf.orientation, rotation)   # normalized Quaternion expected here
end


function randomize!(
    rconf::RigidConf{T}, 
    corner1::Vector{T}, 
    corner2::Vector{T}
) where T
    rconf.position = random_in_box(corner1, corner2)
    rconf.orientation = random_quaternion(T)
end


function too_close(
    rconf1::RigidConf{T}, 
    rconf2::RigidConf{T}, 
    position_cutoff::T, 
    orientation_cutoff::T
) where T
    if sqeuclidean(rconf1.position, rconf2.position) > sqr(position_cutoff)
        return false
    end
    if sqr(quaternion_difference(rconf1.orientation, rconf2.orientation)) > sqr(orientation_cutoff)
        return false
    end
    return true
end


function mutate_position!(rconf::RigidConf{T}, spread::T) where T
    rconf.position += spread * random_inside_sphere(spread)
end


function mutate_orientation!(rconf::RigidConf{T}, spread::T) where T
    tmp = spread * random_inside_sphere(T)
    quaternion_increment(rconf.orientation, tmp)
end


function generate!(
    rconf::RigidConf{T}, 
    position_spread::T, 
    orientation_spread::T,
    rp::T,
    rs::Union{RigidConf{T}, Nothing} = nothing  # might be slow
) where T
    if !isnothing(rs) && rand(T) < rp
        rconf.position = rs.position
    else
        mutate_position!(rconf.position, position_spread)
    end
    if !isnothing(rs) && rand(T) < rp
        rconf.orientation = rs.orientation
    else
        mutate_orientation!(rconf.orientation, orientation_spread)
    end
end

# apply orientation of RigidConf on in::Vector to receive out::Vector
function apply!(
    rconf::RigidConf{T}, 
    in::Vector{T}, # const
    out::Vector{T},
) where T
    @assert size(in) == size(out)
    m = rotation_from_quaternion(rconf.orientation)
    mul!(out, m, in) + rconf.position
end


struct LigandChange{T<:AbstractFloat}
    rigid::RigidChange{T}
    torsions::Vector{T}
end


struct LigandConf{T<:AbstractFloat}
    rigid::RigidConf{T}
    torsions::Vector{T}
end


function set_to_null!(lconf::LigandConf{T}) where T
    set_to_null!(lconf.rigid)
    torsions_set_to_null!(lconf.torsions)
end


function increment!(
    lconf::LigandConf{T}, 
    lchan::LigandChange{T},
    factor::T
) where T
    increment!(lconf.rigid, lchan.rigid, factor)
    torsions_increment!(lconf.torsions, lchan.torsions, factor)
end


function randomize!(
    lconf::LigandConf{T}, 
    corner1::Vector{T}, 
    corner2::Vector{T}
) where T
    randomize!(lconf.rigid, corner1, corner2)
    torsions_randomize!(lconf.torsions)
end


struct ResidueChange{T<:AbstractFloat}
    torsions::Vector{T}
end


struct ResidueConf{T<:AbstractFloat}
    torsions::Vector{T}
end


function set_to_null!(rconf::ResidueConf{T}) where T
    torsions_set_to_null!(rconf.torsions)
end


function increment!(
    rconf::ResidueConf{T}, 
    rchan::ResidueChange{T},
    factor::T
) where T
    torsions_increment!(rconf.torsions, rchan.torsions, factor)
end


function randomize!(rconf::ResidueConf{T}) where T
    torsions_randomize!(rconf.torsions)
end


struct Change{T<:AbstractFloat}
    ligands::Vector{LigandChange{T}}
    flex::Vector{ResidueChange{T}}
end


# This is super ugly. This constructor, when called, also always creates 
# a new LigandChange and RigidChange object, which is not ideal.
function Change(s::ConfSize, type::T=Float64) where T
    # first, initialize a vector of the size of the ConfSize
    ligands = Vector{LigandChange{type}}(undef, length(s.ligands))
    flex = Vector{ResidueChange{type}}(undef, length(s.flex))
    
    # then, resize ligands and flex according to each element in ConfSize.ligands
    # and Confsize.flex
    for i in eachindex(ligands)
        ligands[i] = LigandChange(RigidChange{type}(), fill(type(0), s.ligands[i]))
    end
    
    for i in eachindex(flex)
        flex[i] = ResidueChange(fill(type(0), s.flex[i]))
    end
    
    Change(ligands, flex)
end


function operator(c::Change, index::T) where T
    for i in eachindex(c.ligands)
        lig = c.ligands[i]
        if index < 3 return lig.rigid.position[index] end
        index -= 3
        if index < 3 return lig.rigid.position[index] end
        index -= 3
        if index < length(lig.torsions) return lig.torsions[index] end
        index -= length(lig.torsions)
    end
    for i in eachindex(c.flex)
        res = c.flex[i]
        if index < length(res.torsions) return res.torsions[index] end
        index -= length(res.torsions)
    end
end

struct Conf{T<:AbstractFloat}
    ligands::Vector{LigandChange{T}}
    flex::Vector{ResidueChange{T}}
end

# ...

struct OutputType
    c::Conf
    e
    lb
    ub
    intra
    inter
    conf_independent
    unbound
    total
    coords
end


# Vectors of position, orientation, torsions
pos = Float64[]
ori = Float64[]
tor = Float64[]
rchan = RigidChange{Float64}(pos, ori)
lchan = LigandChange{Float64}(rchan, tor)
vlc = Vector{LigandChange{Float64}}([lchan])
vrc = Vector{ResidueChange{Float64}}(tor)

# Change(vlc, vrc) # unfinished

# Vectors of ligands and flex
lig = [1,2,3,4,5]
flex = [1,1,1]
cs = ConfSize(lig, flex)

change = Change(cs)


Vector{LigandChange{Float64}}()