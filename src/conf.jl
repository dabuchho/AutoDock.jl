# configuration

struct Scale{T<:AbstractFloat}
    position::T       
    orientation::T  
    torsion::T
end

# typedef std::size_t sz;           # => size of any object
# typedef std::vector<sz> szv;      # => vector of object sizes ? 
struct ConfSize              # TODO: option between Int64 and Int32
    ligands::Vector{Int64}   # AutoDock typed this szv (number of elements in a vector)
    flex::Vector{Int64}      # szv
    degrees_of_freedom::Int64
    function ConfSize(ligands, flex)
        new(ligands,
            flex,
            sum(ligands) + sum(flex) + 6*length(ligands) # ?
        ) 
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
    rs::Vector{T}       # rs as kwargs... ?
) where T
    @assert size(rs) == size(torsions)    # rs not always present => add nothing condition?
    for i in eachindex(torsions)
        if rand(T) < rp                     # rs condition
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

RigidChange{T}() where T = RigidChange(T[0,0,0], T[0,0,0])  # is this bad practice?


mutable struct RigidConf{T<:AbstractFloat}
    position::Vector{T}
    orientation::Quaternion{T}
end

RigidConf{T}() where T = RigidConf(
    Vector{T}(0,0,0), 
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


struct LigandChange
    rigid::RigidChange
    torsions::AbstractVector
end

struct LigandConf
    rigid::RigidChange
    torsions::AbstractVector
end

# ...

struct ResidueChange
    torsions::AbstractVector
end

# ...

struct Change
    ligands::Vector{LigandChange}
    flex::Vector{ResidueChange}
end

# ...

struct Conf
    ligands::Vector{LigandChange}
    flex::Vector{ResidueChange}
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
