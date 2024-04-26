# configuration
export 
    Scale, 
    ConfSize, # is ConfSize necessary?
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
    Conf



struct Scale{T<:Real}
    position::T       
    orientation::T  
    torsion::T
end

Scale{T}() where T<:Real = Scale(zero(T), zero(T), zero(T))

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

ConfSize() = ConfSize(Int64[], Int64[])

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
        torsions[i] = normalize_angle(torsions[i])
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
    rconf.position = Quaternion{T}(1,0,0,0)
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
    in::Vector{Vector{T}},
    out::Vector{Vector{T}}
    #begin::Integer
    #end::Integer
) where T
    @assert size(in) == size(out)
    m = rotation_from_quaternion(rconf.orientation)
    for i in eachindex(in) # begin:end
        mul!(out[i], m, in[i]) + rconf.position
    end
end


struct LigandChange{T<:AbstractFloat}
    rigid::RigidChange{T}
    torsions::Vector{T}
end

LigandChange{T}() where T = LigandChange(RigidChange{T}(), T[])


struct LigandConf{T<:AbstractFloat}
    rigid::RigidConf{T}
    torsions::Vector{T}
end

LigandConf{T}() where T = LigandConf(RigidConf{T}(), T[])

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

ResidueChange{T}() where T = ResidueChange(T[])

struct ResidueConf{T<:AbstractFloat}
    torsions::Vector{T}
end

ResidueConf{T}() where T = ResidueConf(T[])

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


Change{T}() where T = Change(LigandChange{T}[], ResidueChange{T}[])


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


# AutoDock includes a pointer method as well as an identical reference method?
function (c::Change)(index::T) where T
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


function num_floats(c::Change)
    tmp = 0
    for i in eachindex(c.ligands)
        tmp += 6 + length(c.ligands[i].torsions)
    end
    for i in eachindex(c.flex)
        tmp += length(c.flex[i].torsions)
    end
    tmp
end

# does probably not need to be mutable, because we change elements inside of Vector
struct Conf{T<:AbstractFloat}
    ligands::Vector{LigandConf{T}}
    flex::Vector{ResidueConf{T}}
end


Conf{T}() where T = Conf(LigandConf{T}[], ResidueConf{T}[])


function Conf(s::ConfSize, type::T=Float64) where T
    # almost identical to Change
    ligands = Vector{LigandConf{type}}(undef, length(s.ligands))
    flex = Vector{ResidueConf{type}}(undef, length(s.flex))
    
    for i in eachindex(ligands)
        ligands[i] = LigandConf(RigidConf{type}(), fill(type(0), s.ligands[i]))
    end
    
    for i in eachindex(flex)
        flex[i] = ResidueConf(fill(type(0), s.flex[i]))
    end
    
    Conf(ligands, flex)
end


function set_to_null!(conf::Conf{T}) where T
    for i in eachindex(conf.ligands)
        set_to_null!(conf.ligands[i]) # set_to_null! called on each LigandConf
    end
    for i in eachindex(conf.flex)
        set_to_null!(conf.flex[i]) # set_to_null! called on each ResidueConf
    end
end


function increment!(conf::Conf{T}, c::Change{T}, factor::T) where T
    for i in eachindex(conf.ligands)
        increment!(conf.ligands[i], c.ligands[i], factor)
    end
    for i in eachindex(conf.flex)
        increment!(conf.flex[i], c.flex[i], factor)
    end
end


function internal_too_close(
    conf1::Conf{T}, 
    conf2::Conf{T}, 
    torsions_cutoff::T
) where T
    @assert length(conf1.ligands) == length(conf2.ligands)
    for i in eachindex(conf1.ligands) # could be written more elegently with all()
        if !torsions_too_close(
            conf1.ligands[i].torsions, 
            conf2.ligands[i].torsions, 
            torsions_cutoff
            )
            return false
        end
    end
    return true
end


function external_too_close(
    conf1::Conf{T}, 
    conf2::Conf{T}, 
    cutoff::Scale{T}
) where T
    @assert length(conf1.ligands) == length(conf2.ligands)
    for i in eachindex(conf1.ligands) # could be written more elegently with all()
        if !too_close(
            conf1.ligands[i].rigid, 
            conf2.ligands[i].rigid, 
            cutoff.position, 
            cutoff.orientation
            )
            return false
        end
    end
    @assert length(conf1.flex) == length(conf2.flex)
    for i in eachindex(conf1.flex)
        if !torsions_too_close(
            conf1.flex[i].torsions, 
            conf2.flex[i].torsions, 
            cutoff.torsions
            )
            return false
        end
    end
    return true
end


function too_close(conf1::Conf{T}, conf2::Conf{T}, cutoff::Scale{T}) where T
    (internal_too_close(conf1, conf2, cutoff.torsions) && 
     external_too_close(conf1, conf2, cutoff))
end


function generate_internal!(
    conf::Conf{T}, 
    torsion_spread::T, 
    rp::T, 
    rs::Union{Conf{T}, Nothing} = nothing   # rs does NOT need to exist
) where T
    for i in eachindex(conf.ligands)
        zero(conf.ligands[i].rigid.position)
        # check, if qt_identity works for different types
        conf.ligands[i].rigid.orientation = qt_identity(T) # Quaternion{T}(1,0,0,0)
        # check wether rs->ligands[i].rigid is referring to conf or rs!!!
        isnothing(rs) ? torsions_rs = rs : torsions_rs = rs.ligands[i].torsions
        torsions_generate!(conf.ligands[i].torsions, torsion_spread, rp, torsions_rs)
    end
end


function generate_external!(
    conf::Conf{T},
    spread::Scale{T}, 
    rp::T,
    rs::Union{Conf{T}, Nothing} = nothing
) where T
    for i in eachindex(conf.ligands)
        isnothing(rs) ? rigid_conf_rs = rs : rigid_conf_rs = rs.ligands[i].rigid
        generate!(
            conf.ligands[i].rigid, 
            spread.position, 
            spread.orientation, 
            rp, 
            rigid_conf_rs
        )
    end
    for i in eachindex(conf.flex)
        isnothing(rs) ? torsions_rs = rs : torsions_rs = rs.flex[i].torsions
        torsions_generate!(
            conf.flex[i].torsions, 
            spread.torsion, 
            rp, 
            torsions_rs
        )
    end
end


function randomize!(
    conf::Conf{T}, 
    corner1::Vector{T}, 
    corner2::Vector{T}
) where T
    for i in eachindex(conf.ligands)
        randomize!(conf.ligands[i], corner1, corner2)
    end
    for i in eachindex(conf.flex)
        randomize!(conf.flex[i])
    end
end


struct OutputType{T}
    c::Conf{T}
    # what are all these member variables?
    e::T
    lb::T
    ub::T
    intra::T
    inter::T
    conf_independent::T
    unbound::T
    total::T
    coords::Vector{T}
end


function OutputType{T}() where T
    OutputType(
        Conf{T}(), 
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        zero(T),
        T[]
    )
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