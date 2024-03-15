# configuration

struct Scale{T<:Real}
    position::T       
    orientation::T  
    torsion::T
end


struct ConfSize
    ligands::Real   # AutoDock typed this szv (number of elements in a vector)
    flex::Real      # szv
    degrees_of_freedom # = sum(ligands) + sum(flex) + 6*length(ligands) # ?
end

# ...

struct RigidChange{T<:AbstractVector}
    position
    orientation
end

struct RigidConf
    position::AbstractVector{Float64}
    orientation::Quaternion
    RigidConf(position, orientation) = new(position, orientation)
end

RigidConf() = RigidConf(
    SVector(0,0,0), 
    Quaternion(1,0,0,0)     # identity Quaternion
)

# ...

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