using Quaternions

# example testing
# create H₂O Molecule


rand(Float64, 3)
zeros(Float64, 3)

torsions = Int64[1,2,3]
torsions_set_to_null!(torsions)
torsions == zeros(Int64, length(torsions))

torsions = Float64[10,2,3]
c = Float64[1,2,3]
factor = 1.0
# Conf torsions vector, torsion Changes vector, factor
torsions_increment!(torsions, c, factor)
torsions



# general testing
@testitem "Configurations" begin
    # Empty instances
    for T in [Float32, Float64]
        @test typeof(Scale{T}())         == Scale{T}
        @test typeof(ConfSize())         == ConfSize
        @test typeof(RigidChange{T}())   == RigidChange{T}
        @test typeof(RigidConf{T}())     == RigidConf{T}
        @test typeof(LigandChange{T}())  == LigandChange{T}
        @test typeof(LigandConf{T}())    == LigandConf{T}
        @test typeof(ResidueChange{T}()) == ResidueChange{T}
        @test typeof(ResidueConf{T}())   == ResidueConf{T}
        @test typeof(Change{T}())        == Change{T}
        @test typeof(Conf{T}())          == Conf{T}
    end
    
    # (methods that don't require composite types as input)
    # move this to LigandConf and ResidueConf tests
    for T in [Float32, Float64]
        n = rand(1:10)
        torsions = rand(T, n)
        torsions_set_to_null!(torsions)
        @test torsions == zeros(T, n) # LigandConf and ResidueConf
        #torsions_increment!()  # LigandConf and ResidueConf
        #torsions_randomize!()  # LigandConf and ResidueConf
        #torsions_too_close()   # Conf
        #torsions_generate!()   # Conf
    end

    # Scale
    Scale{Float64}(10, 10, 10) # only used for cutoffs. Set very high here

    # ConfSize
    ConfSize([3], [2]) # 3 Atoms? (ligands), 2 Bonds? (flex)

    for T in [Float32, Float64]
        # RigidChange
        rchan_position = rand(T, 3)    # always single atom position?
        rchan_orientation = rand(T, 3) # always single atom orientation?
        rchan = RigidChange(rchan_position, rchan_orientation)

        # RigidConf
        rconf_position = rand(T, 3)              # always single atom position?
        rconf_orientation = random_quaternion(T) # always single atom orientation?
        rconf = RigidConf(rconf_position, rconf_orientation)
        
        factor1 = T(1.0) # does not effect increment!
        increment!(rconf, rchan, factor1)
        @test rconf.position == rconf_position + rchan_position
        # TODO: find way to predict the effect of the rotation
        @test typeof(rconf.orientation) == Quaternion{T}
        # different factor
        # ...

        # randomize!()
        # too_close()
        # mutate_position!()
        # mutate_orientation!()
        # generate!()
        # apply!()

        set_to_null!(rconf)
        @test rconf.position    == zeros(T, 3)
        @test rconf.orientation == Quaternion{T}(1,0,0,0)
    end



    # Create a mockup molecule
    # ...   
end