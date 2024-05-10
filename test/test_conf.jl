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
    using Quaternions
    using Distances: euclidean
    # Constructors and Type Conservation
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
    ##############################################################

    # Scale
    Scale{Float64}(10, 10, 10) # only used for cutoffs. Set very high here

    # ConfSize
    ConfSize([3], [2]) # 3 Atoms? (ligands), 2 Bonds? (flex)

    for T in [Float32, Float64]
        # RigidChange
        rig_chan_position = rand(T, 3)    # always single atom position?
        rig_chan_orientation = rand(T, 3) # always single atom orientation?
        rig_chan = RigidChange(rig_chan_position, rig_chan_orientation)

        # RigidConf
        rig_conf_position = rand(T, 3)              # always single atom position?
        rig_conf_orientation = random_quaternion(T) # always single atom orientation?
        rig_conf = RigidConf(rig_conf_position, rig_conf_orientation)
        
        factor1 = T(1.0) # does not effect increment!
        increment!(rig_conf, rig_chan, factor1)
        @test rig_conf.position == rig_conf_position + rig_chan_position
        # TODO: find way to predict the effect of the rotation
        @test typeof(rig_conf.orientation) == Quaternion{T}
        
        factor2 = T(2.0) # factor ≠ 1
        increment!(rig_conf, rig_chan, factor2)
        #@test rig_conf.position == rig_conf_position + rig_chan_position

        corner1 = T[0,0,0]
        corner2 = T[1,1,1]
        randomize!(rig_conf, corner1, corner2)
        @test euclidean(rig_conf.position, corner1) < T(sqrt(3))
        @test euclidean(rig_conf.position, corner2) < T(sqrt(3))
        #@test rig_conf.orientation # ...
        

        rig_conf1 = RigidConf{T}(T[0,0,0], Quaternion{T}(1,0,0,0))
        rig_conf2 = RigidConf{T}(T[1,0,0], Quaternion{T}(0,0,0,1))
        pos_cutoff = T(0)   # low cutoff
        ori_cutoff = T(100) # high cutoff
        @test too_close(rig_conf1, rig_conf2, pos_cutoff, ori_cutoff) == false
        pos_cutoff = T(100) # high cutoff
        ori_cutoff = T(0)   # low cutoff
        @test too_close(rig_conf1, rig_conf2, pos_cutoff, ori_cutoff) == false
        rig_conf3 = RigidConf{T}(T[0,0,0], Quaternion{T}(1,0,0,0))
        pos_cutoff = T(0)
        # too_close should return true for overlaps, independent of the cutoff
        @test too_close(rig_conf1, rig_conf3, pos_cutoff, ori_cutoff) == true
        
        current_pos = rig_conf1.position
        current_ori = rig_conf1.orientation
        spread = T(2.0)
        mutate_position!(rig_conf1, spread)
        @test current_pos != rig_conf1.position
        @test euclidean(current_pos, rig_conf1.position) <= spread
        # TODO: find better test
        mutate_orientation!(rig_conf1, spread)
        @test current_ori != rig_conf1.orientation # this should almost never fail
        # generate!()
        # apply!()

        set_to_null!(rig_conf)
        @test rig_conf.position    == zeros(T, 3)
        @test rig_conf.orientation == Quaternion{T}(1,0,0,0)

        # LigandChange
        # LigandConf

        # ResidueChange
        # ResidueConf

        # Change
        # Conf

        # OutputType

    end



    # Create a mockup molecule
    # ...   
end
using Quaternions
using AutoDock
quaternion_difference(Quaternion{Float64}(1.0,0,0,0),Quaternion{Float64}(1,0,0,0))

rc = RigidConf{Float64}()
old_ori = rc.orientation
mutate_orientation!(rc, 2.0)
curr_ori = rc.orientation