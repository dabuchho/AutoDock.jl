# AutoDock

ad4cache                # switch case statments, mostly for atom types and errors
array3d                 # resizable 3d Array
atom                    # basically BALL Atom()
atom_base               
atom_constants          # data sheet and atom properties
atom_type               # atom_type struct 
bfgs                    # Broyden–Fletcher–Goldfarb–Shanno algorithm
brick                   # brick distances (?)
cache                   # many switch cases for caching + mapping
common                  # masic maths / linear algebra
conf                    # Atom configurations, torsions, ligands, residues
conf_independent        # similar to conf
convert_substring       # 3 important functions (rmsd, find_closest, to_container)
coords                  
curl                    # curl! function
file                    # in and out file handling
grid                    # grid class: linear algebra, voxels (= BALL System?)
grid_dim                # grid stuff
igrid                   # grid stuff
incrementable           
int_pow
macros  
matrix                  # matrix, triangular matrix
model                   # atom range, branch, appender, beads, bonds
                        # over 1000 lines of code. Eval how much is in BALL
monte_carlo             # monte carlo operator (solver?)
                        # https://qutip.org/docs/latest/guide/dynamics/dynamics-monte.html, 
                        # https://de.wikipedia.org/wiki/Metropolis-Algorithmus, 
                        # https://de.wikipedia.org/wiki/MCMC-Verfahren
mutate                  # mutating configuration randomly
non_cache               # many switch cases about atom types
parallel                # parallel computing (parallel for loops)
parallel_mc             # parallel monte carlo tasks
parallel_progress       
parse_error             
parse_pdbqt             # 600+ lines of code to parse AutoDock structure files
                        # https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/AutoDock4.2.6_UserGuide.pdf (page 28)
potentials              # Potentials class: here we could include NESSie
precalculate            # precalculations: described in paper
quasi_newton            # https://en.wikipedia.org/wiki/Quasi-Newton_method
quaternion              # quaternion logic
random                  # pseudo random numbers generators
scoring_function        # Scoring function (based on Potentials)
szv_grid                # uses a lot from model to create grid (BALL System?)
tree                    # atom structure? (tree is used in model.h)
triangular_matrix_index 
utils                   
vina                    # console output