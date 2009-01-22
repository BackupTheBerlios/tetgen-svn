*******************************************
-- The Error Message
*******************************************
    Located (4) tet (3001, 5730, 4566, -1).
    Size of the cavity: 22 faces 22 tets.
    Insert point 3479
    Walk distance (# tets): 6
    Located (4) tet (2375, 4566, 2316, -1).
    Size of the cavity: 12 faces 8 tets.
    Insert point 3478
    Walk distance (# tets): 7
    Located (4) tet (2375, 4253, 4566, -1).
    Size of the cavity: 20 faces 14 tets.
    Insert point 5569
    Walk distance (# tets): 7
    Located (4) tet (3479, 4566, 4188, -1).
    Size of the cavity: 10 faces 6 tets.
    Insert point 5102
    Walk distance (# tets): 4
    Located (4) tet (5569, 3478, 4566, -1).
    Size of the cavity: 16 faces 10 tets.
    Insert point 5267
    Walk distance (# tets): 6
    Located (4) tet (5569, 4566, 4188, -1).
    Size of the cavity: 10 faces 6 tets.
    Insert point 5978
    Walk distance (# tets): 4
    Located (4) tet (5267, 5102, 4566, -1).
    Size of the cavity: 8 faces 4 tets.
    Insert point 224
    Walk distance (# tets): 5
    Located (4) tet (5978, 5102, 4566, -1).
    Size of the cavity: 10 faces 4 tets.
    Insert point 2969
    Walk distance (# tets): 5
    Located (4) tet (224, 4566, 5978, -1).
    Size of the cavity: 8 faces 3 tets.
    Scout subface (3017, 6071, 6069) (13620).
    Scout subface (6073, 3017, 6069) (13620).
    Scout subface (5343, 6073, 6069) (13620).
    Scout subface (3025, 5343, 6069) (13620).
    Scout subface (6071, 3025, 6069) (13620).
    Scout subface (6071, 3115, 3025) (13620).
    Scout subface (3001, 5343, 3025) (13620).
    Scout subface (5730, 3001, 3025) (13620).
constrain.cxx:1729: failed assertion `0'
Abort trap
tetgen:~/sync/trunk/data si$ open -e cavity.lua 
tetgen:~/sync/trunk/data si$ cp cavity.lua dump-cavity-case5.lua
tetgen:~/sync/trunk/data si$ mv dump-cavity-case4.lua dump-cavity-case5a.lua

*******************************************
-- How to reproduce this example
*******************************************

tetgen:~/sync/trunk/data si$ ../src/tetgen -pVV thomas.poly 

tetgen:~/sync/trunk/src si$ svn status -u   
sihang@svn.berlios.de's password: 
sihang@svn.berlios.de's password: 
?                   tetgen
?                   libtet.a
Status against revision:    324

