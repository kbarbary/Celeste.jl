using Celeste: LinearWCS
using WCS

# check that we implemented vector-matrix multiplication correctly
m = rand(2, 2)
v = rand(2)
ans = m * v

mtup = ((m[1, 1], m[2, 1]), (m[1, 2], m[2, 2]))
vtup = (v...)
result = mtup * vtup
@assert result[1] == ans[1] && result[2] == ans[2]

# Some arbitrary transform
wcs = WCSTransform(2;
                   cdelt = [-0.066667, 0.066667],
                   ctype = ["RA---AIR", "DEC--AIR"],
                   crpix = [-234.75, 8.3393],
                   crval = [0., -90],
                   pv    = [(2, 1, 45.0)])

# approximate it at some location
pixcoords = [-250., 25.]
worldcoords = pix_to_world(wcs, pixcoords)
linwcs = LinearWCSTransform2D(wcs, pixcoords)

# check how good the approximation is 10 pixels away
testwc = pix_to_world(wcs, pixcoords .+ 1.0)
println(world_to_pix(wcs, testwc))
println(world_to_pix(linwcs, testwc))
