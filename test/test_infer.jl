## test the main entry point in Celeste: the `infer` function
import Celeste
import JLD

"""
test infer with a single (run, camcol, field).
This is basically just to make sure it runs at all.
"""
function test_infer_single()
    # very small patch of sky that turns out to have 4 sources.
    # We checked that this patch is in the given field.
    ra_range = (164.39, 164.41)
    dec_range = (39.11, 39.13)
    fieldids = [(3900, 6, 269)]
    dirs = [datadir]

    result = Celeste.infer(ra_range, dec_range, fieldids, dirs)

    fname = @sprintf("%s/celeste-%.4f-%.4f-%.4f-%.4f.jld",
                     datadir, ramin, ramax, decmin, decmax) 
    JLD.save(fname, "result", result)
end

test_infer_single()