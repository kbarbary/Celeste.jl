# linear approximation to WCS transformations

module LinearWCS

export LinearWCSTransform2D

using WCS
import Base: convert, *
import WCS: world_to_pix

# We could use FixedSizeArrays here, but it's not a lot of functionality,
# so we'll just use plain NTuples and define a few things we need.

convert{N, T}(::Type{Vector{T}}, x::NTuple{N, T}) = [x...]
convert{N, T}(::Type{Matrix{T}}, x::NTuple{N, NTuple{N, T}}) =
    [x[j][i] for i=1:N, j=1:N]

# define matrix-vector multiplication for tuples
function *{T, S}(x::NTuple{2, NTuple{2, T}}, y::NTuple{2, S})
    (x[1][1] * y[1] + x[2][1] * y[2], x[1][2] * y[1] + x[2][2] * y[2])
end

immutable LinearWCSTransform2D{T}
    world_offset::NTuple{2, T}
    pix_offset::NTuple{2, T}
    jacobian::NTuple{2, NTuple{2, T}}  # transforms world -> pixel
end

function LinearWCSTransform2D(wcs::WCSTransform, pixcoords; dpix=0.5)
    wc = pix_to_world(wcs, convert(Vector{Float64}, pixcoords))

    # We'll determine the local rotation matrix via finite differences
    # in pixel space.
    pixcoords1 = [pixcoords[1] + dpix, pixcoords[2]       ]
    pixcoords2 = [pixcoords[1]       , pixcoords[2] + dpix]
    wc1 = pix_to_world(wcs, pixcoords1)
    wc2 = pix_to_world(wcs, pixcoords2)

    # Construct pixel-to-world rotation matrix and take the inverse
    # to get the world-to-pixel rotation matrix.
    J = inv([(wc1[1] - wc[1]) / dpix  (wc2[1] - wc[1]) / dpix;
             (wc1[2] - wc[2]) / dpix  (wc2[2] - wc[2]) / dpix])

    # convert J to a tuple
    Jtup = ((J[1, 1], J[2, 1]), (J[1, 2], J[2, 2]))

    return LinearWCSTransform2D((wc[1], wc[2]),
                                (pixcoords[1], pixcoords[2]),
                                Jtup)
    end


function world_to_pix(linwcs::LinearWCSTransform2D, worldcoords)
    pixoff = linwcs.jacobian * (worldcoords[1] - linwcs.world_offset[1],
                                worldcoords[2] - linwcs.world_offset[2])
    return (pixoff[1] + linwcs.pix_offset[1],
            pixoff[2] + linwcs.pix_offset[2])
end

end  # module
