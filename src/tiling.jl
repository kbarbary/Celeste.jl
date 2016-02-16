

"""
    circle_rectangle_intersect(xctr, yctr, r, xmin, xmax, ymin, ymax)

Determine if a circle and rectangle intersect. It is assumed that
`xmax` >= `xmin` and `ymax` >= `ymin`.
"""
function circle_rectangle_intersect(xctr, yctr, r, xmin, xmax, ymin, ymax)

    # shift coordinates so that circle center is at (0, 0)
    xmin -= xctr
    xmax -= xctr
    ymin -= yctr
    ymax -= yctr

    # There are 2^4 = 16 configurations for the signs
    # of (xmin, xmax, ymin, ymax). We assume the inputs obey xmax >= xmin
    # and ymax >= ymin, leaving 9 cases.

    # The case where coordinates are both negative is symmetric with
    # the case where they're both positive. This eliminates 5 more cases.
    if xmin < 0.0 && xmax < 0.0
        xmin = -xmax
        xmax = -xmin
    end
    if ymin < 0.0 && ymax < 0.0
        ymin = -ymax
        ymax = -ymin
    end

    # The remaining 4 cases for the signs of (xmin, xmax, ymin, ymax):
    if xmin > 0.0
        (ymin > 0.0) && return (xmin*xmin + ymin*ymin < r*r)  # (+, +, +, +)
        return (xmin < r)  # (+, +, -, +)
    else
        (ymin > 0.0) && return (ymin < r)  # (-, +, +, +)
        return true  # (-, +, -, +)
    end
end

"""
Return the range of image pixels in an ImageTile.

Args:
  - hh: The tile row index (in 1:number of tile rows)
  - ww: The tile column index (in 1:number of tile columns)
  - H: The number of pixel rows in the image
  - W: The number of pixel columns in the image
  - tile_width: The width and height of a tile in pixels
"""
function tile_range(hh::Int, ww::Int, H::Int, W::Int, tile_width::Int)
    h1 = 1 + (hh - 1) * tile_width
    h2 = min(hh * tile_width, H)
    w1 = 1 + (ww - 1) * tile_width
    w2 = min(ww * tile_width, W)
    h1:h2, w1:w2
end

"""
    tile_indicies(xctr, yctr, r, h, w, tile_width)

Return indicies of tiles that intersect the circle defined by
`xctr`, `yctr`, and `r`. Result is a `Vector{Tuple{Int, Int}}`
where each entry gives the coordinates of the tile.
"""
function tile_indicies(xctr, yctr, r, h, w, tile_width)

    # extent of pixel indicies in the image that circle could possibly
    # subtend (inclusive)
    xmin = max(round(Int, xctr - r), 1) 
    xmax = min(round(Int, xctr + r), h)
    ymin = max(round(Int, yctr - r), 1)
    ymax = min(round(Int, yctr + r), w)

    # extent of tile indicies that circle could possibly subtend
    # (inclusive)
    txmin = div(xmin-1, tile_width) + 1 
    txmax = div(xmax-1, tile_width) + 1
    tymin = div(ymin-1, tile_width) + 1
    tymax = div(ymax-1, tile_width) + 1

    # check whether each tile overlaps the circle; if so, append its index
    out = Tuple{Int, Int}[]
    for j=tymin:tymax, i=txmin:txmax
        hr, wr = tile_range(i, j, h, w, tile_width)
        if circle_rectangle_intersect(xctr, yctr, r, first(hr), last(hr),
                                      first(wr), last(wr))
            push!(out, (i, j))
        end
    end

    return out
end
