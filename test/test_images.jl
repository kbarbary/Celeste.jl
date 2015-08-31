using Celeste
using Base.Test
using CelesteTypes
using SampleData
using DataFrames

import ModelInit
import Images
import SloanDigitalSkySurvey: SDSS

println("Running Images tests.")

field_dir = joinpath(dat_dir, "sample_field")
run_num = "003900"
camcol_num = "6"
field_num = "0269"

function test_blob()
  # A lot of tests are in a single function to avoid having to reload
  # the full image multiple times.

  blob = Images.load_sdss_blob(field_dir, run_num, camcol_num, field_num);
  cat_df = SDSS.load_catalog_df(field_dir, run_num, camcol_num, field_num);
  cat_entries = Images.convert_catalog_to_celeste(cat_df, blob);

  # Just check some basic facts about the catalog.
  @test size(cat_df)[1] == 805
  @test length(cat_entries) == size(cat_df)[1]
  @test sum([ cat_entry.is_star for cat_entry in cat_entries ]) ==
        sum(cat_df[:is_star])

  # Find an object near the middle of the image.
  #img_center = Float64[ median(cat_df[:ra]), median(cat_df[:dec]) ]
  img_center = WCS.pixel_to_world(blob[3].wcs, Float64[blob[3].H / 2, blob[3].W / 2])
  dist = by(cat_df, :objid,
     df -> DataFrame(dist=(df[:ra] - img_center[1]).^2 +
                          (df[:dec] - img_center[2]).^2))
  obj_rows = dist[:dist] .== minimum(dist[:dist])
  @assert sum(obj_rows) == 1
  obj_loc = Float64[ cat_df[obj_rows, :ra][1], cat_df[obj_rows, :dec][1]]
  objid = cat_df[obj_rows, :objid][1]

  # # Test cropping.  TODO -- delete this, it messes with the original center.
  # cropped_blob = deepcopy(blob)
  # width = 5.0
  # Images.crop_image!(cropped_blob, width, obj_loc)
  # for b=1:5
  #   @test 2 * width <= cropped_blob[b].H <= 2 * (width + 1)
  #   @test 2 * width <= cropped_blob[b].W <= 2 * (width + 1)
  # end
  # @test Images.test_catalog_entry_in_image(cropped_blob, obj_loc)

  # Test get_source_psf at point while we have the blob loaded.
  test_b = 3
  img = blob[test_b];
  obj_index = find(obj_rows)
  mp = ModelInit.cat_init(cat_entries[obj_index]);
  pixel_loc = WCS.world_to_pixel(img.wcs, obj_loc);
  original_psf_val =
    PSF.get_psf_at_point(pixel_loc[1], pixel_loc[2], img.raw_psf_comp);
  original_psf_gmm, scale = PSF.fit_psf_gaussians(original_psf_val);
  original_psf_celeste = Images.convert_gmm_to_celeste(original_psf_gmm, scale);
  fit_original_psf_val = PSF.get_psf_at_point(original_psf_celeste);

  obj_psf = Images.get_source_psf(mp.vp[1], img);
  obj_psf_val = PSF.get_psf_at_point(obj_psf);

  # The fits should match exactly.
  @test_approx_eq_eps(obj_psf_val[:], fit_original_psf_val[:], 1e-6)

  # The raw psf will not be as good.
  @test_approx_eq_eps(obj_psf_val[:], original_psf_val[:], 1e-2)

  mp_several =
    ModelInit.cat_init([cat_entries[1], cat_entries[obj_index]]);
  Images.set_patch_psfs!(blob, mp_several);

  # The second set of vp is the object of interest
  point_patch_psf = PSF.get_psf_at_point(mp_several.patches[2, test_b].psf);
  @test_approx_eq_eps(obj_psf_val[:], point_patch_psf[:], 1e-6)
end


function test_stamp_get_object_psf()
  stamp_blob, stamp_mp, body = gen_sample_star_dataset();
  img = stamp_blob[3];
  obj_loc =  stamp_mp.vp[1][ids.u]
  pixel_loc = WCS.world_to_pixel(img.wcs, obj_loc)
  original_psf_val = PSF.get_psf_at_point(img.psf);

  obj_psf_val =
    PSF.get_psf_at_point(Images.get_source_psf(stamp_mp.vp[1], img))
  @test_approx_eq_eps(obj_psf_val[:], original_psf_val[:], 1e-6)
end


function test_get_tiled_image_source()
  # Test that an object only occurs the appropriate tile's local sources.
  blob, mp, body = gen_sample_star_dataset();
  img = blob[3];

  [ mp.patches[1, b].radius = 1e-6 for b=1:5 ]

  tiled_img = Images.break_image_into_tiles(img, 10);
  for hh in 1:size(tiled_img)[1], ww in 1:size(tiled_img)[2]
    tile = tiled_img[hh, ww]
    loc = Float64[mean(tile.h_range), mean(tile.w_range)]
    [ mp.vp[1][ids.u] = mp.patches[1, b].center = loc for b=1:5 ]
    local_sources = ModelInit.get_tiled_image_sources(tiled_img, img.wcs, mp)
    @test local_sources[hh, ww] == Int64[1]
    for hh2 in 1:size(tiled_img)[1], ww2 in 1:size(tiled_img)[2]
      if (hh2 != hh) || (ww2 != ww)
        @test local_sources[hh2, ww2] == Int64[]
      end
    end
  end
end

test_blob()
test_stamp_get_object_psf()
test_get_tiled_image_source()
