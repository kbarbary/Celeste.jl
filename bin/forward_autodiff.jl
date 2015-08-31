using Celeste
using CelesteTypes

using DataFrames
using SampleData

using ForwardDiff
using DualNumbers
import Transform
import Optim
import JLD

import NLsolve # Try the trust region method here?

include("src/interpolating_linesearch.jl")
include("src/NewtonsMethod.jl")

# Note that the u hessians are no good.
#omitted_ids = Int64[ids_free.u, ids_free.k[:], ids_free.c2[:], ids_free.r2];
#omitted_ids = ids_free.u;

galaxy_ids = union(ids_free.c1[:,2],
                   ids_free.c2[:,2],
                   ids_free.r1[2],
                   ids_free.r2[2],
                   ids_free.k[:,2],
                   ids_free.e_dev, ids_free.e_axis, ids_free.e_angle, ids_free.e_scale);

star_ids = union(ids_free.c1[:,1],
                   ids_free.c2[:,1],
                   ids_free.r1[1],
                   ids_free.r2[1],
                   ids_free.k[:,1]);


jld_file = "$dat_dir/SDSS_blob.jld"


simulation = true
if simulation
    blob, mp_original, body = gen_sample_star_dataset(perturb=false);
    #blob, mp_original, body = gen_sample_galaxy_dataset(perturb=true);
    #blob, mp_original, body = gen_three_body_dataset(perturb=true); # Too slow.
    transform = Transform.get_mp_transform(mp_original, loc_width=1.0)
else
    # An actual celestial body.
    field_dir = joinpath(dat_dir, "sample_field")
    run_num = "003900"
    camcol_num = "6"
    field_num = "0269"

    original_blob =
      Images.load_sdss_blob(field_dir, run_num, camcol_num, field_num,
                            mask_planes=Set());
    # Can't write a WCS pointer to a JLD file.
    #JLD.save(jld_file, "original_blob", original_blob)

    # Need to do this until WCS has an actual deep copy.
    original_crpix_band = Float64[unsafe_load(original_blob[b].wcs.crpix, i) for i=1:2, b=1:5];
    function reset_crpix!(blob)
        for b=1:5
            unsafe_store!(blob[b].wcs.crpix, original_crpix_band[1, b], 1)
            unsafe_store!(blob[b].wcs.crpix, original_crpix_band[2, b], 2)
        end
    end

    original_cat_df = SDSS.load_catalog_df(field_dir, run_num, camcol_num, field_num);
    cat_loc = convert(Array{Float64}, original_cat_df[[:ra, :dec]]);

    obj_cols = [:objid, :is_star, :is_gal, :psfflux_r, :compflux_r, :ra, :dec];
    sort(original_cat_df[original_cat_df[:is_gal] .== true, obj_cols], cols=:compflux_r, rev=true)
    sort(original_cat_df[original_cat_df[:is_gal] .== false, obj_cols], cols=:psfflux_r, rev=true)

    objid = "1237662226208063541" # A bright star with bad pixels
    #objid = "1237662226208063576" # A galaxy
    #objid = "1237662226208063565" # A brightish star but with good pixels.
    obj_row = original_cat_df[:objid] .== objid;
    obj_loc = Float64[original_cat_df[obj_row, :ra][1], original_cat_df[obj_row, :dec][1]]

    blob = deepcopy(original_blob);
    reset_crpix!(blob)
    WCS.world_to_pixel(blob[3].wcs, obj_loc)
    width = 15.
    Images.crop_image!(blob, width, obj_loc);
    @assert Images.test_catalog_entry_in_image(blob, obj_loc)
    entry_in_image = [Images.test_catalog_entry_in_image(blob, cat_loc[i,:][:]) for i=1:size(cat_loc, 1)];
    original_cat_df[entry_in_image, obj_cols]
    cat_entries = Images.convert_catalog_to_celeste(original_cat_df[entry_in_image, :], blob)
    mp_original = ModelInit.cat_init(cat_entries, patch_radius=20.0, tile_width=5);
    transform = Transform.get_mp_transform(mp_original)
end

function fit_only_type!(obj_type::Symbol, mp::ModelParams)
  valid_types = Symbol[:star, :galaxy, :both, :a]
  if !any(obj_type .== valid_types)
    error("obj_type must be in $(valid_types)")
  end
  epsilon = 0.006
  if obj_type == :star
    for s=1:mp.S
        mp.vp[s][ids.a] = [ 1.0 - epsilon, epsilon ]
    end
    omitted_ids = sort(unique(union(galaxy_ids, ids_free.a, ids_free.u)));
  elseif obj_type == :galaxy
    for s=1:mp.S
        mp.vp[s][ids.a] = [ epsilon, 1.0 - epsilon ]
    end
    omitted_ids = sort(unique(union(star_ids, ids_free.a, ids_free.u)));
  elseif obj_type == :both
    omitted_ids = Int64[];
  elseif obj_type == :a
    omitted_ids = setdiff(1:length(ids_free), ids_free.a)
    for s=1:mp.S
        mp.vp[s][ids.a] = [ 0.5, 0.5 ]
    end
  else
    error("obj_type must be in $(valid_types)")
  end
  omitted_ids
end


##############
# Get a BFGS fit for comparison
iter_count = NaN
bfgs_v = NaN

function bfgs_fit_params(mp_original::ModelParams, omitted_ids::Array{Int64})
  mp_bfgs = deepcopy(mp_original);
  iter_count, max_f, max_x, ret =
    OptimizeElbo.maximize_f(ElboDeriv.elbo, blob, mp_bfgs, transform,
                            omitted_ids=omitted_ids, verbose=true);
  mp_bfgs, iter_count, max_f
end
bfgs_v = ElboDeriv.elbo(blob, mp_bfgs).v;


####################
# Newton's method with our own hessian regularization

# For newton's method.
max_iters = 20;

function newton_fit_params(mp_original::ModelParams, omitted_ids::Array{Int64})
  mp_optim = deepcopy(mp_original);
  iter_count, max_f, max_x, ret =
    maximize_f_newton(mp -> ElboDeriv.elbo(blob, mp), mp_optim, transform,
                      omitted_ids=omitted_ids, verbose=true, max_iters=max_iters,
                      hess_reg=10.0)
  mp_optim, iter_count, max_f
end


function fit_type(obj_type::Symbol, mp_original::ModelParams, fit_fun::Function)
  mp_type = deepcopy(mp_original)
  omitted_ids = fit_only_type!(obj_type, mp_type);
  mp_type, iter_count, max_f = fit_fun(mp_type, omitted_ids)
  mp_type, iter_count, max_f
end

function combine_star_gal(mp_star::ModelParams, mp_gal::ModelParams)
  mp_combined = deepcopy(mp_original);
  for s=1:mp_combined.S
    mp_combined.vp[s][galaxy_ids] = mp_gal.vp[s][galaxy_ids]
    mp_combined.vp[s][star_ids] = mp_gal.vp[s][star_ids]
    mp_combined.vp[s][ids.a] = [0.5, 0.5]
  end
  mp_combined
end

# Try fitting one type at a time.
bfgs_results = Dict()
nm_results = Dict()

nm_results[:star], nm_star_iters, nm_star_v = fit_type(:star, mp_original, newton_fit_params)
bfgs_results[:star], bfgs_star_iters, bfgs_star_v = fit_type(:star, mp_original, bfgs_fit_params)

nm_results[:galaxy], nm_gal_iters, nm_gal_v = fit_type(:galaxy, mp_original, newton_fit_params)
bfgs_results[:galaxy], bfgs_gal_iters, bfgs_gal_v = fit_type(:galaxy, mp_original, bfgs_fit_params)

println(nm_star_v, ", ", nm_star_iters)
println(bfgs_star_v, ", ", bfgs_star_iters)

println(nm_gal_v, ", ", nm_gal_iters)
println(bfgs_gal_v, ", ", bfgs_gal_iters)

nm_results[:combined] = combine_star_gal(nm_results[:star], nm_results[:galaxy])
bfgs_results[:combined] = combine_star_gal(bfgs_results[:star], bfgs_results[:galaxy])
print_params(nm_results[:combined], bfgs_results[:combined])
ElboDeriv.get_brightness(nm_results[:combined])
ElboDeriv.get_brightness(bfgs_results[:combined])

nm_results[:a] = fit_type(:a, nm_results[:combined], newton_fit_params)[1]
bfgs_results[:a] = fit_type(:a, bfgs_results[:combined], bfgs_fit_params)[1]
print_params(nm_results[:a], bfgs_results[:a])

# For some reason this sets a to be 0.5...
nm_results[:both], nm_both_iters, nm_both_v =
fit_type(:both, nm_results[:a], newton_fit_params)
bfgs_results[:both], nm_both_iters, nm_both_v =
  fit_type(:a, bfgs_results[:a], bfgs_fit_params)



######################
# Print results
print_params(mp_original, mp_optim, mp_bfgs)

ElboDeriv.get_brightness(mp_original)[1]
ElboDeriv.get_brightness(mp_optim)[1]
ElboDeriv.get_brightness(mp_bfgs)[1]


##########################
# Simpler tests

function verify_sample_star(vs, pos)
    @test vs[ids.a[2]] <= 0.011

    @test_approx_eq_eps vs[ids.u[1]] pos[1] 0.1
    @test_approx_eq_eps vs[ids.u[2]] pos[2] 0.1

    brightness_hat = vs[ids.r1[1]] * vs[ids.r2[1]]
    @test_approx_eq_eps brightness_hat / sample_star_fluxes[3] 1. 0.01

    true_colors = log(sample_star_fluxes[2:5] ./ sample_star_fluxes[1:4])
    for b in 1:4
        @test_approx_eq_eps vs[ids.c1[b, 1]] true_colors[b] 0.2
    end
end

function verify_sample_galaxy(vs, pos)
    @test vs[ids.a[2]] >= 0.98

    @test_approx_eq_eps vs[ids.u[1]] pos[1] 0.1
    @test_approx_eq_eps vs[ids.u[2]] pos[2] 0.1

    @test_approx_eq_eps vs[ids.e_axis] .7 0.05
    @test_approx_eq_eps vs[ids.e_dev] 0.1 0.08
    @test_approx_eq_eps vs[ids.e_scale] 4. 0.2

    phi_hat = vs[ids.e_angle]
    phi_hat -= floor(phi_hat / pi) * pi
    five_deg = 5 * pi/180
    @test_approx_eq_eps phi_hat pi/4 five_deg

    brightness_hat = vs[ids.r1[2]] * vs[ids.r2[2]]
    @test_approx_eq_eps brightness_hat / sample_galaxy_fluxes[3] 1. 0.01

    true_colors = log(sample_galaxy_fluxes[2:5] ./ sample_galaxy_fluxes[1:4])
    for b in 1:4
        @test_approx_eq_eps vs[ids.c1[b, 2]] true_colors[b] 0.2
    end
end

omitted_ids = [ids_free.k[:], ids_free.c2[:], ids_free.r2]
blob, mp_original, body = gen_sample_star_dataset();

mp_bfgs = deepcopy(mp_original);
OptimizeElbo.maximize_likelihood(blob, mp_bfgs, transform);
verify_sample_star(mp_bfgs.vp[1], [10.1, 12.2]);

mp = deepcopy(mp_original);
iter_count, max_f, max_x, ret =
    OptimizeElbo.maximize_f_newton(mp -> ElboDeriv.elbo(blob, mp), mp, transform,
                               omitted_ids=omitted_ids, verbose=false, max_iters=10);
print_params(mp, mp_bfgs)
verify_sample_star(mp.vp[1], [10.1, 12.2])











##########################
#########################
# Newton's method by hand, probably obsolete.

obj_wrap = OptimizeElbo.ObjectiveWrapperFunctions(
    mp -> ElboDeriv.elbo(blob, mp), deepcopy(mp_original), transform, kept_ids, omitted_ids);
obj_wrap.state.scale = -1.0 # For minimization, which is required by the linesearch algorithm.

x0 = transform.vp_to_vector(mp_original.vp, omitted_ids);
elbo_grad = zeros(Float64, length(x0));
elbo_hess = zeros(Float64, length(x0), length(x0));


function f_grad!(x, grad)
  grad[:] = obj_wrap.f_grad(x)
end
d = Optim.DifferentiableFunction(obj_wrap.f_value, f_grad!, obj_wrap.f_value_grad!);

f_vals = zeros(Float64, max_iters);
cumulative_iters = zeros(Int64, max_iters);
x_vals = [ zeros(Float64, length(x0)) for iter=1:max_iters ];
grads = [ zeros(Float64, length(x0)) for iter=1:max_iters ];
hesses = [ zeros(Float64, length(x0), length(x0)) for iter=1:max_iters ];

# warm start with BFGS
warm_start = false
if warm_start
    mp_start = deepcopy(mp_original)
    start_iter_count, start_f, x_start = OptimizeElbo.maximize_f(ElboDeriv.elbo, blob, mp_start, Transform.free_transform, omitted_ids=omitted_ids, ftol_abs=1);
    obj_wrap.state.f_evals = start_iter_count;
    x_new = deepcopy(x_start); # For quick restarts while debugging
    new_val = old_val = -start_f;
else
    x_new = transform.vp_to_vector(mp_original.vp, omitted_ids);
    obj_wrap.state.f_evals = 0
    new_val = old_val = obj_wrap.f_value(x_new);
end

for iter in 1:max_iters
    println("-------------------$iter")
    x_old = deepcopy(x_new);
    old_val = new_val;

    elbo_hess = obj_wrap.f_ad_hessian(x_new);
    hesses[iter] = elbo_hess
    hess_ev = eig(elbo_hess)[1]
    min_ev = minimum(hess_ev)
    max_ev = maximum(hess_ev)
    println("========= Eigenvalues: $(max_ev), $(min_ev)")
    if min_ev < 0
        println("========== Warning -- non-convex, $(min_ev)")
        elbo_hess += eye(length(x_new)) * abs(min_ev) * 2
        hess_ev = eig(elbo_hess)[1]
        min_ev = minimum(hess_ev)
        max_ev = maximum(hess_ev)
        println("========= New eigenvalues: $(max_ev), $(min_ev)")
    end
    # if abs(max_ev) / abs(min_ev) > 1e6
    #     println("Regularizing hessian")
    #     elbo_hess += eye(length(x_new)) * (abs(max_ev) / 1e6)
    #     hess_ev = eig(elbo_hess)[1]
    #     min_ev = minimum(hess_ev)
    #     max_ev = maximum(hess_ev)
    #     println("========= New eigenvalues: $(max_ev), $(min_ev)")
    # end
    f_val, gr_new = obj_wrap.f_value_grad(x_old);
    x_direction = -(elbo_hess \ gr_new)

    lsr = Optim.LineSearchResults(Float64); # Not used
    c = -1.; # Not used
    mayterminate = true; # Not used
    pre_linesearch_iters = obj_wrap.state.f_evals
    interpolating_linesearch!(d, x_old, x_direction,
                              x_new, gr_new,
                              lsr, c, mayterminate;
                              c1 = 1e-4,
                              c2 = 0.9,
                              rho = 2.0, verbose=false);
    new_val, gr_new = obj_wrap.f_value_grad(x_new)
    println("Spent $(obj_wrap.state.f_evals - pre_linesearch_iters) iterations on linesearch for an extra $(f_val - new_val).")
    val_diff = new_val / old_val - 1
    f_vals[iter] = new_val;
    x_vals[iter] = deepcopy(x_new)
    grads[iter] = deepcopy(gr_new)
    cumulative_iters[iter] = obj_wrap.state.f_evals
    println(">>>>>>  Current value after $(obj_wrap.state.f_evals) evaluations: $(new_val) (BFGS got $(-bfgs_v) in $(iter_count) iters)")
    mp_nm = deepcopy(mp_original);
    transform.vector_to_vp!(x_new, mp_nm.vp, omitted_ids);
    #println(ElboDeriv.get_brightness(mp_nm))
    println("\n\n")
end

# f_vals are negative because it's minimization
println("Newton objective - BFGS objective (higher is better)")
println("Cumulative fuction evaluation ratio:")
hcat(((-f_vals) - bfgs_v) / abs(bfgs_v), cumulative_iters ./ iter_count)

reduce(hcat, [ x_diff ./ x_vals[1] for x_diff in diff(x_vals) ])

mp_nm = deepcopy(mp_original);
transform.vector_to_vp!(x_new, mp_nm.vp, omitted_ids);