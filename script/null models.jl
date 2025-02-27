using ArchGDAL, Rasters, DataFrames, Plots, JLD2
using AvianRangeShapes


#loading geographical domain
dd = Raster("data/sf1_mainland.tif")  
dom_master = .!ismissing.(dd)

#loading topographical raster
top=resample(Raster("data/top_q_proj.tif"); to = dom_master)
replace!(top,missing=>NaN32)

#loading species elevational range limits
elv=load_object("data/elevational range limits.jld2")
ele_range = Dict(row.Species => (min = row.minimum_elevation, max = row.Maximu_elevation) for row in eachrow(elv))


#loading the geographic geographic ranges of six example species
dis_elv=load_object("data/bird ranges.jld2")
nam=["Acropternis orthonyx","Aglaeactis castelnaudii","Coeligena lutetiae","Diglossa mystacalis","Phlogophilus harterti","Scytalopus griseicollis"]
geo_range = Dict(zip(nam, dis_elv))

#loading standardized range sizes
formated_rs=load_object("data/standardized_range_sizes.jld2")

############################ run the null model 

example_species="Phlogophilus harterti" # name of one of the example species from the nam object
rs_std=false # should the null model use the standardized range size (true) or the empirical range size (false)
nrep=10 # nuber of repetitions

#### empirical range
map_emp = falses(dims(dom_master))
map_emp[geo_range[example_species]] .= true
map_emp=prep_map(geo_range[example_species],trim_map=true)
ex=ArchGDAL.extent(map_emp)
plot(map_emp)

#### running null model 1
res_nm1=null_models(example_species, # state name of example species
                    "nm1",# state which of the four null models (nm1,nm2,nm3, or nm4)
                    geo_range, #grid cell ids comprising the species empirical range
                    rs_std, # should the null model use the standardized range size (true) or the empirical range size (false)
                    copy(dom_master), #biogeographical domain
                    top, #topographical raster
                    elv, #data frame with the species' elevational range limits (only used in null model nm2 and nm4)
                    formated_rs, #data frame with standardized range sizes (only used if rs_std=true)
                    nrep # nuber of repetitions
                    ; patchy_ranges = false, bounded_dispersal = false)

#plotting results from the null model iteration
m1=prep_map(res_nm1[1],crop_to_ext=ex)
plot(m1)


#### running null model 2
res_nm2=null_models(example_species, 
                    "nm2",
                    geo_range, 
                    rs_std, 
                    copy(dom_master), 
                    top,
                    elv, 
                    formated_rs, 
                    nrep
                    )

#plotting results from the null model iteration
m2=prep_map(res_nm2[1],crop_to_ext=ex)
plot(m2)


#### running null model 3
res_nm3=null_models(example_species, 
                    "nm3",
                    geo_range, 
                    rs_std, 
                    copy(dom_master), 
                    top,
                    elv, 
                    formated_rs, 
                    nrep
                    )
#plotting results from the null model iteration
m3=prep_map(res_nm3[1],crop_to_ext=ex)
plot(m3)


#### running null model 4
@btime res_nm4=null_models(example_species, 
                    "nm4",
                    geo_range, 
                    rs_std, 
                    copy(dom_master), 
                    top,
                    elv, 
                    formated_rs, 
                    1
                    )
 #plotting results from the null model iteration
 
 
 m4=prep_map(res_nm4[1],crop_to_ext=ex)
 plot(first(res_nm4))
 
