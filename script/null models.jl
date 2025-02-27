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
map_emp=crop_map(map_emp; trim_map=true)
ex=ArchGDAL.extent(map_emp)
plot(map_emp)

results = [null_models(example_species,           # state name of example species
                       geo_range,                 # grid cell ids comprising the species empirical range
                       copy(dom_master),          # biogeographical domain
                       top,                       # topographical raster
                       ele_range,                 # data frame with the species' elevational range limits
                       formated_rs,               # data frame with standardized range sizes (only used if rs_std=true)
                       nrep,                      # number of repetitions
                       rs_std;                    # should the null model use the standardized range size (true) or the empirical range size (false)
                       bounded_dispersal=disp,    # true for bounded dispersal, false for free dispersal
                       range_coherence=cohe,      # true for range coherence, false for patchy ranges 
                      ) 
               for disp in (false,true), cohe in (true,false)]
    
maps=[compile_map(results[i,j][1],dom_master) for i in (1,2),j in (1,2)]
maps=[crop_map(maps[i,j]; crop_to_ext=ex) for i in (1,2),j in (1,2)]
plot([plot(res, title="Null model $i", ticks=false) for (i, res) in enumerate(maps)]..., layout=(2,2))

