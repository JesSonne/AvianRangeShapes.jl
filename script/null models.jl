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


results = [null_models(example_species,
                       geo_range,
                       copy(dom_master),
                       top, 
                       ele_range,
                       formated_rs,
                       nrep,
                       rs_std;
                       bounded_dispersal=disp,
                       range_coherence=cohe,
                      ) 
               for disp in (false,true), cohe in (true,false)]
    
maps=[compile_map(results[i,j][1],dom_master) for i in (1,2),j in (1,2)]
maps=[crop_map(maps[i,j]; crop_to_ext=ex) for i in (1,2),j in (1,2)]
plot([plot(res, title="Null model $i", ticks=false) for (i, res) in enumerate(maps)]..., layout=(2,2))

