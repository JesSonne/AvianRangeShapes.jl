#add https://github.com/mkborregaard/SpreadingDye.jl
using ArchGDAL,Rasters, DataFrames,Plots, SpreadingDye, NearestNeighbors, Images, JLD2, SkipNan, VerySimpleRasters

include("functions.jl")

#loading geographical domain
dd = Raster("Data files/sf1_mainland.tif")  
dom = dd 
dom=dom./dom
dom=dom.==1
dom_master=copy(dom)
all_true=copy(dom_master)
all_true[:].=true
nas=findall(dom[:].==false) 


#empty raster object
zero=copy(dom_master);zero[:].=false
zero_na=copy(dd);zero_na[:].=NaN


#loading topographical raster
top=Raster("Data files/top_q_proj.tif")
top=resample(top;to=dom)


#loading species elevational range limits
elv=load_object("Data files/elevational range limits.jld2")

#loading the geographic geographic ranges of six example species
dis_elv=load_object("Data files/bird ranges.jld2")
nam=["Acropternis orthonyx","Aglaeactis castelnaudii","Coeligena lutetiae","Diglossa mystacalis","Phlogophilus harterti","Scytalopus griseicollis"]

#loading standardized range sizes
formated_rs=load_object("Data files/standardized_range_sizes.jld2")

############################ run the null model 

example_species="Phlogophilus harterti" # name of one of the example species from the nam object
rs_std=false # should the null model use the standardized range size (true) or the empirical range size (false)


#### empirical range
map_emp=copy(zero_na);map_emp[dis_elv[findall(nam.==example_species)][1]].=1
map_emp=Rasters.trim(map_emp,pad=10)
bd=bounds(map_emp)
plot(map_emp)


#### running null model 1
res_nm1=null_models(example_species, # state name of example species
                    "nm1",# state which of the four null models (nm1,nm2,nm3, or nm4)
                    rs_std
                    )

map_nm1=copy(zero_na);map_nm1[res_nm1].=1
map_nm1=Rasters.crop(map_nm1,to=map_emp)
plot(map_nm1)


#### running null model 2
res_nm2=null_models(example_species, 
                    "nm2",
                    rs_std
                    )

map_nm2=copy(zero_na);map_nm2[res_nm2].=1
map_nm2=Rasters.crop(map_nm2,to=map_emp)
plot(map_nm2)


#### running null model 3
res_nm3=null_models(example_species, 
                    "nm3",
                    rs_std
                    )

map_nm3=copy(zero_na);map_nm3[res_nm3].=1
map_nm3=Rasters.crop(map_nm3,to=map_emp)
plot(map_nm3)


#### running null model 4
res_nm4=null_models(example_species, 
                    "nm4",
                    rs_std
                    )

map_nm4=copy(zero_na);map_nm4[res_nm4].=1
map_nm4=Rasters.crop(map_nm4,to=map_emp)
plot(map_nm4)
