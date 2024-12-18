
using ArchGDAL,Rasters, DataFrames,Plots, SpreadingDye, NearestNeighbors, Images, JLD2, SkipNan, VerySimpleRasters
using AvianRangeShapes


cd("/Users/jespersonne/Documents/GitHub/Avian_Range_Shapes")

#loading geographical domain
dd = Raster("data/sf1_mainland.tif")  
dom = dd 
dom=dom./dom
dom=dom.==1
dom_master=copy(dom)
all_true=copy(dom_master)
all_true[:].=true
nas=findall(dom[:].==false) 

#empty raster object
zero_na=copy(dd);zero_na[:].=NaN


#loading topographical raster
top=Raster("data/top_q_proj.tif")
top=resample(top;to=dom)


#loading species elevational range limits
elv=load_object("data/elevational range limits.jld2")

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
map_emp=copy(zero_na);map_emp[dis_elv[findall(nam.==example_species)][1]].=1
map_emp=Rasters.trim(map_emp,pad=10)
bd=bounds(map_emp)
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
                    )

                    

#plotting results from the null model iteration
map_nm1=copy(zero_na);map_nm1[res_nm1[1]].=1 
map_nm1=Rasters.crop(map_nm1,to=map_emp)
plot(map_nm1)


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
map_nm2=copy(zero_na);map_nm2[res_nm2[1]].=1
map_nm2=Rasters.crop(map_nm2,to=map_emp)
plot(map_nm2)


#### running null model 3
res_nm3=null_models(example_species, 
                    "nm3",
                    geo_range, 
                    rs_std, 
                    copy(dom_master), 
                    elv, 
                    formated_rs, 
                    nrep
                    )
#plotting results from the null model iteration
map_nm3=copy(zero_na);map_nm3[res_nm3[1]].=1
map_nm3=Rasters.crop(map_nm3,to=map_emp)
plot(map_nm3)


#### running null model 4
res_nm4=null_models(example_species, 
                    "nm4",
                    geo_range, 
                    rs_std, 
                    copy(dom_master), 
                    elv, 
                    formated_rs, 
                    nrep
                    )
                    #plotting results from the null model iteration
map_nm4=copy(zero_na);map_nm4[res_nm4[1]].=1
map_nm4=Rasters.crop(map_nm4,to=map_emp)
plot(map_nm4)
