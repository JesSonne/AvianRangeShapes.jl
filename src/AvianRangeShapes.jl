module AvianRangeShapes

include("functions.jl")

function null_models(example_species="Phlogophilus harterti", # state name of example species
    Anal_nam="nm4",# state which of the four null models (nm1,nm2,nm3, or nm4)
    rs_std=false # should the null model use the standardized range size (true) or the empirical range size (false)
    )

    i=findall(nam.==example_species)
    dom=copy(dom_master)



    #filtering the geographic domain by the species elevational range limits     
    if Anal_nam in ["nm2","nm4"]
        if example_species in elv.Species
            sp_elv_sub=elv[elv.Species.==example_species,:]
            elv_rm=[findall(top[Band=1][:].>sp_elv_sub.Maximu_elevation);findall(top[Band=2][:].<sp_elv_sub.minimum_elevation)]
            dom[elv_rm].=false
        end   
    end

    #list of grid cells outside the geographic domain
    nas=findall(dom[:].==false) 


    #grid cells comprising the species empirical range
    ab=dis_elv[i][1]


    #constructing raster of the species empirical range
    emp=copy(dd)
    emp.=NaN
    emp[ab].=true
    emp2= (!isnan).(emp)  

    if Anal_nam in ["nm1","nm2"]
        group_size=[length(ab)]
        groups=[ab]
    end

    if Anal_nam in ["nm3","nm4"]
        #Identify groups of isolated range patches 
        groups = collect_groups(emp2)
        groups=groups[findall(length.(groups).>0)]

        #### algorithm treating tiny range patches as extensions of the agacent larger coherent range rather than independent range patches
        groups=join_neighbours(groups; 
        max_dist=5, #minimum patch size relative to the species’ largest patch 
        min_prop=0.1) #minimum distance between range patches


        #extracting the size of each patch and the entire range
        group_size=length.(groups)
    end    

    total_rangesize=sum(group_size)

    ######### ## standardizing the range size frequency distribution
    if rs_std
        #if the standardized range size is smaller than the empirical, randomly subtract grid cells from the groups in stepwise fashion, weighted by the patch size
        #i.e. large patches have greater chance of beeing modified by the standardization
        new_range=formated_rs[formated_rs.nam.==example_species,:rank_range][1]
        if new_range<total_rangesize 
            for x in 1:(total_rangesize-new_range)
                subt=sample(collect(1:length(group_size)),Weights(group_size)) 
                group_size[subt]=group_size[subt]-1  
            end
        end 

        #if the standardized range size is larger than the empirical, randomly add grid cells from the groups in stepwise fashion, weighted by the patch size
        #i.e. large patches have greater chance of beeing modified by the standardization
        if new_range>total_rangesize 
            for x in 1:(new_range-total_rangesize)
                subt=sample(collect(1:length(group_size)),Weights(group_size)) 
                group_size[subt]=group_size[subt]+1  
            end
        end 

        total_rangesize =new_range
        ppo=findall(group_size.>0)
        group_size=group_size[ppo]
        groups=groups[ppo]

    end



    #running the spreading dye algorithm for each range patch
    sd_out=copy(zero)
    for x in 1:length(groups)
        sd_sub=copy(zero)
        sd_sub[groups[x]].=true
        ppo=Tuple.(collect(CartesianIndices(sd_sub))[sd_sub[:]])
        rangesize = group_size[x]
        start=(rand(ppo,1))[1] # select random starting posision
        sd = spreading_dye(rangesize, dom, start) 
        id=findall(sd[:].==true)
        sd_out[id].=true
    end


    #checking for overlaps in the simulated range patches (i.e. if simulated range size is less than the empirical)
    ids=findall(sd_out[:].==true)    
    rs_check=total_rangesize-length(ids)

    if rs_check > 0
        sd_out=expand_spreading(sd_out,rs_check,dom)
    end    
    ids=findall(sd_out[:].==true)     
    ids

end



export null_models
end
