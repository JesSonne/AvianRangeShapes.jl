using SpreadingDye, NearestNeighbors, SkipNan, StatsBase, Rasters, ImageMorphology


   #filtering the geographic domain by the species elevational range limits     
function update_dom_to_elevation!(dom, species, elv, top) 
    if !(species in elv.Species)
        warn("No elevational range data found for $(species)")
    else
        sp_elv_sub=elv[elv.Species.==species,:]
        dom[top[Band=1] .> sp_elv_sub.Maximu_elevation .|| top[Band=2] .< sp_elv_sub.minimum_elevation] .= false
    end
end

#Identify groups of isolated range patches 
function find_groups(emp2,dom,
    max_dist=5, #minimum patch size relative to the species’ largest patch 
    min_prop=0.1) #minimum distance between range patches) 
   
    #Identify groups of isolated range patches 
   groups = collect_groups(emp2)
   groups=groups[findall(length.(groups).>0)]

   #### algorithm treating tiny range patches as extensions of the agacent larger coherent range rather than independent range patches
   groups=join_neighbours(groups,dom; 
   max_dist, #minimum patch size relative to the species’ largest patch 
   min_prop) #minimum distance between range patches

   groups
end

#stadardizing range sizes of groups
function update_group_size!(new_range,total_rangesize,group_size)        
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
end

#compiling arguments and running the spreading die model nrep times
function Run_SpreadingDye(groups,group_size,dom,nrep,zero)
    total_rangesize=sum(group_size)
    output=Any[]
    for i in 1:nrep
        #running the spreading dye algorithm for each range patch
        sd_out=copy(zero)
        for x in 1:length(groups)
            sd_sub=copy(zero)
            sd_sub[groups[x]].=true
            ppo=Tuple.(collect(CartesianIndices(sd_sub))[sd_sub[:]])
            rangesize = group_size[x]
            start=(rand(ppo,1))[1] # select random starting posision
            sd = SpreadingDye.spreading_dye(rangesize, dom, start) 
            id=findall(sd[:])
            sd_out[id].=true
        end

        #checking for overlaps in the simulated range patches (i.e. if simulated range size is less than the empirical) 
        ids=findall(sd_out[:])
        rs_check=total_rangesize-length(ids)

        if rs_check > 0
            sd_out=expand_spreading(sd_out,rs_check,dom)
        end    
        ids=findall(sd_out[:])
        push!(output,ids)     
    end

    output
end
    


#the outward_facing function to run each null model
function null_models(
    species::String, # state name of example species
    Anal_nam::String,# state which of the four null models (nm1,nm2,nm3, or nm4)
    geo_range::Dict, #grid cell ids comprising the species empirical range
    rs_std::Bool, # should the null model use the standardized range size (true) or the empirical range size (false)
    dom_master::Any, #biogeographical domain
    top::Any, #topographical raster
    elv::Any, #data frame with the species' elevational range limits
    formated_rs::Any, #data frame with standardized range sizes (only used if rs_std=true)
    nrep::Int64 # nuber of repetitions
    )

    dom = copy(dom_master)
    #filtering the geographic domain by the species elevational range limits     
    if Anal_nam in ["nm2","nm4"]
        update_dom_to_elevation!(dom, species, elv, top)
    end

    #list of grid cells outside the geographic domain
    nas=findall(dom[:].==false) 

    #empty raster object
    zero=copy(dom);zero[:].=false

    #grid cell ids comprising the species' empirical range
    ab=geo_range[species]


    #constructing raster of the species empirical range
    emp=Float64.(dom)
    emp.=NaN
    emp[ab].=true
    emp2= (!isnan).(emp)  

    if Anal_nam in ["nm1","nm2"]
        group_size=[length(ab)]
        groups=[ab]
    end

    if Anal_nam in ["nm3","nm4"]
        groups=find_groups(emp2,dom)
        #extracting the size of each patch and the entire range
        group_size=length.(groups)
    end    

    total_rangesize=sum(group_size)

    ######### ## standardizing the range size frequency distribution
    if rs_std
        #if the standardized range size is smaller than the empirical, randomly subtract grid cells from the groups in stepwise fashion, weighted by the patch size
        #i.e. large patches have greater chance of beeing modified by the standardization
        new_range=formated_rs[formated_rs.nam.==species,:rank_range][1]
        
        update_group_size!(new_range,total_rangesize,group_size)    
        println(group_size)

        total_rangesize =new_range
        ppo=findall(group_size.>0)
        
        group_size=group_size[ppo]

        groups=groups[ppo]

    end

   Run_SpreadingDye(groups,group_size,dom,nrep,zero)
end



#Checking if selected grid cell is on the geographical domain
on_domain(dom::AbstractMatrix{Bool}, point::Tuple{Int, Int}) = within_edges(dom, point) && dom[point...]
within_edges(dom::AbstractMatrix, point::Tuple{Int, Int}) = min(point...) > 0 && first(point) <= size(dom, 1) && last(point) <= size(dom, 2)

#compute the number of isolated range patches along with the patches' grid cell ids    
function collect_groups(emp2)
    labels = label_components(emp2,strel_box((3, 3))) # queen style neighbourhood
    groups = [Int[] for i = 1:maximum(labels)]
    for (i,l) in enumerate(labels)
        if l != 0
            push!(groups[l], i)
        end
    end
    groups
end

## stating neigbborhood type
const algos = Dict(
    :rook => ((-1,0), (0, 1), (1, 0), (0, -1)),
    :queen => Tuple((x,y) for x in -1:1, y in -1:1 if !(x == y == 0))    
)
algo=:rook


### constructing neigbborhood matrix for the range patches
function reldists(point)
    a = fill(Inf, length(point), length(point))
    sort!(point, by = length, rev = true)
    for i in 1:length(point)-1
        m1 = Float64.(stack(point[i]))
        kdtree = KDTree(m1; leafsize = 8)
        for j in i+1:length(point)
            m2 = Float64.(stack(point[j]))
            ind, dist = knn(kdtree, m2, 1)
            a[i, j] = only(minimum(dist)); a[j, i] = a[i,j]
        end
    end
    a
end


#parameters: maximum distance between patches and minimum percentage size of patches
function join_neighbours(groups,dom;max_dist::Int64=5,min_prop::Float64=0.1)
    #organizing patch groups
    groups=sort(groups, by = length, rev = true)
    group_size=length.(groups)
    group_size_prop=group_size./maximum(group_size)

    ## convert raster ID's to matrix coordinates
    point=Any[]
    for t in 1:length(groups)
        zz=copy(Float64.(dom));zz.=NaN;zz[groups[t]].=1;push!(point,Tuple.(collect(CartesianIndices(zz))[isfinite.(zz)]))
    end

    ###check if the following is faster than the above: point = [Tuple.(CartesianIndices[group]) for group in groups]

    #construct neigbborhood matrix
    dm=reldists(point)
    p = Set.(1:length(groups))
    ll = group_size
    id=copy(dm);id[:].=1:length(id[:])

    # identify small patches
    small=reverse(findall(group_size_prop.<min_prop)) # reverse to start merging from the smallest patch to the largest


    while length(small)>0
        po=small
        mer=Tuple.(argmin(dm[po,:])) # picks first element (will prioritize merging with the largest patch)
        mer2=Tuple.(findall(id.==id[po,:][mer...])[1])
        mi=minimum(mer2);ma=maximum(mer2)

        if dm[mer2...]<max_dist
            p[mi] = union(p[mi], p[ma])
            empty!(p[ma])
            ll[mi]=ll[ma]+ll[mi]; ll[ma]=0
            
            comp=[reshape(dm[mi,:], 1, :); reshape(dm[ma,:], 1, :)]
            mins=minimum.(eachcol(comp))
            infs=.!isfinite.(dm)
            dm[ma,:]=mins
            dm[:,mi]=mins
            dm[infs].=Inf
        else ll[ma]=0 end

        dm[mi,ma]=Inf
        dm[ma,mi]=Inf
        group_size_prop=ll./maximum(ll)
        small= reverse(findall((group_size_prop.<min_prop) .& (group_size_prop.>0)))# reverse to start merging from the smallest patch to the largest
    end

    out=Any[]
    ss=findall(length.(p).>0)
    for a in 1:length(ss)
        push!(out,reduce(vcat, groups[collect(p[ss[a]])]))
    end

    out
end

function find_edges(georange::AbstractMatrix{Bool}, dom::AbstractMatrix{Bool}, algo::Symbol)
    Set(((i,j), (i,j).+nb) 
    for i in axes(dom, 1), j in axes(dom, 2), nb in algos[algo] 
        if within_edges(dom, (i,j).+nb) && georange[i,j] && !georange[((i,j).+nb)...]
    )
end

function expand_spreading(georange::AbstractMatrix{Bool}, add_cells::Int, dom::AbstractMatrix{Bool}; algo::Symbol = :rook)
    edges = find_edges(georange, dom, algo)
    for i in 1:add_cells
        grow!(georange, edges, dom, algo)
    end
    georange
end

function grow!(georange::AbstractMatrix{Bool}, edges::Set, dom::AbstractMatrix{Bool}, algo::Symbol; ignore_domain = false)
    newcell = isempty(edges) ? pick_random(georange) : first(rand(edges))
    while georange[newcell...]
        isempty(edges) && push!(edges, jump(georange, dom, algo)) # allows for patchy ranges
        edge, newcell = rand(edges)
        pop!(edges, (edge, newcell))
    end
    georange[newcell...] = true
    for nb in algos[algo]
        neighbor = newcell .+ nb
        if within_edges(dom, neighbor) && (ignore_domain || dom[neighbor...]) && !georange[neighbor...]
            push!(edges, (newcell, neighbor))
        end
    end
    newcell
end

function prep_map(res_nm,dom=dd;trim_map=true,crop_to_ext=nothing)
    map_nm=copy(dom)
    map_nm[:].=NaN
    map_nm[res_nm].=1
    if crop_to_ext === nothing
        if trim_map
            map_nm=Rasters.trim(map_nm,pad=10)
        end
    end
    if !(crop_to_ext === nothing)
        map_nm = Rasters.crop(map_nm, to=crop_to_ext)
    end
    plot(map_nm)
    map_nm
end

