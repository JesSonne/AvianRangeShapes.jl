using SpreadingDye, NearestNeighbors, SkipNan, StatsBase, Rasters, ImageMorphology


# ——————————————————————————————————————————————
# 3) Edge‐on‐segment test
function point_on_edge(x, y, x1, y1, x2, y2; tol=1e-10)
    # collinearity?
    if abs((x2-x1)*(y-y1) - (y2-y1)*(x-x1)) > tol
        return false
    end
    # within bounding box?
    if x < min(x1,x2)-tol || x > max(x1,x2)+tol ||
       y < min(y1,y2)-tol || y > max(y1,y2)+tol
        return false
    end
    return true
end

function points_outside_hull(hull_pts::AbstractMatrix{<:Real},
                              query_pts::AbstractMatrix{<:Union{Real,Missing}};
                              tol::Real=1e-8)

    # build and close the polygon
    hidx   = convex_hull_indices(hull_pts)
    px     = hull_pts[hidx, 1]
    py     = hull_pts[hidx, 2]
    px_c   = [px; px[1]]
    py_c   = [py; py[1]]

    # for fast vertex‐lookup
    vertset = Set(zip(px, py))

    M = size(query_pts, 1)
    out = Vector{Union{Bool,Missing}}(undef, M)

    for k in 1:M
        row = query_pts[k, :]
        if any(ismissing, row)
            out[k] = missing
            continue
        end
        x, y = row[1], row[2]

        # 1) exact vertex?
        if (x,y) in vertset
            out[k] = false
            continue
        end

        # 2) on any edge?
        hit_edge = false
        for i in 1:length(px)
            if point_on_edge(x, y, px_c[i], py_c[i], px_c[i+1], py_c[i+1]; tol=tol)
                out[k] = false
                hit_edge = true
                break
            end
        end
        if hit_edge
            continue
        end

        # 3) ray‐casting
        inside = point_in_polygon(x, y, px, py)
        out[k] = !inside
    end

    return out
end


   #filtering the geographic domain by the species elevational range limits     
function update_dom_to_elevation!(dom, species, ele_range, top) 
    if !(species in keys(ele_range))
        warn("No elevational range data found for $(species)")
    else
        dom[top[Band=1] .> ele_range[species][2] .|| top[Band=2] .< ele_range[species][1]] .= false
    end
end

function update_dom_to_climate_volume!(dom, ab, clim_mat) 
    if length(ab)>10
        cl=clim_mat[ab,1:2]
        out = points_outside_hull(cl, clim_mat[:,1:2])
        dom[findall(out)].=false
    
        cl=clim_mat[ab,3:4]
        out = points_outside_hull(cl, clim_mat[:,3:4])
        dom[findall(out)].=false
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
    point = [Tuple.(CartesianIndices(dom)[group]) for group in groups]
    
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

#the outward_facing function to run each null model
function null_models(
    species::String, # state name of example species
    geo_range::Dict, #grid cell ids comprising the species empirical range
    dom_master::Any, #biogeographical domain
    top::Any, #topographical raster
    ele_range::Any, #data frame with the species' elevational range limits
    clim_mat::Any #raster stack with climate variables   
    formated_rs::Any, #data frame with standardized range sizes (only used if rs_std=true)
    nrep::Int64, # nuber of repetitions
    rs_std::Bool; # should the null model use the standardized range size (true) or the empirical range size (false)
    bounded_dispersal::Bool, #true for bounded dispersal, false for free dispersal
    range_coherence::Bool, #true for range coherence, false for patchy ranges 
    )

    dom = copy(dom_master)
    #filtering the geographic domain by the species elevational range limits     
    if bounded_dispersal
        update_dom_to_elevation!(dom, species, ele_range, top)
    end

    #list of grid cells outside the geographic domain
    nas=findall(dom[:].==false) 

    #empty raster object
    zero=copy(dom);zero[:].=false

    #grid cell ids comprising the species' empirical range
    ab=geo_range[species]

    if bounded_dispersal
        update_dom_to_climate_volume!(dom,ab,clim_mat)
    end
    #constructing raster of the species empirical range
    emp=Float64.(dom)
    emp.=NaN
    emp[ab].=true
    emp2= (!isnan).(emp)  

    if range_coherence
        group_size=[length(ab)]
        groups=[ab]
    end

    if !range_coherence
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

function crop_map(map;trim_map=true,crop_to_ext=nothing)
    if crop_to_ext === nothing
                if trim_map
                    map_nm=trim_raster(map,10)
                end
            end
            if !(crop_to_ext === nothing)
                map_nm = Rasters.crop(map; to=crop_to_ext)
    end
   map_nm 
end


function compile_map(res_nm,dom=dd)
    map_nm=copy(dom)
    map_nm[:].=0
    map_nm[res_nm].=1
    map_nm
end

function trim_raster(r::AbstractMatrix{Bool}, pad::Int=0)
    # Find row indices with at least one true value.
    row_indices = findall(row -> any(row), eachrow(r))
    # Find column indices with at least one true value.
    col_indices = findall(col -> any(col), eachcol(r))
    
    # If no true values exist, return an empty array.
    if isempty(row_indices) || isempty(col_indices)
        return similar(r, 0, 0)
    end
    
    # Determine the raw bounds.
    row_min, row_max = minimum(row_indices), maximum(row_indices)
    col_min, col_max = minimum(col_indices), maximum(col_indices)
    
    # Expand the bounds by the specified padding, ensuring they remain within the raster.
    row_min_p = max(1, row_min - pad)
    row_max_p = min(size(r, 1), row_max + pad)
    col_min_p = max(1, col_min - pad)
    col_max_p = min(size(r, 2), col_max + pad)
    
    return r[row_min_p:row_max_p, col_min_p:col_max_p]
end
