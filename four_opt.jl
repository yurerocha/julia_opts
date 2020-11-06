# 4-Opt neighborhood move in O(n²) time.

using Random

mutable struct Coord
    x::Float64
    y::Float64
end

mutable struct DistMatrix
    customer::Dict{Int, Coord} # Customer id and its coordinates.
    cost::Dict{Tuple{Int, Int}, Float64} # Cost between two customers.
end

mutable struct Solution
    route::Vector{Int}
    cost::Float64
end

mutable struct FourOptMoveData
    x1::Char
    x2::Char
    i1::Int
    j1::Int
    i2::Int
    j2::Int
    cost_var::Float64
end

function gen_dist_matrix(n::Int, max_dist::Float64)
    v = Dict{Int, Coord}()

    # Generate random customers' coordinates and the depot coordinates.
    for i in 1:(n + 1)
        v[i] = Coord(floor(rand(Float64)*max_dist), floor(rand(Float64)*max_dist))
    end

    c = Dict{Tuple{Int, Int}, Float64}()
    # Compute the Euclidean distance between all customers.
    for i in 0:n, j in 0:n
        c[i, j] = floor(sqrt((v[j + 1].x - v[i + 1].x)^2 + (v[j + 1].y - v[i + 1].y)^2))
    end

    return DistMatrix(v, c)
end

function gen_init_sol(D::DistMatrix)
    v = deepcopy(D.customer)

    sol = Solution(Vector{Int}(), 0.0)
    push!(sol.route, 0)
    pop!(v, 1)
    # Create the initial route solution randomly.
    while !isempty(v)
        el = rand(v)
        push!(sol.route, el.first - 1)
        pop!(v, el.first)
    end
    push!(sol.route, 0)
    sol.cost = compute_cost(D, sol)

    return sol
end

function compute_cost(D::DistMatrix, sol::Solution)
    c = 0.0

    for i in 2:length(sol.route)
        c = c + D.cost[sol.route[i - 1], sol.route[i]]
    end

    return c
end

# 4-Opt neighborhood move
function deltax(D::DistMatrix, sol::Solution, x::Char, i::Int, j::Int)
    r = sol.route
    cost_add = 0.0
    if x == 'c'
        cost_add = D.cost[r[i], r[j]] + D.cost[r[i + 1], r[j + 1]] # Cost of a connecting 2-opt.
    else
        cost_add = D.cost[r[i], r[j + 1]] + D.cost[r[i + 1], r[j]] # Cost of a disconnecting 2-opt.
    end

    return cost_add - D.cost[r[i], r[i + 1]] - D.cost[r[j], r[j + 1]]
end

function perform_four_opt!(s::Solution, d::FourOptMoveData)
    pi1 = s.route[1: d.i1]
    pi2 = s.route[d.i1 + 1: d.i2]
    pi3 = s.route[d.i2 + 1: d.j1]
    pi4 = s.route[d.j1 + 1: d.j2]
    pi5 = s.route[d.j2 + 1: length(s.route)]

    nr = Vector{Int}()
    if (d.x1, d.x2) == ('d', 'd')
        nr = cat(pi1, pi4, pi3, pi2, pi5, dims=1)
    elseif (d.x1, d.x2) == ('c', 'd')
        nr = cat(pi1, reverse(pi3), reverse(pi4), pi2, pi5, dims=1)
    elseif (d.x1, d.x2) == ('d', 'c')
        nr = cat(pi1, pi4, reverse(pi2), reverse(pi3), pi5, dims=1)
    end

    s.route = nr
    s.cost = s.cost + d.cost_var
end

function four_opt!(D::DistMatrix, s::Solution)
    # Char is either 'c' cor connecting or 'd' for disconnecting;
    # Int, Int are the i and j indexes;
    #  Float64 is the cost.
    phi_sub = Dict{Tuple{Char, Int, Int}, Tuple{Int, Int, Float64}}()
    phi = Dict{Tuple{Char, Int, Int}, Tuple{Int, Int, Float64}}()
    # move 1 type, move 2 type, i2, j2, i1, j1, cost.
    delta_best = FourOptMoveData('c', 'c', -1, -1, -1, -1, 0.0)
    n = length(s.route)
    for j1 in 3:(n - 1), i2 in 2:(j1 - 1), x in ('c', 'd')
        if i2 == 2
            phi_sub[x, i2, j1] = 1, j1, deltax(D, s, x, 1, j1)
        else
            if phi_sub[x, i2 - 1, j1][3] < deltax(D, s, x, i2 - 1, j1)
                phi_sub[x, i2, j1] = phi_sub[x, i2 - 1, j1]
            else
                phi_sub[x, i2, j1] = i2 - 1, j1, deltax(D, s, x, i2 - 1, j1)
            end
        end
    end
    for i2 in 2:(n - 2), j2 in (i2 + 2):(n - 1)
        for x in ('c', 'd')
            if j2 == i2 + 2
                phi[x, i2, j2] = phi_sub[x, i2, i2 + 1]
            else
                if phi[x, i2, j2 - 1][3] < phi_sub[x, i2, j2 - 1][3]
                    phi[x, i2, j2] = phi[x, i2, j2 - 1]
                else
                    phi[x, i2, j2] = phi_sub[x, i2, j2 - 1]
                end
            end
        end
        # 4-opt complete move costs.
        delta_type1 = deltax(D, s, 'd', i2, j2) + phi['d', i2, j2][3] # Two disconnecting 2-opts.
        delta_type2A = deltax(D, s, 'd', i2, j2) + phi['c', i2, j2][3] # A connecting and a disconnecting 2-opt;
        delta_type2B = deltax(D, s, 'c', i2, j2) + phi['d', i2, j2][3] # A connecting and a disconnecting 2-opt.

        # Check and record the overall best (feasible?) improving move.
        best_cost_var = min(delta_best.cost_var, delta_type1, delta_type2A, delta_type2B)
        if best_cost_var == delta_type1
            delta_best = FourOptMoveData('d', 'd', phi['d', i2, j2][1], phi['d', i2, j2][2], i2, j2, delta_type1)
        elseif best_cost_var == delta_type2A
            delta_best = FourOptMoveData('c', 'd', phi['c', i2, j2][1], phi['c', i2, j2][2], i2, j2, delta_type2A)
        elseif best_cost_var == delta_type2B
            delta_best = FourOptMoveData('d', 'c', phi['d', i2, j2][1], phi['d', i2, j2][2], i2, j2, delta_type2B)
        end
    end

    if delta_best.cost_var < 0.0
        @show delta_best
        perform_four_opt!(s, delta_best)
        return true
    end

    return false
end

function run(ncusts::Int, max_dist::Float64)
    D = gen_dist_matrix(ncusts, max_dist)
    s = gen_init_sol(D)
    @show s
    println("")

    while four_opt!(D, s)
        @show s
        @assert s.cost == compute_cost(D, s)
        println("------------------------------")
    end
end
