# Atomic weight data
const A = Float64[
  1.008,
  4.0026022,
  6.94,
  9.01218315,
  10.81,
  12.011,
  14.007,
  15.999,
  18.9984031636,
  20.17976,
  22.989769282,
  24.305,
  26.98153857,
  28.085,
  30.9737619985,
  32.06,
  35.45,
  39.9481,
  39.09831,
  40.0784,
  44.9559085,
  47.8671,
  50.94151,
  51.99616,
  54.9380443,
  55.8452,
  58.9331944,
  58.69344,
  63.5463,
  65.382,
  69.7231,
  72.6308,
  74.9215956,
  78.9718,
  79.904,
  83.7982,
  85.46783,
  87.621,
  88.905842,
  91.2242,
  92.906372,
  95.951,
  98.0,
  101.072,
  102.905502,
  106.421,
  107.86822,
  112.4144,
  114.8181,
  118.7107,
  121.7601,
  127.603,
  126.904473,
  131.2936,
  132.905451966,
  137.3277,
  138.905477,
  140.1161,
  140.907662,
  144.2423,
  145.0,
  150.362,
  151.9641,
  157.253,
  158.925352,
  162.5001,
  164.930332,
  167.2593,
  168.934222,
  173.0451,
  174.96681,
  178.492,
  180.947882,
  183.841,
  186.2071,
  190.233,
  192.2173,
  195.0849,
  196.9665695,
  200.5923,
  204.38,
  207.21,
  208.980401,
  209.0,
  210.0,
  222.0,
  223.0,
  226.0,
  227.0,
  232.03774,
  231.035882,
  238.028913,
  237.0,
  244.0,
  243.0,
  247.0,
  247.0,
  251.0,
  252.0,
]

function readdata(z::Int)
  fn = joinpath(@__DIR__, "..", "data", "phxs[$z].dat")
  (ee, df) = open(fn) do f
    ee = Dict{Tuple{Int,Int},Float64}()
    energy, muorho = Float64[], Float64[]
    nsh, npts = -1, -1
    for (i, l) in enumerate(eachline(f))
      if startswith("# ", l[1])
        if i == 3
          nsh = parse(Int, strip(l[2:6]))
        elseif i >= 7 && i <= 6 + nsh
          sh = parse(Int, strip(l[3:4]))
          e = parse(Float64, l[19:31])
          ee[(z, sh)] = e
        elseif nsh != -1 && i == nsh + 9
          npts = parse(Int, strip(l[3:8]))
        end
      end
      if nsh != -1 && i >= nsh + 13 && i <= nsh + 12 + npts
        s = parse(Float64, l[3:15])
        x = parse(Float64, l[18:30])
        push!(energy, s)
        push!(muorho, 6.02214076E+023 * x / A[z])
      elseif i > nsh + 12 + npts
        break
      end
    end
    (ee, (energy, muorho))
  end
end

const MACData = Vector{Tuple{Vector{Float64}, Vector{Float64}}}()
const EdgeData = Dict{Tuple{Int,Int}, Float64}()

function __init__()
  for z = 1:99
    ee, mac = readdata(z)
    merge!(EdgeData, ee)
    push!(MACData, mac)
  end
end

"""
    eachelement()

The range of available elements.
"""
eachelement() = eachindex(MACData)

"""
    eachedge(z::Int)::Set{Integer}

Returns a set containing the shells for which there is an edge energy in the database
for the specified element.
"""
eachedge(z::Int) = Set{Integer}(filter(sh->haskey(EdgeData, (z, sh)), 1:29))

"""
    hasedge(::Type{SabbatucciMAC}, z::Int, shell::Int)::Bool

Is a value available for the specific shell's edge energy for the element identified by atomic number, z.
"""
hasedge(z::Int, shell::Int) = haskey(EdgeData, (z, shell))

"""
    atomicweight(::Type{SabbatucciMAC}, z::Int)::Float64

The mean atomic weight for the specified element.
"""
atomicweight(z::Int) = A[z]


"""
    binarysearch(lst::Vector{T}, val::T)::Int

If `val` exists in the sorted Vector `lst`, returns the index of `val` in `lst`.  If `val` does not exist in `lst`, returns `i` where
`lst[-i]<val` and `lst[-i+1]>val`.
"""
function binarysearch(lst::Vector{T}, val::T)::Int where {T}
  low, high = 1, length(lst)
  while low ≤ high
    mid = (low + high) ÷ 2
    if lst[mid] > val
      high = mid - 1
    elseif lst[mid] < val
      low = mid + 1
    else
      return mid
    end
  end
  return -low + 1
end

"""
    edgeenergy(::Type{SabbatucciMAC}, z::Int, shell::Int)::Float64

The source of the edge energy data is Williams, G.P., 2011. Electron binding energies of the elements.
In: Haynes, W.M., Lide, D.R. (Eds.), CRC Handbook of Chemistry and Physics, 91st edition CRC, Boca Raton, pp. 221–226.

Throws an error if (z,sh) does not exist.
"""
edgeenergy(z::Int, sh::Int)::Float64 = EdgeData[(z, sh)]

function mac(z::Int, energy::Float64)::Float64
  energies, macs = MACData[z]
  i = binarysearch(energies, energy)
  if i < 1
    if i == 0 || i == -length(energies)
      return NaN
    else
      return macs[-i] + (energy - energies[-i]) * ((macs[-i+1] - macs[-i]) / (energies[-i+1] - energies[-i]))
    end
  else
    return macs[i]
  end
end
