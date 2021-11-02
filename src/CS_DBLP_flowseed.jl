include("CS_DBLP.jl")
include("CP_FlowSeed.jl")

function GetFSSet(R)
    return LocalCond(B, R)[1]
end
