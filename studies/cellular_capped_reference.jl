# Generate a private reference from the exact loaded production source. This
# changes no production method and avoids maintaining a second copy of its law.
function install_cellular_capped_reference!()
    path=joinpath(dirname(pathof(@__MODULE__)),"cellular_contacts.jl")
    source=read(path,String);digest=bytes2hex(SHA.sha256(source))
    arena_sources()[basename(path)]==digest||throw(ArgumentError("contacts source differs from loaded module; reload before generating reference"))
    start=findfirst("function cellular_mechanics!(",source)
    start===nothing&&error("mechanics declaration not found")
    finish=findnext("\n\"\"\"Conservative, positive exchange",source,last(start))
    finish===nothing&&error("mechanics end marker not found")
    generated=source[first(start):prevind(source,first(finish))]
    replacements=("function cellular_mechanics!("=>"function _cellular_capped_reference!(",
        "reuse_contact_workspace::Bool=true)"=>"reuse_contact_workspace::Bool=true,max_step=Inf)",
        "h=min(remaining,speed>"=>"h=min(remaining,max_step,speed>")
    for (before,after) in replacements
        count(before,generated)==1||error("production reference structure changed; inspect before generating")
        generated=replace(generated,before=>after)
    end
    Base.include_string(@__MODULE__,generated,"study-capped-cellular-reference.jl")
    (;method=Base.invokelatest(getproperty,@__MODULE__,:_cellular_capped_reference!),source=generated,contacts_sha256=digest)
end
