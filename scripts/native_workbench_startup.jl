# Lightweight publication before Interior/factory compilation. This creates only
# a new directory; a token from this process is required to finish initialization.
module NativeWorkbenchStartup
using JSON3, SHA
function publish(root, value)
    path=joinpath(root,"worker.json");temp=path*".pending"
    open(io->JSON3.write(io,value),temp,"w");mv(temp,path;force=true)
end
function begin_startup(root, driver)
    root=abspath(root);ispath(root)&&error("fresh owned output directory required")
    mkpath(dirname(root));mkdir(root)
    started=time();token=bytes2hex(sha256(string(getpid(),time_ns(),root)))
    owner=(pid=getpid(),token=token,started_unix=started)
    open(io->JSON3.write(io,owner),joinpath(root,".startup-owner.json"),"w")
    publish(root,(schema="interior-native-workbench-1",status="starting",phase="Loading Julia and preparing paired native states",
        started_unix=started,driver_sha256=bytes2hex(sha256(read(driver)))))
    token
end
function validate(root,token)
    token isa String || error("fresh owned output directory required")
    islink(root)&&error("startup directory must be owned")
    owner=JSON3.read(read(joinpath(root,".startup-owner.json"),String))
    owner.pid==getpid()&&owner.token==token || error("startup ownership token mismatch")
    Set(readdir(root))<=Set((".startup-owner.json","worker.json")) || error("startup directory contains existing state")
    Float64(owner.started_unix)
end
function failed(root,token,e)
    owner=JSON3.read(read(joinpath(root,".startup-owner.json"),String))
    owner.pid==getpid()&&owner.token==token || error("startup ownership token mismatch")
    prior=JSON3.read(read(joinpath(root,"worker.json"),String),Dict{String,Any})
    prior["status"]="stopped";prior["stop_reason"]="initialization failed";prior["error"]=sprint(showerror,e)
    publish(root,prior)
end
end
