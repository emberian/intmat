# Offline, CPU-only selection audit. Reads the exact archived propensity law.
using Interior, SHA
function electrical_weight_law(path=joinpath(@__DIR__,"electrothermal_archive.jl"))
 source=read(path,String)
 a=findfirst("function choose_archive(",source);b=findfirst("function archive_geometry(",source)
 a!==nothing&&b!==nothing||error("archive propensity source markers changed")
 body=source[first(a):first(b)-1]
 cutoff=findfirst(" target=Interior.electrical_uniform",body)
 cutoff!==nothing||error("archive selection source marker changed")
 prefix=body[1:first(cutoff)-1]
 generated=replace(prefix,"function choose_archive("=>"function diagnostic_archive_weights(";count=1)*
  " (sites=sites,weights=weights,total=total)\nend\n"
 module_name=Symbol("ElectricalWeightLaw_",bytes2hex(sha256(source))[1:12])
 owner=Module(module_name);Core.eval(owner,:(using Interior))
 Base.include_string(owner,generated,path*":propensity-prefix")
 f=Base.invokelatest(getproperty,owner,:diagnostic_archive_weights)
 (weights=f,source_sha256=bytes2hex(sha256(source)),generated_source=generated,
  generated_sha256=bytes2hex(sha256(generated)))
end
function electrical_cdf_audit(law,p,sites,values,doses,foods,c)
 data=(values,zeros(Float32,length(sites)),doses,foods)
 distribution=Base.invokelatest(law.weights,p,sites,data,c)
 distribution===nothing&&error("no viable frontier")
 weights=distribution.weights;total=distribution.total
 u=Interior.electrical_uniform(p.seed,p.events+1);target=u*total
 acc=0.;cdf=Float64[]
 for weight in weights
  acc+=weight;push!(cdf,acc/total)
 end
 # Use the production's strict cumulative comparison, without sentinel ambiguity.
 acc=0.;chosen=length(sites);before=0.;after=total
 for (i,weight) in enumerate(weights)
  previous=acc;acc+=weight
  if acc>target;chosen=i;before=previous;after=acc;break;end
 end
 (site=Tuple(sites[chosen]),chosen=chosen,uniform=u,weights=weights,
  probabilities=weights./total,cdf=cdf,total=total,
  selected_probability=weights[chosen]/total,
  lower_margin=(target-before)/total,upper_margin=(after-target)/total,
  cdf_margin=min(target-before,after-target)/total)
end
function electrical_cdf_compare(a,b)
 length(a.weights)==length(b.weights)||throw(DimensionMismatch("common ordered frontier required"))
 cdf_difference=maximum(abs.(a.cdf.-b.cdf))
 (same_choice=a.site==b.site,cdf_maximum_difference=cdf_difference,
  total_variation=.5sum(abs.(a.probabilities.-b.probabilities)),
  a_margin=a.cdf_margin,b_margin=b.cdf_margin,
  a_choice_certified_against_b=cdf_difference<min(a.lower_margin,a.upper_margin))
end
