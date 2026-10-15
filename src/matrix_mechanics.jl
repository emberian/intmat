"""A paid catalyst matrix increases local drag without storing mechanical energy.

At the beginning of a constructive step, each cell samples total catalyst
concentration E over its physical exchange support. The held multiplier
ζ=1+resistance*E increases motor drag work and reduces passive mobility by 1/ζ.
The same sampled multiplier informs respiratory work forecasting. Newly secreted
matrix first affects mechanics on the next step; no ancestry or adhesion reward
is introduced. Resistance has units of physical volume per catalyst mass.
"""
struct MatrixDrag
    resistance::Float64
    function MatrixDrag(resistance::Real)
        value=Float64(resistance)
        isfinite(value) && value>=0 || throw(ArgumentError("matrix drag resistance must be finite nonnegative"))
        new(value)
    end
end
MatrixDrag(;resistance=20.)=MatrixDrag(resistance)

# Exceptional ratio ranges use bounded mantissas instead of overflowing or
# underflowing resistance/volume before multiplication by a finite amount.
function _matrix_scaled_amount(resistance,amount,volume)
    amount==0 && return 0.
    mr,er=frexp(resistance);ma,ea=frexp(amount);mv,ev=frexp(volume)
    ldexp((mr*ma)/mv,er+ea-ev)
end
_validate_matrix_mechanics(w,h,::Nothing)=nothing
function _validate_matrix_mechanics(w,h,law::MatrixDrag)
    law.resistance==0 && return nothing
    scale=law.resistance/h.cells.stencil.cell_volume
    for site in eachindex(h.enzymes[1])
        first=h.enzymes[1][site];second=h.enzymes[2][site]
        z=if isfinite(scale) && scale>=floatmin(Float64)
            1+scale*first+scale*second
        else
            1+_matrix_scaled_amount(law.resistance,first,h.cells.stencil.cell_volume)+
                _matrix_scaled_amount(law.resistance,second,h.cells.stencil.cell_volume)
        end
        isfinite(z) && z>=1 || throw(ArgumentError("matrix drag exceeds finite representable range"))
    end
    nothing
end
_matrix_drag_coefficients(::Nothing,environment)=nothing
function _matrix_drag_coefficients(law::MatrixDrag,environment)
    law.resistance==0 && return ones(size(environment,2))
    drag=[1+law.resistance*environment[3,i]+law.resistance*environment[4,i] for i in axes(environment,2)]
    all(z->isfinite(z)&&z>=1,drag) || throw(ArgumentError("sampled matrix drag must be finite and at least one"))
    drag
end
_matrix_drag_forecast(planned,w,c,i,delta,q,::Nothing)=planned
function _matrix_drag_forecast(planned,w,c,i,delta,q,drag)
    # The original forecast is retained; this adds only the extra work for
    # physically active axes (the 2D body's z action does not move it).
    nd=_cellular_ndim(w.grid)
    planned+delta*q.motor_cost*q.speed^2*(drag[i]-1)*sum(abs2,view(c.last_actions,1:nd,i))
end
