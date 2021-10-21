##
# Generate getter functions for struct fields
macro records(T)
    ret = Expr(:tuple)
    for name in fieldnames(Core.eval(__module__, T))
        Expr(:(function $(esc(name))(obj :: $T) getproperty(obj, Symbol($(esc(name)))) end))
    end
    return ret
end

##
struct Cuckold
    who
    whom
end

##
@records Cuckold

##
emily = Cuckold("Cal", "David")

##
who(emily)

##
whom(emily)

## Generated functions
@generated function Base.+:(a :: T, b :: T) where {T}
    if isempty(fieldnames(T))
        Expr(:call, T, :(a + b))
    else
        Expr(:call, T, map(field -> :(a.$field + b.$field), fieldnames(t))...)
    end
end