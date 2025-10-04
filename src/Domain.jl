macro domain(ex, msg::String, val)
    esc(:($ex ? nothing : throw(DomainError($val, $msg))))
end

macro domain(ex, msg::String, vals...)
    esc(:(@domain($ex, $msg, ($(vals...),))))
end

macro domain(ex, val)
    esc(:(@domain($ex, $(string("\nDomain:\n", ex, "\n")), $val)))
end

macro domain(ex, vals...)
    esc(:(
        @domain($ex, join($(string.(vals)) .* " = " .* string.(($(vals...),)), ", "))
    ))
end

macro domain_check(ex, m, n)
    quote
        @domain($ex,
            $("""
            \nDomain:
            $n ≥ 0
            |$m| ≤ $n
            $n ≡ $m (mod 2)
            """),
            $m, $n
        )
    end |> esc
end

macro domain_check(m, n)
    esc(:(@domain_check($n ≥ 0 && μ ≤ $n && iseven($n - μ), $m, $n)))
end

macro domain_check(j)
    esc(:(@domain($j ≥ 0, $j)))
end
