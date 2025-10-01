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

macro domain_check_mn()
    quote
        @domain(n ≥ 0 && μ ≤ n && iseven(n - μ),
            """
            \nDomain:
            n ≥ 0
            |m| ≤ n
            n ≡ m (mod 2)
            """,
            m, n
        )
    end |> esc
end

macro domain_check_j()
    esc(:(@domain(j ≥ 0, j)))
end

macro domain_check_mn_max()
    quote
        @domain(m_max ≥ 0 && n_max ≥ 0 && m_max ≤ n_max && iseven(n_max - m_max),
            """
            \nDomain:
            m_max ≥ 0
            n_max ≥ 0
            m_max ≤ n_max
            n_max ≡ m_max (mod 2)
            """,
            m_max, n_max
        )
    end |> esc
end
