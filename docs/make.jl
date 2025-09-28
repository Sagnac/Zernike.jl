using Pkg
Pkg.activate("docs")
Pkg.instantiate()

using Documenter, Zernike

module MakeDocs

const api_links = [
    "Zernike\\.metrics",
    "Zernike\\.grad",
    "Zernike\\.lap",
    "plotconfig",
    "zplot",
    "format_strings",
    "print_strings",
    "radial_coefficients",
    "transform_coefficients",
]

const api_regex = Regex("(`(" * join(api_links, "|(Zernike\\.)??") * ").*?`)")

const link_paths = Dict(
    "docs/src/assets/images/image.png" => "assets/images/image.png",
    r".*github\.io.*"s => "",
    api_regex => s"[\1](@ref \2)",
    "(precompile)" => "(https://github.com/Sagnac/Zernike.jl/tree/master/precompile)"
)

const primaries = ["zernike", "wavefront", "transform"]

const prim_regex = Regex(join(primaries, "|"))

const meta_block = """
```@meta
CurrentModule = Zernike
DocTestSetup = :(using Zernike)
EditURL = "https://github.com/Sagnac/Zernike.jl/blob/master/README.md"
```
"""

check(m) = m isa RegexMatch

function trimlines(file, io)
    close(io)
    text = read(file, String)
    text = replace(text, r"(\n)\n\n+$" => s"\1", r"(\n\n)\n+" => s"\1")
    write(file, text)
    return
end

function readme(path)
    pages = Any["Home" => "index.md"]
    section = []
    file = path * pages[1][2]
    write(file, meta_block)
    io = open(file, "a")
    for line in eachline("README.md"; keep = true)
        line = replace(line, link_paths...)
        m = match(r"(##+ )(?<header>.+)", line)
        matched = check(m)
        if matched && m.captures[1] == "## "
            trimlines(file, io)
            header = m[:header]
            src_page = replace(header, r"[\\/:*?\"<>|]" => ",", "`" => "") * ".md"
            file = path * src_page
            subsection = check(match(prim_regex, header))
            push!(subsection ? section : pages, src_page)
            if length(section) == length(primaries)
                push!(pages, "Primary functions" => section)
                empty!(primaries)
            end
            write(file, "# " * header * "\n" * meta_block)
            io = open(file, "a")
            continue
        end
        write(io, replace(line, (matched ? "#" : r"-{3,}") => ""; count = 1))
    end
    trimlines(file, io)
    push!(pages, "API" => "api.md")
    return pages
end

const branch = get(ENV, "GITHUB_REF_NAME", "")

end # module MakeDocs

DocMeta.setdocmeta!(Zernike, :DocTestSetup, :(using Zernike); recursive = true)

makedocs(
    sitename = "Zernike.jl",
    pages = MakeDocs.readme("docs/src/"),
    modules = [Zernike],
    format = Documenter.HTML(warn_outdated = false)
)

if !isempty(MakeDocs.branch)
    deploydocs(
        repo = "github.com/Sagnac/Zernike.jl.git",
        devbranch = MakeDocs.branch,
        devurl = MakeDocs.branch,
        versions = ["stable" => "master", "v^", "dev" => "dev"]
    )
end
