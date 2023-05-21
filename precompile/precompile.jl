using PackageCompiler

sysimage_path = "Zernike." * Libc.Libdl.dlext

precompile_execution_file = "precompile/precompile_script.jl"

create_sysimage(:Zernike; sysimage_path, precompile_execution_file)
