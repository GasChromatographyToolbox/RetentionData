using Documenter
using CSV
using RetentionData

makedocs(
    sitename = "RetentionData",
    #format = Documenter.HTML(),
    #modules = [RetentionData]
    pages = Any[
                "Home" => "index.md",
                "File structure" => "filestructure.md",
                "Naming convention" => "nameconvention.md",
                "Docstrings" => "docstrings.md",
                "Supplemental Materials" => "suppmat.md"
            ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/GasChromatographyToolbox/RetentionData",
    devbranch = "main"
)
