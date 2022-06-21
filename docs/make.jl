using Documenter
using TestRay

makedocs(
    sitename = "TestRay",
    format = Documenter.HTML(),
    modules = [TestRay]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
