push!(LOAD_PATH, "../src")

using Documenter

makedocs(
    sitename = "Thistle Documentation",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Numerics" => "numerics.md",
    ],
)
