# Load the data
data(laja)
# Make a 'comparative community dataset' that merges the datatypes together
c.data <- comparative.comm(invert.tree, river.sites, invert.traits, river.env)

# Format all the data to be loaded into stan
laja.pres <- as.numeric(river.sites)
as.numeric(matrix(1:4, 2))
laja.traits <- rep(c.data$data$length, each=length(species(c.data)))
laja.env <- rep(c.data$env$temperature, length(sites(c.data)))

# Run a stan model (replace '...' with all the other parts of a stan call)
stan(..., pres=laja.pres, traits=laja.traits, env=laja.env)
