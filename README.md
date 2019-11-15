# optimal_covering

An attempt to solve the following problem:

We have a region on the Earth we want to observe with a satellite in a finite
number of observations. Where should we look?

Or more generally, we have some shape we want to cover with a finite number of
other, smaller shapes. Where should we place those smaller shapes?

# Getting Started

First, fork this repository. Or use Docker, [TODO: however exactly that's done]

All dependencies are listed in the `environment.yml` file.

# Python Modules

## Ones you should take a look at:

- `tests.py`: run this first, to test that everything's working well
- `genetics.py`: performs genetic algorithm

## Hopefully you shouldn't need to open these:

- `constants.py`: stores constant parameters / paths / ...
- `earth_data.py`: downloads earth data, finds observable spots
- `plot_funcs.py`: plotting functions
- `geometry.py`: deals with satellite / Earth system geometry
- `sampling.py`: a few functions related to random sampling
- `timer.py`: a handy-dandy timer function to use as a decorator

# Important Notes

- For consistency, always write geometric coordinates as (lon, lat)
- Recall that longitudes range from -180 to 180, latitudes from -90 to 90



