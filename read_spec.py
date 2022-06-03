def read_spec(filename, nspec, coords, coord_sys='ga', bandwidth=12e6):
    """
    Inputs:
    - filename: Name of the output FITS file.
    - nspec: Number of spectra collected.
    - coords: Coordinates of thetarget.
        Format: (l/ra, b/dec)
    - coord_sys: Coordinate system used for ''coords''.
        Default is galactic coordinates. Can use galactic ('ga') 
        or equatorial ('eq') coordinate systems.
    Returns:
    -
