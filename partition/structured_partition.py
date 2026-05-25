
def structured_quad_element_partitions(Nx, Ny, nx, ny, overlap_width=0):
    """Build structured quad element partitions in NGSolve element numbering.

    ``Nx`` x ``Ny`` is the coarse subdomain grid. Each subdomain contains
    ``nx`` x ``ny`` fine quadrilateral elements. For ``MakeQuadMesh(nelx,
    nely)``, NGSolve numbers elements row by row:

        elnr = ix + iy * nelx

    The overlapping partition expands each structured element block by
    ``overlap_width`` element layers in the structured grid.
    """
    if overlap_width < 0:
        raise ValueError("overlap_width must be >= 0")

    nelx = Nx * nx
    nely = Ny * ny

    def elements_in_box(ix0, ix1, iy0, iy1):
        return [
            int(ix + iy * nelx)
            for iy in range(iy0, iy1 + 1)
            for ix in range(ix0, ix1 + 1)
        ]

    core_partition = []
    overlapping_partition = []

    for py in range(Ny):
        for px in range(Nx):
            ix0 = px * nx
            ix1 = (px + 1) * nx - 1
            iy0 = py * ny
            iy1 = (py + 1) * ny - 1

            core_partition.append(elements_in_box(ix0, ix1, iy0, iy1))

            ox0 = max(0, ix0 - overlap_width)
            ox1 = min(nelx - 1, ix1 + overlap_width)
            oy0 = max(0, iy0 - overlap_width)
            oy1 = min(nely - 1, iy1 + overlap_width)
            overlapping_partition.append(elements_in_box(ox0, ox1, oy0, oy1))

    return core_partition, overlapping_partition
