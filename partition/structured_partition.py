def structured_quad_partitions(Nx, Ny, nx, ny, overlap_size=1):
    """Return structured non-overlapping and overlapping DoF partitions.

    This helper intentionally knows nothing about finite elements or NGSolve
    meshes. It only partitions row-major structured DoF indices.

    Assumptions for this first version:

    - structured quadrilateral grid,
    - scalar H1 order=1,
    - row-major structured DoF numbering

          dof = ix + iy * ndof_x

    - ``Nx`` x ``Ny`` is the structured subdomain grid,
    - ``nx`` x ``ny`` is the number of nodes/DoFs per non-overlapping
      subdomain,
    - ``overlap_size`` is the number of structured DoF-grid layers used to
      expand each non-overlapping block.

    Returns:

        nonoverlapping_partition, overlapping_partition

    Both are ``list[list[int]]``. The non-overlapping partition is disjoint and
    used for bookkeeping/visualization. The overlapping partition gives the
    actual local unknowns passed to C++.
    """
    if Nx <= 0 or Ny <= 0:
        raise ValueError("Nx and Ny must be positive")
    if nx <= 0 or ny <= 0:
        raise ValueError("nx and ny must be positive")
    if overlap_size < 0:
        raise ValueError("overlap_size must be >= 0")

    ndof_x = Nx * nx
    ndof_y = Ny * ny

    def box_dofs(ix0, ix1, iy0, iy1):
        return [
            int(ix + iy * ndof_x)
            for iy in range(iy0, iy1 + 1)
            for ix in range(ix0, ix1 + 1)
        ]

    nonoverlapping_partition = []
    overlapping_partition = []

    for py in range(Ny):
        for px in range(Nx):
            ix0 = px * nx
            ix1 = (px + 1) * nx - 1
            iy0 = py * ny
            iy1 = (py + 1) * ny - 1

            nonoverlapping_partition.append(box_dofs(ix0, ix1, iy0, iy1))

            ox0 = max(0, ix0 - overlap_size)
            ox1 = min(ndof_x - 1, ix1 + overlap_size)
            oy0 = max(0, iy0 - overlap_size)
            oy1 = min(ndof_y - 1, iy1 + overlap_size)
            overlapping_partition.append(box_dofs(ox0, ox1, oy0, oy1))

    return nonoverlapping_partition, overlapping_partition
