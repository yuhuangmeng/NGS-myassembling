from pymetis import part_graph
from ngsolve import ElementId, VOL


def metis_partition_from_fes(fes, nparts: int, overlap_width: int = 0, free_dofs_only: bool = False):
    """Partition volume elements using a shared-DoF element graph.

    Returns ``core_partition, overlapping_partition, cutcount``. The optional
    ``overlap_width`` expands METIS parts in the element graph and is separate
    from exact support-closure construction.
    """
    if overlap_width < 0:
        raise ValueError("overlap_width must be >= 0")

    mesh = fes.mesh
    elements = list(mesh.Elements(VOL))
    ne = len(elements)
    freedofs = fes.FreeDofs() if free_dofs_only else None

    dof_to_positions: dict[int, list[int]] = {}
    for pos, el in enumerate(elements):
        for d in fes.GetDofNrs(ElementId(VOL, int(el.nr))):
            d = int(d)
            if d < 0:
                continue
            if free_dofs_only and not freedofs[d]:
                continue
            dof_to_positions.setdefault(d, []).append(pos)

    adjacency_sets = [set() for _ in range(ne)]
    for positions in dof_to_positions.values():
        for a, i in enumerate(positions):
            for j in positions[a + 1:]:
                adjacency_sets[i].add(j)
                adjacency_sets[j].add(i)

    adjacency = [sorted(neighbors) for neighbors in adjacency_sets]
    cutcount, part_ids = part_graph(nparts, adjacency=adjacency)

    core_partition_pos = [[] for _ in range(nparts)]
    for pos, part_id in enumerate(part_ids):
        core_partition_pos[part_id].append(pos)

    overlapping_partition_pos = []
    for part in core_partition_pos:
        patch = set(part)
        for _ in range(overlap_width):
            expanded = set(patch)
            for pos in patch:
                expanded.update(adjacency_sets[pos])
            patch = expanded
        overlapping_partition_pos.append(sorted(patch))

    core_partition = [
        [int(elements[pos].nr) for pos in part]
        for part in core_partition_pos
    ]
    overlapping_partition = [
        [int(elements[pos].nr) for pos in part]
        for part in overlapping_partition_pos
    ]

    return core_partition, overlapping_partition, cutcount
