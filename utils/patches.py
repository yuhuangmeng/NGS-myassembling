from dataclasses import dataclass

from ngsolve import ElementId, GridFunction, L2, VOL
from ngsolve.webgui import Draw


@dataclass
class LocalSupportPatch:
    """Index data for exact support-closure local assembly."""

    core_elements: list[int]
    support_elements: list[int]
    core_dofs: list[int]
    support_dofs: list[int]
    core_in_support: list[int]


def element_dofs(fes, elnr: int) -> list[int]:
    """Return non-negative global dofs on one volume element."""
    dnums = fes.GetDofNrs(ElementId(VOL, int(elnr)))
    return [int(d) for d in dnums if int(d) >= 0]


def dofs_of_elements(fes, elements) -> list[int]:
    """Return sorted unique non-negative global dofs on elements."""
    dofs = set()
    for elnr in elements:
        dofs.update(element_dofs(fes, int(elnr)))
    return sorted(dofs)


def build_dof_to_elements(fes) -> dict[int, list[int]]:
    """Map each global dof to all volume elements whose dof list contains it."""
    dof_to_elements: dict[int, list[int]] = {}
    for el in fes.mesh.Elements(VOL):
        elnr = int(el.nr)
        for d in element_dofs(fes, elnr):
            dof_to_elements.setdefault(d, []).append(elnr)
    return dof_to_elements


def build_support_patch(fes, core_elements) -> LocalSupportPatch:
    """Build the exact support closure of a local cell partition.

    ``core_dofs`` are all global DoFs appearing on ``core_elements``.
    ``support_elements`` are all volume elements whose DoF list intersects
    ``core_dofs``. Reassembling on this support closure exactly reproduces the
    global core-core matrix block:

        A_support[core_in_support, core_in_support]
            == A_global[core_dofs, core_dofs]

    up to roundoff. The full support matrix is not required to match the
    algebraic global submatrix on ``support_dofs``. The same support closure
    also gives the corresponding exact core entries for a linear right-hand
    side:

        b_support[core_in_support] == b_global[core_dofs]

    The same support closure is used in the nonlinear verification notebook for
    restricted residual and Jacobian checks. There the nonlinear objects are
    still global-size NGSolve forms, restricted to ``support_elements`` for
    verification.
    """
    core_elements = sorted({int(elnr) for elnr in core_elements})
    core_dofs = dofs_of_elements(fes, core_elements)

    dof_to_elements = build_dof_to_elements(fes)
    support_element_set = set()
    for d in core_dofs:
        support_element_set.update(dof_to_elements.get(d, []))

    support_elements = sorted(support_element_set)
    support_dofs = dofs_of_elements(fes, support_elements)
    support_index = {d: i for i, d in enumerate(support_dofs)}
    core_in_support = [support_index[d] for d in core_dofs]

    for i, d in enumerate(core_dofs):
        assert support_dofs[core_in_support[i]] == d

    return LocalSupportPatch(
        core_elements=core_elements,
        support_elements=support_elements,
        core_dofs=core_dofs,
        support_dofs=support_dofs,
        core_in_support=core_in_support,
    )


def build_support_patches(fes, partition) -> list[LocalSupportPatch]:
    """Build support closures for all parts in an element partition."""
    return [build_support_patch(fes, core_elements) for core_elements in partition]


def print_patch_summary(patch: LocalSupportPatch, max_entries: int = 20) -> None:
    """Print compact support-patch information for notebook debugging."""
    def preview(values):
        values = list(values)
        suffix = " ..." if len(values) > max_entries else ""
        return f"{values[:max_entries]}{suffix}"

    print("core_elements   =", len(patch.core_elements), preview(patch.core_elements))
    print("support_elements=", len(patch.support_elements), preview(patch.support_elements))
    print("core_dofs       =", len(patch.core_dofs), preview(patch.core_dofs))
    print("support_dofs    =", len(patch.support_dofs), preview(patch.support_dofs))
    print("core_in_support =", len(patch.core_in_support), preview(patch.core_in_support))


def draw_patch(mesh, patch: LocalSupportPatch, name: str = "support patch"):
    """Draw 0 outside, 1 core elements, and 2 support-only elements."""
    l2 = L2(mesh, order=0)
    gf = GridFunction(l2, name=name)
    gf.vec[:] = 0
    for elnr in patch.support_elements:
        gf.vec[l2.GetDofNrs(ElementId(VOL, elnr))[0]] = 2
    for elnr in patch.core_elements:
        gf.vec[l2.GetDofNrs(ElementId(VOL, elnr))[0]] = 1
    return Draw(gf, mesh, name)
