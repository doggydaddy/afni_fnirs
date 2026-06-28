"""
render_figure.py — dorsal two-hemisphere SUMA-style render (headless).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Optional

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm as mpl_cm
from matplotlib.colors import Normalize
import nibabel as nib
import numpy as np
from nilearn import datasets, plotting, surface

import warnings
warnings.filterwarnings(
    "ignore",
    message=".*non integer threshold but configured the colorbar.*",
    category=UserWarning,
)

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DEFAULT_BRAIN = REPO / "brain.nii"
DEFAULT_TEMPLATE_TXT = REPO / "fNIRS_template.txt"
DEFAULT_TEMPLATE_NII = REPO / "fNIRS_template.nii"


def _abs_max(img):
    d = np.asarray(img.dataobj, dtype=float)
    d = d[np.isfinite(d)]
    if d.size == 0:
        return 1.0
    m = float(np.max(np.abs(d)))
    return m if m > 0 else 1.0


def _signed_threshold(img, eps=1e-12):
    d = np.asarray(img.dataobj, dtype=float)
    nz = np.abs(d[np.isfinite(d) & (d != 0)])
    if nz.size == 0:
        return eps
    return max(float(nz.min()) * 0.99, eps)


def channel_mni_centres(template_txt=DEFAULT_TEMPLATE_TXT,
                        template_nii=DEFAULT_TEMPLATE_NII):
    """Return {channel: MNI centre (x, y, z)} averaged over the 2 template voxels."""
    affine = nib.load(str(template_nii)).affine
    rows = np.loadtxt(template_txt)
    bucket = {}
    for i, j, k, ch in rows:
        bucket.setdefault(int(ch), []).append(np.array([i, j, k, 1.0]))
    out = {}
    for ch, ijks in bucket.items():
        ijk = np.mean(np.stack(ijks, axis=0), axis=0)
        out[ch] = (affine @ ijk)[:3]
    return out


def significant_channels(overlay, centres, template_nii=DEFAULT_TEMPLATE_NII,
                         eps=1e-12):
    """Return {channel: value} for channels whose overlay value is non-zero."""
    tpl = nib.load(str(template_nii))
    inv = np.linalg.inv(tpl.affine)
    data = np.asarray(overlay.dataobj)
    out = {}
    for ch, mni in centres.items():
        ijk = inv @ np.array([mni[0], mni[1], mni[2], 1.0])
        i, j, k = np.round(ijk[:3]).astype(int)
        if (0 <= i < data.shape[0] and 0 <= j < data.shape[1]
                and 0 <= k < data.shape[2]):
            v = float(data[i, j, k])
            if abs(v) > eps:
                out[ch] = v
    return out


def _align_to_pial(pial_verts, other_verts):
    """
    Translate + uniform-scale ``other_verts`` so its centroid and RMS extent
    match ``pial_verts``. Returns a copy. Used to align the inflated mesh
    onto pial for blending — inflated fsaverage is each-hemi-at-origin with
    ~40 mm radius, while pial is in MNI-ish space, so they need rigid+scale
    alignment before any per-vertex interpolation makes geometric sense.
    """
    pc = pial_verts.mean(axis=0)
    oc = other_verts.mean(axis=0)
    ps = np.sqrt(((pial_verts - pc) ** 2).sum(axis=1).mean())
    os_ = np.sqrt(((other_verts - oc) ** 2).sum(axis=1).mean())
    if os_ < 1e-9:
        return pial_verts.copy()
    return (other_verts - oc) * (ps / os_) + pc


def _rotate_z_around_pivot(verts, theta_deg, pivot_xy):
    """
    Rotate vertices around a vertical (Z) axis that passes through
    ``pivot_xy = (px, py)``, by ``theta_deg`` degrees. Positive theta is
    counter-clockwise when looking DOWN from +Z (right-hand rule).

    Used to swing the two hemispheres around the frontal pole so the back
    fans out and the dorsal-lateral cortex rotates toward an anterior
    camera, like opening a book hinged at the front.
    """
    if theta_deg == 0:
        return verts
    theta = np.deg2rad(theta_deg)
    c, s = np.cos(theta), np.sin(theta)
    px, py = pivot_xy
    # Row-vector convention: verts_new = (verts - pivot) @ R.T + pivot
    R_T = np.array([[c,  s, 0],
                    [-s, c, 0],
                    [0,  0, 1]])
    pivot = np.array([px, py, 0.0])
    return (verts - pivot) @ R_T + pivot


def _push_out_radial(point, brain_centre=(0.0, 0.0, 0.0), distance=10.0):
    """
    Move a point radially outward from ``brain_centre`` by ``distance`` mm.

    Kept for back-compat / simple cases. For accurate scalp-like placement
    that respects the cortical surface, use ``_vertex_normals`` +
    ``_channel_scalp_positions`` instead.
    """
    p = np.asarray(point, dtype=float)
    c = np.asarray(brain_centre, dtype=float)
    direction = p - c
    norm = np.linalg.norm(direction)
    if norm < 1e-9:
        return p
    return p + (distance / norm) * direction


def _vertex_normals(verts, faces):
    """
    Compute outward-pointing per-vertex normals for a triangular mesh.

    Each vertex normal is the (length-normalised) sum of the unit face
    normals of all incident faces. Any normal pointing inward (toward the
    mesh centroid) is flipped so the result is consistently outward.
    """
    v0 = verts[faces[:, 0]]
    v1 = verts[faces[:, 1]]
    v2 = verts[faces[:, 2]]
    face_n = np.cross(v1 - v0, v2 - v0)
    face_n /= (np.linalg.norm(face_n, axis=1, keepdims=True) + 1e-12)

    vn = np.zeros_like(verts)
    np.add.at(vn, faces[:, 0], face_n)
    np.add.at(vn, faces[:, 1], face_n)
    np.add.at(vn, faces[:, 2], face_n)
    vn /= (np.linalg.norm(vn, axis=1, keepdims=True) + 1e-12)

    # Flip any normals that point inward.
    centroid = verts.mean(axis=0)
    out_dir = verts - centroid
    dots = (vn * out_dir).sum(axis=1)
    vn[dots < 0] *= -1.0
    return vn


def _channel_scalp_positions(centres, pial_verts, pial_normals,
                             render_verts, render_normals, pushout_mm=10.0,
                             search_radius_mm=20.0, outer_fraction=0.3):
    """
    For each channel centre in MNI space, return the cube position on the
    'scalp' just outside the render mesh.

    Algorithm
    ---------
    1. Collect all pial vertices within ``search_radius_mm`` of the channel
       centre (channel positions come from a pial-based template).
    2. Filter to the top ``outer_fraction`` of those by distance from the
       brain centroid — this keeps only vertices on the outer hull of the
       cortex (excluding medial walls and vertices tucked inside sulci).
    3. Among the survivors, pick the one closest to the channel centre.
       That's the cortex point directly under the fNIRS sensor.
    4. Look up the same vertex index in the *render* mesh and push out by
       ``pushout_mm`` along the render-mesh outward normal. fsaverage
       pial / inflated / white share one-to-one vertex correspondence, so
       the same index points to the same logical cortex location regardless
       of how the mesh has been deformed.

    The ``pial_normals`` argument is currently unused (kept for API
    stability).

    Returns ``{channel: (x, y, z)}`` in the render mesh coordinate frame.
    """
    centroid = pial_verts.mean(axis=0)
    out = {}
    for ch, mni in centres.items():
        mni = np.asarray(mni, dtype=float)
        diffs = pial_verts - mni
        d2 = (diffs ** 2).sum(axis=1)
        nearby = np.where(d2 < search_radius_mm ** 2)[0]
        if nearby.size == 0:
            nearby = np.array([int(np.argmin(d2))])

        # Keep only the most "outer" vertices (top fraction by distance
        # from brain centroid). These are the scalp-facing vertices —
        # excludes medial-wall and sulcus-bottom vertices.
        dist_centroid = np.linalg.norm(pial_verts[nearby] - centroid, axis=1)
        cutoff = np.quantile(dist_centroid, 1.0 - outer_fraction)
        outer_mask = dist_centroid >= cutoff
        outer = nearby[outer_mask]
        if outer.size == 0:
            outer = nearby

        # Among outer vertices, pick the one closest to the channel.
        diffs_outer = pial_verts[outer] - mni
        d2_outer = (diffs_outer ** 2).sum(axis=1)
        idx = int(outer[np.argmin(d2_outer)])

        out[ch] = render_verts[idx] + pushout_mm * render_normals[idx]
    return out


def _draw_cube(ax, center, side, face_color, alpha=1.0,
               edge_color="white", edge_width=1.2, fill=True):
    """
    Draw an axis-aligned cube centred at ``center`` (mm in MNI/render space)
    with the given side length. Uses ``Poly3DCollection`` so matplotlib
    handles depth sorting against the brain mesh.

    Parameters
    ----------
    ax : matplotlib 3-D Axes
    center : (cx, cy, cz)
    side : float
        Edge length in mm.
    face_color : RGB(A) tuple or matplotlib color string
        Face fill colour. Pass an alpha here for per-face transparency, or
        set ``alpha`` to apply to the whole collection.
    alpha : float
        Face alpha (matplotlib applies this on top of any alpha embedded
        in face_color).
    edge_color : matplotlib color
        Edge (wireframe) colour. Keep contrasting (white on black bg) so
        the cube reads clearly.
    edge_width : float
        Edge line width in points.
    fill : bool
        If False, draw only the wireframe (edges) — face_color is ignored.
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    cx, cy, cz = center
    s = side / 2.0
    # 8 cube vertices
    v = np.array([
        [cx - s, cy - s, cz - s], [cx + s, cy - s, cz - s],
        [cx + s, cy + s, cz - s], [cx - s, cy + s, cz - s],
        [cx - s, cy - s, cz + s], [cx + s, cy - s, cz + s],
        [cx + s, cy + s, cz + s], [cx - s, cy + s, cz + s],
    ])
    # 6 quad faces, vertex order chosen so outward normals point outward
    faces = [
        [v[0], v[1], v[2], v[3]],  # bottom (z = -s)
        [v[4], v[5], v[6], v[7]],  # top    (z = +s)
        [v[0], v[1], v[5], v[4]],  # front  (y = -s)
        [v[3], v[2], v[6], v[7]],  # back   (y = +s)
        [v[0], v[3], v[7], v[4]],  # left   (x = -s)
        [v[1], v[2], v[6], v[5]],  # right  (x = +s)
    ]
    if fill:
        fc = face_color
    else:
        fc = (1.0, 1.0, 1.0, 0.0)  # fully transparent fill = wireframe only
    pc = Poly3DCollection(faces, facecolors=fc, alpha=alpha,
                          edgecolors=edge_color, linewidths=edge_width)
    ax.add_collection3d(pc)


def _combine_hemispheres(fs, surface_kind="pial", hemisphere_gap=5.0,
                         inflate_alpha=0.0, hemisphere_tilt=0.0):
    """
    Load both hemispheres of the requested fsaverage surface and return a
    single mesh plus the concatenated sulc background.

    Parameters
    ----------
    fs : nilearn fsaverage bunch
    surface_kind : {"pial", "inflated", "white"}
        Which cortical surface to load. "inflated" exposes the channels that
        would otherwise be buried in sulci (good for fNIRS forehead-cap
        figures where channels 1/2/15/16 sit on the curving lateral wall).
    hemisphere_gap : float
        Gap (mm) between the two hemispheres along the X axis. Used whenever
        the natural separation isn't sufficient (inflated meshes, tilted /
        puffed pial meshes); raw pial is left as-is.
    inflate_alpha : float in [0, 1]
        Only meaningful for ``surface_kind="pial"``. Linearly blend the
        pial vertices toward the inflated mesh (aligned-and-rescaled to
        match pial). 0.0 = pure pial (default); 0.3 ≈ "puffed pial" that
        fills shallow sulci so lateral-wall channels show better from an
        anterior view; 1.0 ≈ fully inflated (use ``surface_kind="inflated"``
        instead). Topology is preserved, so patch assignment indices remain
        valid.
    hemisphere_tilt : float
        Degrees to swing each hemisphere around a vertical (Z) axis through
        the brain's frontal pole (X=0, Y=Y_max). Positive tilt rotates LH
        clockwise and RH counter-clockwise looking down from above — the
        posterior ends fan OUT laterally while the anterior pole stays
        anchored. This rotates dorsal-lateral cortex (where fNIRS channels
        13-16 sit) toward an anterior camera. 5-15° usually exposes the
        dorsal-most channels without looking unnatural. 0 = no tilt
        (default). Topology is preserved.

    Returns
    -------
    mesh : (vertices, faces)
        Combined mesh; RH faces are re-indexed so they refer to the bottom
        half of ``vertices``.
    sulc : (n_vertices,) ndarray
        Concatenated sulc depth map (LH then RH).
    n_lh : int
        Number of vertices in the LH half (so callers can split textures
        again if needed).
    """
    surf_attr = {
        "pial": ("pial_left", "pial_right"),
        "inflated": ("infl_left", "infl_right"),
        "white": ("white_left", "white_right"),
    }
    if surface_kind not in surf_attr:
        raise ValueError(
            f"surface_kind must be one of {list(surf_attr)}, got {surface_kind!r}")
    lh_key, rh_key = surf_attr[surface_kind]
    lh_verts, lh_faces = surface.load_surf_mesh(getattr(fs, lh_key))
    rh_verts, rh_faces = surface.load_surf_mesh(getattr(fs, rh_key))

    # "Puff up" pial by blending with aligned inflated. Done per-hemisphere
    # so each lobe is puffed within its own coordinate frame.
    if surface_kind == "pial" and inflate_alpha > 0:
        lh_infl, _ = surface.load_surf_mesh(fs.infl_left)
        rh_infl, _ = surface.load_surf_mesh(fs.infl_right)
        lh_infl_aligned = _align_to_pial(lh_verts, lh_infl)
        rh_infl_aligned = _align_to_pial(rh_verts, rh_infl)
        a = float(inflate_alpha)
        lh_verts = (1.0 - a) * lh_verts + a * lh_infl_aligned
        rh_verts = (1.0 - a) * rh_verts + a * rh_infl_aligned

    # Swing each hemisphere around a vertical (Z) axis through the brain's
    # frontal pole (X=0, Y=Y_max), so posterior ends fan outward and the
    # anterior pole stays anchored. LH gets -theta (CW from above), RH gets
    # +theta (CCW from above).
    if hemisphere_tilt != 0:
        y_max = float(max(lh_verts[:, 1].max(), rh_verts[:, 1].max()))
        pivot_xy = (0.0, y_max)
        lh_verts = _rotate_z_around_pivot(lh_verts, -hemisphere_tilt, pivot_xy)
        rh_verts = _rotate_z_around_pivot(rh_verts, +hemisphere_tilt, pivot_xy)

    # fsaverage inflated (and to a lesser extent white) hemispheres are each
    # centred around the origin, so concatenating them naively overlaps the
    # two blobs. Tilted / puffed pial can also push the dorsal medial edges
    # past the midline. In all these cases, shift X so the two hemispheres
    # abut with ``hemisphere_gap`` mm of breathing room.
    needs_shift = (surface_kind != "pial") or hemisphere_tilt != 0 \
                  or inflate_alpha > 0
    if needs_shift:
        lh_shift = -(lh_verts[:, 0].max() + hemisphere_gap / 2.0)
        rh_shift = -(rh_verts[:, 0].min() - hemisphere_gap / 2.0)
        lh_verts = lh_verts.copy()
        rh_verts = rh_verts.copy()
        lh_verts[:, 0] += lh_shift
        rh_verts[:, 0] += rh_shift

    n_lh = lh_verts.shape[0]
    verts = np.vstack([lh_verts, rh_verts])
    faces = np.vstack([lh_faces, rh_faces + n_lh])

    sulc_lh = surface.load_surf_data(fs.sulc_left)
    sulc_rh = surface.load_surf_data(fs.sulc_right)
    sulc = np.concatenate([sulc_lh, sulc_rh])
    return (verts, faces), sulc, n_lh


# Named camera presets (elev, azim) for matplotlib 3-D axes.
# Anterior pole points toward +y in fsaverage. azim=90 looks from +y toward
# the origin (i.e. from in front of the brain).
VIEW_PRESETS = {
    "anterior":  (0.0,   90.0),   # looking from front, superior at top
    "posterior": (0.0,  270.0),   # looking from back
    "dorsal":    (90.0, 270.0),   # top-down, anterior at top
    "ventral":   (-90.0, 90.0),   # bottom-up
    "left":      (0.0,  180.0),   # left lateral
    "right":     (0.0,    0.0),   # right lateral
}


def _resolve_view(view):
    """Accept either a preset name or an (elev, azim) tuple."""
    if isinstance(view, str):
        key = view.lower()
        if key not in VIEW_PRESETS:
            raise ValueError(
                f"Unknown view preset {view!r}. "
                f"Choose from: {sorted(VIEW_PRESETS)} or pass an (elev, azim) tuple."
            )
        return VIEW_PRESETS[key]
    elev, azim = view
    return float(elev), float(azim)


def channel_values_from_overlay(overlay, centres,
                                template_nii=DEFAULT_TEMPLATE_NII,
                                eps=1e-12):
    """
    Sample the overlay NIfTI at each channel's central voxel and return a dict
    ``{channel: value}`` (including zeros).
    """
    tpl = nib.load(str(template_nii))
    inv = np.linalg.inv(tpl.affine)
    data = np.asarray(overlay.dataobj)
    out = {}
    for ch, mni in centres.items():
        ijk = inv @ np.array([mni[0], mni[1], mni[2], 1.0])
        i, j, k = np.round(ijk[:3]).astype(int)
        if (0 <= i < data.shape[0] and 0 <= j < data.shape[1]
                and 0 <= k < data.shape[2]):
            out[ch] = float(data[i, j, k])
        else:
            out[ch] = 0.0
    return out


def patch_texture(mesh_vertices, centres, values, patch_radius=20.0,
                  only_anterior_to=None, return_assignment=False):
    """
    Build a per-vertex texture by Voronoi assignment to channel centres.

    Every cortical vertex within ``patch_radius`` mm of a channel centre is
    coloured with that channel's value; vertices that have no channel within
    range stay at 0 (background). The result looks like discrete patches on
    the cortex — matching the SUMA "all channels" reference figure.

    Parameters
    ----------
    mesh_vertices : (N, 3) ndarray
        Combined LH+RH cortical vertices in MNI/RAS space.
    centres : dict {int: (x, y, z)}
        MNI centres of each channel (see ``channel_mni_centres``).
    values : dict {int: float}
        Per-channel value to render (zeros are kept as zero and threshold-clipped
        away at render time).
    patch_radius : float
        Hard cut-off in mm. Vertices farther than this from any channel get 0.
        Use a value that visually matches your reference (~18-25 mm works well
        for the 16-channel forehead montage).
    only_anterior_to : float or None
        If given, mask out vertices with y < this value (in MNI mm). Use e.g.
        ``-10`` or ``0`` to ensure no spurious patches appear on parietal /
        occipital cortex when a channel centre happens to be near the midline.
    return_assignment : bool
        If True, also return ``{channel: ndarray of vertex indices}``.

    Returns
    -------
    tex : (N,) ndarray
        Per-vertex stat value.
    assignment : dict, only when ``return_assignment=True``
    """
    centres_arr = np.array([centres[c] for c in sorted(centres)])
    channels = sorted(centres.keys())

    # Squared distances to every channel centre, vectorised.
    diffs = mesh_vertices[:, None, :] - centres_arr[None, :, :]  # (N, K, 3)
    d2 = np.einsum("nkc,nkc->nk", diffs, diffs)                  # (N, K)
    nearest = np.argmin(d2, axis=1)                              # (N,)
    nearest_d = np.sqrt(np.take_along_axis(d2, nearest[:, None],
                                           axis=1).ravel())      # (N,)

    in_range = nearest_d <= patch_radius
    if only_anterior_to is not None:
        in_range &= (mesh_vertices[:, 1] >= only_anterior_to)

    tex = np.zeros(mesh_vertices.shape[0], dtype=float)
    assignment = {ch: [] for ch in channels}
    for vi in np.where(in_range)[0]:
        ch = channels[nearest[vi]]
        tex[vi] = values.get(ch, 0.0)
        assignment[ch].append(vi)
    if return_assignment:
        assignment = {ch: np.asarray(idxs, dtype=int)
                      for ch, idxs in assignment.items()}
        return tex, assignment
    return tex


def _patch_label_positions(mesh_vertices, assignment, centres, cam_dir,
                           label_offset_mm=5.0):
    """
    For each channel that has any vertices in its patch, return a label
    position that sits on the **visible** patch surface (the camera-facing
    centroid of the patch), then bumped slightly toward the camera so the
    text isn't z-occluded by the cortex.

    Channels with no visible patch fall back to their MNI centre.
    """
    cam = np.asarray(cam_dir, dtype=float)
    cam = cam / np.linalg.norm(cam)
    out = {}
    for ch, idxs in assignment.items():
        if idxs.size == 0:
            out[ch] = np.asarray(centres[ch], dtype=float)
            continue
        pts = mesh_vertices[idxs]                       # (M, 3)
        # Keep only the half of the patch that faces the camera — i.e. points
        # whose component along cam_dir is above the patch median. This
        # avoids labels landing on the back side of a wrap-around patch.
        proj = pts @ cam
        keep = proj >= np.median(proj)
        front = pts[keep] if keep.any() else pts
        centroid = front.mean(axis=0)
        out[ch] = centroid + label_offset_mm * cam
    return out


def render_brain(overlay_path, out_png, title="", cmap="cold_hot",
                 vmax=None, threshold=None,
                 view="anterior",
                 channel_labels=False, label_only_significant=True,
                 label_offset_mm=3.0,
                 surf_mesh="fsaverage",
                 surface="inflated",
                 inflate_alpha=0.0,
                 hemisphere_tilt=0.0,
                 show_sulc=True,
                 brain_alpha=1.0,
                 surface_mode="patch", patch_radius=20.0,
                 patch_only_anterior_to=-20.0,
                 vol_to_surf_radius=0.0,
                 cube_size=15.0, cube_pushout=12.0,
                 cube_fill=True, cube_alpha=1.0, cube_edge_width=1.2,
                 figsize=(8.0, 8.0),
                 template_txt=DEFAULT_TEMPLATE_TXT,
                 template_nii=DEFAULT_TEMPLATE_NII):
    """
    SUMA-style render of an MNI volumetric overlay on a combined LH+RH pial
    brain. Single 3-D axes, black background, symmetric ±vmax colorbar.

    Parameters
    ----------
    view : str or (elev, azim) tuple
        Camera angle. Common presets: ``"anterior"`` (default — front view,
        superior at top), ``"dorsal"`` (top-down, anterior at top),
        ``"posterior"``, ``"ventral"``, ``"left"``, ``"right"``. Pass a tuple
        like ``(elev, azim)`` for arbitrary angles.
    surface_mode : {"patch", "vol_to_surf"}
        How to map the overlay onto the cortex.

        * ``"patch"`` (default) — discrete Voronoi patches around each
          channel centre with hard cut-off ``patch_radius`` mm. Reproduces
          the "all-channels coloured-tiles" SUMA figure.
        * ``"vol_to_surf"`` — nilearn's radial projection (smooth, blob-like).
          Use ``vol_to_surf_radius`` to control smoothing.

    patch_radius : float
        Cut-off radius for patch mode (mm). 18-25 typically reproduces the
        reference montage.
    patch_only_anterior_to : float or None
        Mask vertices with y < this value in patch mode. Prevents stray
        patches from appearing on the back of the brain when a channel
        centre is close to the dorsal midline. ``None`` disables the mask.
    channel_labels : bool
        Default ``False``. When True, overlay channel-number text on the blobs.
    """
    elev, azim = _resolve_view(view)

    overlay = nib.load(str(overlay_path))
    if vmax is None:
        vmax = _abs_max(overlay)
    if threshold is None:
        threshold = _signed_threshold(overlay)

    fs = datasets.fetch_surf_fsaverage(mesh=surf_mesh)
    # Always load the pial mesh too: it's the one whose vertex coordinates
    # match MNI/RAS, so patch Voronoi assignment and label centroid lookups
    # have to be computed there. The DISPLAY mesh (``mesh``) can be inflated
    # (or pial blended toward inflated via ``inflate_alpha``).
    pial_mesh, _, _ = _combine_hemispheres(fs, surface_kind="pial")
    mesh, sulc, _ = _combine_hemispheres(
        fs, surface_kind=surface,
        inflate_alpha=inflate_alpha,
        hemisphere_tilt=hemisphere_tilt)

    patch_assignment = None
    centres = channel_mni_centres(template_txt, template_nii)
    values = channel_values_from_overlay(overlay, centres, template_nii)
    if surface_mode == "patch":
        # Assignment is done on the PIAL vertices (MNI space). The returned
        # vertex indices are then valid for ``mesh`` too, because fsaverage
        # pial / inflated / white share one-to-one vertex correspondence.
        tex, patch_assignment = patch_texture(
            pial_mesh[0], centres, values,
            patch_radius=patch_radius,
            only_anterior_to=patch_only_anterior_to,
            return_assignment=True)
    elif surface_mode == "vol_to_surf":
        tex_lh = surface.vol_to_surf(overlay, fs.pial_left,
                                     radius=vol_to_surf_radius, kind="auto")
        tex_rh = surface.vol_to_surf(overlay, fs.pial_right,
                                     radius=vol_to_surf_radius, kind="auto")
        tex = np.concatenate([tex_lh, tex_rh])
    elif surface_mode == "cubes":
        # No per-vertex texture — the brain mesh is rendered plain grey and
        # we add 3-D cubes at each channel location below.
        tex = None
    else:
        raise ValueError(f"Unknown surface_mode {surface_mode!r}")

    fig = plt.figure(figsize=figsize, facecolor="black")
    gs = fig.add_gridspec(nrows=1, ncols=2,
                          width_ratios=[1.0, 0.04],
                          wspace=0.02, left=0.02, right=0.93,
                          top=0.92, bottom=0.05)
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    cax = fig.add_subplot(gs[0, 1])
    ax.set_facecolor("black")

    cmap_obj = (plotting.cm.cold_hot if cmap == "cold_hot"
                else mpl_cm.get_cmap(cmap))

    # Cube positions in the render frame are pre-computed: for each
    # channel, find the nearest pial vertex *that faces the same way*
    # (filters out medial-wall vertices for midline channels) then push
    # along the *render mesh's* outward normal at that vertex. Works for
    # any surface_kind / tilt / inflate.
    cube_positions = None
    if surface_mode == "cubes":
        pial_normals = _vertex_normals(pial_mesh[0], pial_mesh[1])
        render_normals = _vertex_normals(mesh[0], mesh[1])
        cube_positions = _channel_scalp_positions(
            centres, pial_mesh[0], pial_normals,
            mesh[0], render_normals,
            pushout_mm=cube_pushout)

    if surface_mode == "cubes":
        # Render the brain as a plain grey backdrop (no stat overlay).
        # Trick: pass a zero stat_map with threshold > vmax → nothing is
        # drawn from the stat layer, only the bg_map (sulc) shading shows.
        plotting.plot_surf_stat_map(
            mesh, stat_map=np.zeros(mesh[0].shape[0]),
            view=(elev, azim),
            bg_map=(sulc if show_sulc else None),
            bg_on_data=show_sulc, alpha=brain_alpha,
            cmap="Greys", vmax=1.0, threshold=2.0,
            colorbar=False, axes=ax, figure=fig)

        # Draw a cube at each pre-computed scalp position.
        norm = Normalize(vmin=-vmax, vmax=vmax)
        for ch in sorted(centres):
            val = values.get(ch, 0.0)
            pos = cube_positions[ch]
            if abs(val) < 1e-12:
                # Zero / sub-threshold channel — render as a faint wireframe
                # so the montage is still visible but doesn't distract.
                _draw_cube(ax, pos, side=cube_size,
                           face_color=(0.5, 0.5, 0.5, 0.0),
                           alpha=0.0,
                           edge_color=(0.6, 0.6, 0.6, 0.4),
                           edge_width=cube_edge_width, fill=False)
            else:
                rgba = cmap_obj(norm(val))
                _draw_cube(ax, pos, side=cube_size,
                           face_color=rgba, alpha=cube_alpha,
                           edge_color="white",
                           edge_width=cube_edge_width, fill=cube_fill)
    else:
        plotting.plot_surf_stat_map(
            mesh, stat_map=tex,
            view=(elev, azim),
            bg_map=(sulc if show_sulc else None),
            bg_on_data=show_sulc, alpha=1.0,
            cmap=cmap, vmax=vmax, threshold=threshold,
            colorbar=False, axes=ax, figure=fig)

    # Orthographic projection (no perspective foreshortening). SUMA renders
    # are orthographic, and it keeps the dorsal-most channels fully visible
    # even from a pure anterior view.
    try:
        ax.set_proj_type("ortho")
    except Exception:
        pass  # very old matplotlib — fall back to default persp

    cmap_obj = (plotting.cm.cold_hot if cmap == "cold_hot"
                else mpl_cm.get_cmap(cmap))
    sm = mpl_cm.ScalarMappable(norm=Normalize(vmin=-vmax, vmax=vmax),
                               cmap=cmap_obj)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.outline.set_edgecolor("white")
    cax.yaxis.set_tick_params(color="white", labelcolor="white")
    for s in cax.spines.values():
        s.set_edgecolor("white")
    ticks = np.linspace(-vmax, vmax, 5)
    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{t:.2f}" if abs(t) > 1e-6 else "0" for t in ticks])

    if title:
        fig.suptitle(title, color="white", fontsize=16, fontweight="bold")

    if channel_labels:
        # ``centres`` and ``values`` were already loaded above for the
        # texture / cube placement step.
        visible = (significant_channels(overlay, centres, template_nii)
                   if label_only_significant
                   else {ch: 0.0 for ch in centres})
        # Camera-pointing unit vector — labels are nudged this far toward
        # the camera so they don't get z-occluded by the cortex they sit on.
        elev_rad = np.deg2rad(elev)
        azim_rad = np.deg2rad(azim)
        cam_dir = np.array([
            np.cos(elev_rad) * np.cos(azim_rad),
            np.cos(elev_rad) * np.sin(azim_rad),
            np.sin(elev_rad),
        ])
        # Pick label positions per surface_mode:
        #   - cubes: re-use the pre-computed scalp positions but push a bit
        #     further along the same outward direction so the label sits
        #     just past the cube's outer face.
        #   - patch: centroid of the camera-facing half of the patch.
        #   - vol_to_surf: MNI centre + camera-direction offset.
        if surface_mode == "cubes":
            label_pushout = cube_pushout + cube_size / 2.0 + label_offset_mm
            label_pos = _channel_scalp_positions(
                centres, pial_mesh[0],
                _vertex_normals(pial_mesh[0], pial_mesh[1]),
                mesh[0], _vertex_normals(mesh[0], mesh[1]),
                pushout_mm=label_pushout)
        elif patch_assignment is not None:
            label_pos = _patch_label_positions(
                mesh[0], patch_assignment, centres, cam_dir,
                label_offset_mm=label_offset_mm)
        else:
            label_pos = {ch: np.asarray(centres[ch], dtype=float)
                         + label_offset_mm * cam_dir
                         for ch in centres}
        for ch in centres:
            if ch not in visible:
                continue
            x, y, z = label_pos[ch]
            ax.text(x, y, z, str(ch),
                    color="white", fontsize=11, fontweight="bold",
                    ha="center", va="center", zorder=100)

    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, facecolor="black", bbox_inches="tight")
    plt.close(fig)
    return out_png


# Back-compat alias so older callers (make_figures.py) keep working.
render_dorsal_pair = render_brain


def render_all(overlay_path, out_dir, stem=None, title="", vmax=None,
               view="anterior", channel_labels=False,
               bg_path=DEFAULT_BRAIN):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = stem or Path(overlay_path).stem
    return [render_brain(Path(overlay_path),
                         out_dir / f"{stem}.{view}.png",
                         title=title, vmax=vmax, view=view,
                         channel_labels=channel_labels)]


def _build_parser():
    p = argparse.ArgumentParser(
        description="Render an overlay NIfTI on a 3-D MNI brain (headless).")
    p.add_argument("--overlay", required=True, type=Path)
    p.add_argument("--out", type=Path)
    p.add_argument("--out-dir", type=Path)
    p.add_argument("--all", action="store_true")
    p.add_argument("--title", default="")
    p.add_argument(
        "--view", default="anterior",
        help="Camera preset (anterior, posterior, dorsal, ventral, left, "
             "right) or 'elev,azim' tuple, e.g. '30,120'. "
             "Default: anterior (front view, superior at top).",
    )
    p.add_argument("--vmax", type=float, default=None)
    p.add_argument(
        "--surface", choices=["inflated", "pial", "white"],
        default="inflated",
        help="Which fsaverage surface to render on. 'inflated' (default) "
             "smooths out the gyri/sulci so channels that sit on the lateral "
             "wall or inside a sulcus are still visible from an anterior "
             "view. 'pial' is anatomically faithful but hides those "
             "channels. Patch assignment is always done on pial.",
    )
    p.add_argument(
        "--inflate-alpha", type=float, default=0.0,
        help="Only with --surface pial: blend pial vertices toward the "
             "inflated mesh by this much, 0..1. 0.0 = raw pial (default), "
             "0.3 = 'puffed pial' that fills shallow sulci so lateral-wall "
             "channels (1/2/15/16) show better from anterior, while "
             "preserving the pial look. 1.0 = fully inflated.",
    )
    p.add_argument(
        "--hemisphere-tilt", type=float, default=0.0,
        help="Degrees to swing each hemisphere around a vertical axis "
             "through the frontal pole. Posterior fans outward, anterior "
             "stays anchored — exposes the dorsal-lateral channels "
             "(13-16) from an anterior view. Try 5-10°. Default 0.",
    )
    p.add_argument(
        "--no-sulc", action="store_true",
        help="Disable sulcal-depth shading. The default sulc bg darkens "
             "the inside of every sulcus, which makes the pial mesh look "
             "extra bumpy. With --no-sulc the cortex is rendered with a "
             "flat tone, which often looks 'fuller'.",
    )
    p.add_argument(
        "--surface-mode", choices=["patch", "vol_to_surf", "cubes"],
        default="patch",
        help="How to map the overlay onto the cortex. 'patch' (default) "
             "draws discrete Voronoi tiles around each channel centre "
             "(matches the SUMA reference). 'vol_to_surf' uses nilearn's "
             "smooth radial projection. 'cubes' draws a 3-D wire/solid "
             "cube at each channel position, floating just outside the "
             "cortex — useful when cortex curvature hides the lateral "
             "channels.",
    )
    p.add_argument(
        "--patch-radius", type=float, default=20.0,
        help="Patch mode: max distance (mm) from a channel centre that a "
             "cortical vertex is allowed to be coloured. Default 20.",
    )
    p.add_argument(
        "--patch-anterior-to", type=float, default=-20.0,
        help="Patch mode: only colour vertices with MNI y >= this. Prevents "
             "stray tiles on the back of the brain. Pass 'none' to disable.",
    )
    p.add_argument(
        "--vol-to-surf-radius", type=float, default=0.0,
        help="vol_to_surf mode only: smoothing radius (mm). 0 = sharpest.",
    )
    p.add_argument(
        "--cube-size", type=float, default=15.0,
        help="Cubes mode: edge length of each cube in mm. Default 15. "
             "Channel spacing is ~25 mm, so 15 keeps neighbours separated.",
    )
    p.add_argument(
        "--cube-pushout", type=float, default=12.0,
        help="Cubes mode: how far to push each cube outward from the brain "
             "centre along the radial direction. 12 mm puts the cubes just "
             "outside the cortex (~scalp depth). Default 12.",
    )
    p.add_argument(
        "--cube-wireframe", action="store_true",
        help="Cubes mode: draw only edges (wireframe), no filled faces. "
             "By default cubes are filled with the channel's value colour.",
    )
    p.add_argument(
        "--cube-alpha", type=float, default=1.0,
        help="Cubes mode: face alpha for filled cubes. Lower if cubes "
             "overlap or obscure each other. Default 1.0.",
    )
    p.add_argument(
        "--cube-edge-width", type=float, default=1.2,
        help="Cubes mode: edge line width in points. Default 1.2.",
    )
    p.add_argument(
        "--brain-alpha", type=float, default=1.0,
        help="Alpha for the brain mesh itself. Use <1 (e.g. 0.6) in cubes "
             "mode if you want to see cubes that sit behind the cortex. "
             "Default 1.0.",
    )
    p.add_argument("--labels", action="store_true",
                   help="Overlay channel-number labels (default: off).")
    p.add_argument("--label-all", action="store_true",
                   help="With --labels, label all 16 channels not just significant ones.")
    p.add_argument(
        "--label-offset", type=float, default=3.0,
        help="mm to nudge each label toward the camera so it sits in front "
             "of the cortex. Patch mode uses this on top of patch-centroid "
             "snapping; smooth mode uses it on top of the channel's MNI "
             "centre. Default 3.",
    )
    p.add_argument(
        "--tilt", type=float, default=0.0,
        help="Extra degrees of elevation (camera tilt UP) on top of the "
             "chosen view. e.g. --view anterior --tilt 15 looks down at "
             "the brain from the front, which makes the dorsal-most "
             "channels (rows 13-16 of the cap) more visible. Default 0.",
    )
    return p


def _parse_view(view_str):
    """Accept a preset name or 'elev,azim' string."""
    s = view_str.strip()
    if "," in s:
        elev, azim = s.split(",", 1)
        return (float(elev), float(azim))
    return s


def main(argv=None):
    args = _build_parser().parse_args(argv)
    view = _parse_view(args.view)
    # Apply --tilt by resolving the view to (elev, azim) and adding to elev.
    if args.tilt != 0.0:
        elev, azim = _resolve_view(view)
        view = (elev + args.tilt, azim)
    anterior_to = (None if str(args.patch_anterior_to).lower() == "none"
                   else float(args.patch_anterior_to))
    common = dict(
        title=args.title, vmax=args.vmax, view=view,
        channel_labels=args.labels,
        label_only_significant=not args.label_all,
        label_offset_mm=args.label_offset,
        surface=args.surface,
        inflate_alpha=args.inflate_alpha,
        hemisphere_tilt=args.hemisphere_tilt,
        show_sulc=not args.no_sulc,
        brain_alpha=args.brain_alpha,
        surface_mode=args.surface_mode,
        patch_radius=args.patch_radius,
        patch_only_anterior_to=anterior_to,
        vol_to_surf_radius=args.vol_to_surf_radius,
        cube_size=args.cube_size,
        cube_pushout=args.cube_pushout,
        cube_fill=not args.cube_wireframe,
        cube_alpha=args.cube_alpha,
        cube_edge_width=args.cube_edge_width,
    )
    if args.all:
        if args.out_dir is None:
            print("--all requires --out-dir", file=sys.stderr)
            return 2
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        stem = Path(args.overlay).stem
        view_name = view if isinstance(view, str) else "custom"
        out = render_brain(Path(args.overlay),
                           out_dir / f"{stem}.{view_name}.png",
                           **common)
        print(f"[render_figure] wrote {out}")
        return 0
    if args.out is None:
        print("Either --all (with --out-dir) or --out is required.",
              file=sys.stderr)
        return 2
    out = render_brain(args.overlay, args.out, **common)
    print(f"[render_figure] wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
