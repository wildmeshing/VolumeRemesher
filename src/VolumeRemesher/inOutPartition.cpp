// Inside/outside classification of the BSP cells with respect to the input surface(s).
//
// Two classifiers:
//  - CellClassifier::WINDING_NUMBER (the default): the winding number of each cell, computed exactly
//    on the BSP itself and, for inputs where that is undefined (open or inconsistently oriented
//    surfaces), the generalized winding number. Only this file and winding_number.cpp, LGPL.
//  - CellClassifier::GRAPH_CUT: the historical min-cut over the dual graph of the cells
//    (inOutPartition_graphcut.cpp). It depends on GCO/maxflow, downloaded on request, whose license is NOT
//    LGPL compatible (GCO: research use only; maxflow: GPL), so it is compiled only with
//    -DVOLUMEREMESHER_WITH_GRAPHCUT=ON.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include "BSP.h"
#include "parallel.h"
#include "winding_number.h"

namespace vol_rem {
// Label = 1 means INTERNAL
// Label = 0 means EXTERNAL

// Returns TRUE if face is part of the skin according to skin_colour
inline bool isSkinFace(const BSPface& face, uint32_t skin_colour)
{
    return face.colour & skin_colour;
}

inline void setInternalCell(BSPcell& c, uint32_t internal_label)
{
    c.place |= internal_label;
}

inline void setExternalCell(BSPcell& c)
{
    c.place = EXTERNAL;
}

// A vertex of cell `cell_ind` that is not on the plane of face `f` (its vertices are marked in
// vrts_visit while searching, and unmarked before returning). The plane is that of the first
// coplanar constraint of f; a cell with no vertex off it is degenerate. Returns UINT32_MAX if
// none is found, which the callers treat as an error.
uint32_t BSPcomplex::cell_vertex_off_face(const BSPface& f, uint64_t cell_ind)
{
    const uint32_t* cid = constraints_vrts.data() + f.coplanar_constraints[0] * 3;
    const genericPoint* pv1 = vertices[cid[0]];
    const genericPoint* pv2 = vertices[cid[1]];
    const genericPoint* pv3 = vertices[cid[2]];

    BSPcell& cell = cells[cell_ind];
    uint64_t num_cellEdges = UINT64_MAX;
    uint32_t num_cellVrts = count_cellVertices(cell, &num_cellEdges);
    vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
    list_cellVertices(cell, num_cellEdges, cell_vrts);

    for (uint64_t ei : f.edges)
        vrts_visit[edges[ei].vertices[0]] = vrts_visit[edges[ei].vertices[1]] = 1;

    uint32_t found = UINT32_MAX;
    for (uint32_t vi : cell_vrts)
        if (!vrts_visit[vi]) {
            const genericPoint* cv = vertices[vi];
            if (genericPoint::orient3D(*cv, *pv1, *pv2, *pv3)) {
                found = vi;
                break;
            }
        }

    for (uint64_t ei : f.edges)
        vrts_visit[edges[ei].vertices[0]] = vrts_visit[edges[ei].vertices[1]] = 0;
    return found;
}

// ---------------------------------------------------------------------------------------------
// Winding-number classifier
// ---------------------------------------------------------------------------------------------

// The constraints of `group` (input surface triangles and hole caps; never edge/point forcing
// or virtual constraints) that cover face `face_ind`: the face lies in their plane and inside
// them. Same tests as blackAB_or_white, which decided the face colour, applied to every
// candidate constraint instead of stopping at the first hit.
void BSPcomplex::covering_constraints(
    uint64_t face_ind,
    CONSTR_GROUP_T group,
    std::vector<uint32_t>& out)
{
    out.clear();
    const BSPface& face = faces[face_ind];

    // coplanar_constraints can list a constraint more than once (a face inherits the lists of
    // the faces it was split from); each constraint must count once.
    std::vector<uint32_t> candidates;
    for (uint32_t constr : face.coplanar_constraints)
        if (constr < first_edge_constraint && constraint_group[constr] == group)
            candidates.push_back(constr);
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    if (candidates.empty()) return;

    const int xyz = face_dominant_normal_component(face);

    // Sample the face at a point strictly inside it that is on the boundary of no candidate:
    // the constraints strictly containing that point are the ones covering the face (a face
    // covered by two coplanar constraints that share an edge is covered once, by the one the
    // sample falls in). The approximate barycenter is the first try; if it lies on the boundary
    // of the face or of a candidate, the midpoints between it and the face vertices are tried.
    double p[3];
    get_approx_faceBaricenterCoord(face, p);
    std::vector<uint32_t> face_vrts(face.edges.size(), 0);
    list_faceVertices(faces[face_ind], face_vrts);

    auto classify_at = [&](const explicitPoint3D& P) -> bool {
        if (!is_baricenter_inFace(face, P, xyz)) return false;
        out.clear();
        for (uint32_t constr : candidates) {
            const uint32_t constr_ID = 3 * constr;
            const int lpt = localizedPointInTriangle(
                P,
                *vertices[constraints_vrts[constr_ID]],
                *vertices[constraints_vrts[constr_ID + 1]],
                *vertices[constraints_vrts[constr_ID + 2]],
                xyz);
            if (lpt == 1) return false; // on the boundary of a candidate: ambiguous sample
            if (lpt == 2) out.push_back(constr);
        }
        return true;
    };

    if (classify_at(explicitPoint3D(p[0], p[1], p[2]))) return;
    for (uint32_t v : face_vrts) {
        double q[3];
        vertices[v]->getApproxXYZCoordinates(q[0], q[1], q[2]);
        const explicitPoint3D mid(0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]), 0.5 * (p[2] + q[2]));
        if (classify_at(mid)) return;
    }

    // No generic sample found: decide from the face vertices. A constraint covers the face if
    // some face vertex is strictly inside it; it does not if some face vertex is strictly
    // outside it; when every face vertex is on its boundary, the constraint covers the face iff
    // its interior meets the face.
    // state: 0 undecided, 1 covers, 2 does not cover
    out.clear();
    std::vector<uint8_t> state(candidates.size(), 0);
    for (uint32_t vid : face_vrts) {
        const genericPoint* face_pt = vertices[vid];
        for (size_t k = 0; k < candidates.size(); k++) {
            if (state[k]) continue;
            const uint32_t constr_ID = 3 * candidates[k];
            const uint32_t vid1 = constraints_vrts[constr_ID];
            const uint32_t vid2 = constraints_vrts[constr_ID + 1];
            const uint32_t vid3 = constraints_vrts[constr_ID + 2];
            if (vid == vid1 || vid == vid2 || vid == vid3) continue; // on the boundary
            const int lpt = localizedPointInTriangle(
                *face_pt,
                *vertices[vid1],
                *vertices[vid2],
                *vertices[vid3],
                xyz);
            if (lpt == 2)
                state[k] = 1;
            else if (lpt == 0)
                state[k] = 2;
        }
    }
    for (size_t k = 0; k < candidates.size(); k++) {
        if (state[k] == 2) continue;
        if (state[k] == 1) {
            out.push_back(candidates[k]);
            continue;
        }
        const uint32_t* constraint = constraints_vrts.data() + 3 * candidates[k];
        if (coplanar_constraint_innerIntersects_face(face.edges, constraint, xyz))
            out.push_back(candidates[k]);
    }
}

// Makes the orientation of the triangles `tris` (constraint indices) consistent for the purpose
// of the winding number: flip[c] is set for the constraints to be read with their orientation
// reversed. Adjacent triangles (sharing an edge that exactly two of them use) are made to
// induce opposite directions on it, propagating from one triangle of each edge-connected
// component; edges used by one triangle (boundary) or more than two (non-manifold) do not
// propagate. A component that is closed (no boundary edge) is then oriented outward, i.e. so
// that its signed volume is positive: overlapping closed components add up instead of
// cancelling, and the winding number inside a single closed component is +1. An open component
// keeps the orientation the input gave it, up to the consistency fixes.
//
// This changes nothing for an input that is already consistently oriented outward, which is
// the common case; for one that is not, it is what makes the winding number meaningful.
static void repair_orientation(
    const std::vector<uint32_t>& tris,
    const std::vector<uint32_t>& constraints_vrts,
    const std::vector<double>& coords,
    std::vector<uint8_t>& flip)
{
    // Undirected edge -> incident triangles (positions in tris)
    struct EdgeRef
    {
        uint32_t a, b; // a < b
        uint32_t tri; // position in tris
    };
    std::vector<EdgeRef> refs;
    refs.reserve(3 * tris.size());
    for (uint32_t i = 0; i < tris.size(); i++) {
        const uint32_t* v = constraints_vrts.data() + 3 * tris[i];
        for (int k = 0; k < 3; k++) {
            const uint32_t a = v[k], b = v[(k + 1) % 3];
            refs.push_back({std::min(a, b), std::max(a, b), i});
        }
    }
    std::sort(refs.begin(), refs.end(), [](const EdgeRef& x, const EdgeRef& y) {
        if (x.a != y.a) return x.a < y.a;
        if (x.b != y.b) return x.b < y.b;
        return x.tri < y.tri;
    });

    // Manifold adjacency: for each triangle, the triangles across its two-sided edges.
    std::vector<std::vector<uint32_t>> adj(tris.size());
    std::vector<uint8_t> has_boundary(tris.size(), 0);
    for (size_t i = 0; i < refs.size();) {
        size_t j = i;
        while (j < refs.size() && refs[j].a == refs[i].a && refs[j].b == refs[i].b) j++;
        if (j - i == 2) {
            adj[refs[i].tri].push_back(refs[i + 1].tri);
            adj[refs[i + 1].tri].push_back(refs[i].tri);
        } else if (j - i == 1)
            has_boundary[refs[i].tri] = 1;
        i = j;
    }

    // Does triangle t, read with its flip, traverse the edge (a, b) in that direction?
    auto traverses = [&](uint32_t t, uint32_t a, uint32_t b) {
        const uint32_t* v = constraints_vrts.data() + 3 * tris[t];
        bool forward =
            (v[0] == a && v[1] == b) || (v[1] == a && v[2] == b) || (v[2] == a && v[0] == b);
        return flip[tris[t]] ? !forward : forward;
    };
    auto shared_edge = [&](uint32_t t, uint32_t u, uint32_t& a, uint32_t& b) {
        const uint32_t* vt = constraints_vrts.data() + 3 * tris[t];
        const uint32_t* vu = constraints_vrts.data() + 3 * tris[u];
        for (int k = 0; k < 3; k++) {
            a = vt[k];
            b = vt[(k + 1) % 3];
            for (int m = 0; m < 3; m++)
                if ((vu[m] == a && vu[(m + 1) % 3] == b) || (vu[m] == b && vu[(m + 1) % 3] == a))
                    return;
        }
        a = b = UINT32_MAX;
    };

    std::vector<uint8_t> visited(tris.size(), 0);
    std::vector<uint32_t> component, queue;
    for (uint32_t seed = 0; seed < tris.size(); seed++) {
        if (visited[seed]) continue;
        component.clear();
        queue.clear();
        queue.push_back(seed);
        visited[seed] = 1;
        bool closed = true;
        for (size_t head = 0; head < queue.size(); head++) {
            const uint32_t t = queue[head];
            component.push_back(t);
            if (has_boundary[t]) closed = false;
            for (uint32_t u : adj[t]) {
                if (visited[u]) continue;
                uint32_t a, b;
                shared_edge(t, u, a, b);
                if (a == UINT32_MAX) continue;
                // Consistent neighbours traverse the shared edge in opposite directions.
                if (traverses(t, a, b) == traverses(u, a, b)) flip[tris[u]] ^= 1;
                visited[u] = 1;
                queue.push_back(u);
            }
        }
        if (!closed) continue;
        double volume = 0.0;
        for (uint32_t t : component) {
            const uint32_t* v = constraints_vrts.data() + 3 * tris[t];
            const double* a = &coords[3 * v[0]];
            const double* b = &coords[3 * v[flip[tris[t]] ? 2 : 1]];
            const double* c = &coords[3 * v[flip[tris[t]] ? 1 : 2]];
            volume += a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) +
                      a[2] * (b[0] * c[1] - b[1] * c[0]);
        }
        if (volume < 0.0)
            for (uint32_t t : component) flip[tris[t]] ^= 1;
    }
}

// The change of the winding number of `group` when crossing face `face_ind` from
// conn_cells[0] to conn_cells[1]: +1 for each covering constraint whose positive side (in the
// sense of orient3D of a point against the constraint's vertices) holds conn_cells[0], -1 for
// each with conn_cells[0] on its negative side (reversed for the constraints flagged in
// `flip`, see repair_orientation). Two coincident constraints of opposite orientation cancel,
// two of the same orientation count twice, exactly as in the winding number of the input
// surface.
int BSPcomplex::face_winding_jump(
    uint64_t face_ind,
    CONSTR_GROUP_T group,
    const std::vector<uint8_t>& flip)
{
    std::vector<uint32_t> covering;
    covering_constraints(face_ind, group, covering);
    if (covering.empty()) return 0;

    const BSPface& f = faces[face_ind];
    const uint32_t cv_id = cell_vertex_off_face(f, f.conn_cells[0]);
    if (cv_id == UINT32_MAX) ip_error("Degenerate cell\n");
    const genericPoint* cv = vertices[cv_id];

    int jump = 0;
    for (uint32_t constr : covering) {
        const uint32_t* cid = constraints_vrts.data() + constr * 3;
        const int o =
            genericPoint::orient3D(*cv, *vertices[cid[0]], *vertices[cid[1]], *vertices[cid[2]]);
        if (o == 0) ip_error("Degenerate cell\n");
        jump += ((o > 0) != (flip[constr] != 0)) ? 1 : -1;
    }
    return jump;
}

// Marks with internal_label the cells inside the surface made of the `group` constraints (the
// faces coloured with skin_colour).
//
// 1. Exact: the winding number is constant on each cell, equal on both sides of a face that no
//    constraint covers and changes by face_winding_jump across one that some do. A breadth-first
//    walk from the outside (the ghost cell beyond the domain boundary, where it is 0) assigns
//    it to every cell using only exact orientation predicates. When the input is closed and
//    consistently oriented the assignment is consistent (the same value is reached whichever
//    way a cell is approached) and a cell is inside iff its winding number is nonzero.
// 2. Fallback: an open or inconsistently oriented surface makes the walk inconsistent (walking
//    around a boundary edge changes the value). Then every cell gets the generalized winding
//    number of the surface at a point inside it (Jacobson et al. 2013; Barill et al. 2018 for
//    the evaluation), and the cell is inside iff that exceeds 1/2 in magnitude. This is what
//    fTetWild does to filter its tetrahedra.
void BSPcomplex::markInternalCells_winding(
    uint32_t skin_colour,
    CONSTR_GROUP_T group,
    uint32_t internal_label,
    bool verbose)
{
    const uint64_t num_cells = cells.size();
    const uint64_t ghost = num_cells; // the outside of the domain

    std::vector<double> coords(vertices.size() * 3);
    for (size_t i = 0; i < vertices.size(); i++)
        vertices[i]->getApproxXYZCoordinates(coords[i * 3], coords[i * 3 + 1], coords[i * 3 + 2]);

    // The group's surface: its input triangles and hole caps, with a consistent orientation.
    std::vector<uint32_t> group_tris;
    for (uint32_t c = 0; c < first_edge_constraint; c++)
        if (constraint_group[c] == group) group_tris.push_back(c);
    std::vector<uint8_t> flip(first_edge_constraint, 0);
    repair_orientation(group_tris, constraints_vrts, coords, flip);
    if (verbose) {
        size_t flipped = 0;
        for (uint32_t c : group_tris) flipped += flip[c];
        if (flipped)
            printf(
                "\tOrientation of %zu of %zu input triangles reversed for the classification\n",
                flipped,
                group_tris.size());
    }

    // Dual graph: one record per face
    struct Arc
    {
        uint64_t c0, c1;
        int jump; // winding number of c1 minus winding number of c0
    };
    std::vector<Arc> arcs;
    arcs.reserve(faces.size());
    std::vector<std::vector<uint32_t>> incident(num_cells + 1); // cell -> arcs (index)
    for (uint64_t fi = 0; fi < faces.size(); fi++) {
        const BSPface& f = faces[fi];
        Arc a;
        a.c0 = f.conn_cells[0];
        a.c1 = (f.conn_cells[1] == UINT64_MAX) ? ghost : f.conn_cells[1];
        a.jump = isSkinFace(f, skin_colour) ? face_winding_jump(fi, group, flip) : 0;
        incident[a.c0].push_back((uint32_t)arcs.size());
        incident[a.c1].push_back((uint32_t)arcs.size());
        arcs.push_back(a);
    }

    // Breadth-first walk from the ghost cell
    const int64_t UNSET = INT64_MIN;
    std::vector<int64_t> w(num_cells + 1, UNSET);
    std::vector<uint64_t> queue;
    queue.reserve(num_cells + 1);
    w[ghost] = 0;
    queue.push_back(ghost);
    for (size_t head = 0; head < queue.size(); head++) {
        const uint64_t c = queue[head];
        for (uint32_t ai : incident[c]) {
            const Arc& a = arcs[ai];
            const uint64_t other = (a.c0 == c) ? a.c1 : a.c0;
            if (w[other] != UNSET) continue;
            w[other] = (a.c0 == c) ? w[c] + a.jump : w[c] - a.jump;
            queue.push_back(other);
        }
    }

    uint64_t unreached = 0, conflicts = 0;
    for (uint64_t c = 0; c <= num_cells; c++)
        if (w[c] == UNSET) unreached++;
    for (const Arc& a : arcs)
        if (w[a.c0] != UNSET && w[a.c1] != UNSET && w[a.c1] - w[a.c0] != a.jump) conflicts++;
    const bool consistent = unreached == 0 && conflicts == 0;
    if (verbose)
        printf(
            "\tWinding number walk: %llu cells unreached, %llu inconsistent faces\n",
            (unsigned long long)unreached,
            (unsigned long long)conflicts);

    if (consistent) {
        if (verbose) printf("\tExact winding number classification\n");
        for (uint64_t c = 0; c < num_cells; c++)
            if (w[c] != 0) setInternalCell(cells[c], internal_label);
        return;
    }

    // Fallback: generalized winding number of the group's triangles at each cell's barycenter.
    if (verbose)
        printf(
            "\tSurface is open or inconsistently oriented: generalized winding number "
            "classification\n");

    std::vector<uint32_t> tris;
    tris.reserve(3 * group_tris.size());
    for (uint32_t c : group_tris) {
        tris.push_back(constraints_vrts[3 * c]);
        tris.push_back(constraints_vrts[3 * c + (flip[c] ? 2 : 1)]);
        tris.push_back(constraints_vrts[3 * c + (flip[c] ? 1 : 2)]);
    }
    const FastWindingNumber fwn(coords.data(), tris.data(), tris.size() / 3);

    std::vector<double> gwn(num_cells, 0.0);
    parallel_blocks(num_cells, [&](uint64_t begin, uint64_t end) {
        std::vector<uint32_t> cell_vrts;
        for (uint64_t c = begin; c < end; c++) {
            cell_vrts.clear();
            list_cellVertices_ts(cells[c], cell_vrts);
            double q[3] = {0.0, 0.0, 0.0};
            for (uint32_t v : cell_vrts)
                for (int k = 0; k < 3; k++) q[k] += coords[3 * v + k];
            const double inv = 1.0 / (double)cell_vrts.size();
            for (int k = 0; k < 3; k++) q[k] *= inv;
            gwn[c] = fwn(q);
        }
    });
    for (uint64_t c = 0; c < num_cells; c++)
        if (std::fabs(gwn[c]) > 0.5) setInternalCell(cells[c], internal_label);
}

// ---------------------------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------------------------
void BSPcomplex::constraintsSurface_complexPartition(
    bool two_files,
    CellClassifier classifier,
    bool verbose)
{
    // Make all cells external
    for (size_t i = 0; i < cells.size(); i++) setExternalCell(cells[i]);

    // Clear vrts_visit for use in cell_vertex_off_face()
    for (size_t i = 0; i < vertices.size(); i++) vrts_visit[i] = 0;

    if (classifier == CellClassifier::GRAPH_CUT) {
#ifdef VOLUMEREMESHER_HAVE_GRAPHCUT
        partition_graphcut(two_files);
#else
        ip_error(
            "The graph-cut classifier is not available: this build was configured without "
            "VOLUMEREMESHER_WITH_GRAPHCUT\n");
#endif
        return;
    }

    markInternalCells_winding(BLACK_A, CONSTR_A, INTERNAL_A, verbose);
    if (two_files) markInternalCells_winding(BLACK_B, CONSTR_B, INTERNAL_B, verbose);
}
} // namespace vol_rem
