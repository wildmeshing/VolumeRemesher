#include "BSP.h"
#include <time.h>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <set>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include "delaunay.h"

#define OPPOSITE_SIGNE(a, b) (a < 0 && b > 0) || (a > 0 && b < 0)
#define MIN_VECT_ELEM(v, n, it, i_min)    \
    it = 1;                               \
    i_min = 0;                            \
    do {                                  \
        if (v[it] < v[i_min]) i_min = it; \
        it++;                             \
    } while (it < n)
#define SAME_EDGE_ENDPTS(e1, e2, E1, E2) ((e1 == E1 && e2 == E2) || (e1 == E2 && e2 == E1))
#define CONSECUTIVE_EDGES(v1, v2, u1, u2) \
    (v1 == u1 || v1 == u2 || v2 == u1 || v2 == u2) // assuming <u1,u2> != <v1,v2>
#define FIND_VECT_POS(e, v, pos) \
    pos = 0;                     \
    do {                         \
        if (v[pos] == e) break;  \
        pos++;                   \
    } while (pos < v.size()) // !! assumes e in v. !!
#define REMOVE_ELEM_VECT(e, v) v.erase(std::find(v.begin(), v.end(), e))
#define IS_GHOST_CELL(c) (c == UINT64_MAX)
#define IS_GHOST_TET(t) (mesh->tet_node[4 * t + 3] == UINT32_MAX)

namespace vol_rem {
//----------------
// General purpose
//----------------

// Input: vector of uint64_t type elements: vect,
//        number of elements to shift: num_shift.
// Output: by using vect returns the original vector shifted-down (i.e.
//         according to increasing order) of num_shift position.
// EX. vect = {3,7,2,5,2} num_shift = 2 -> vect = {5,2,3,7,2}.
inline void UINT64_vect_down_shift(vector<uint64_t>& vect, uint64_t num_shift)
{
    vector<uint64_t> tmp(num_shift, UINT64_MAX);
    uint64_t shift_start_pos = vect.size() - num_shift;
    for (uint64_t pos = 0; pos < num_shift; pos++) tmp[pos] = vect[shift_start_pos + pos];
    for (uint64_t pos = shift_start_pos - 1; pos >= 0; pos--) {
        vect[pos + num_shift] = vect[pos];
        if (pos == 0) break; // Otherwise pos(uint64_t) becomes negative, i.e. error!!
    }
    for (uint64_t pos = 0; pos < num_shift; pos++) vect[pos] = tmp[pos];
}

//  Input: vector of uint64_t type elements: vect,
//        number of elements to shift: num_shift.
// Output: by using vect returns the original vector shifted-up (i.e.
//         according to decreasing order) of num_shift position.
// EX. vect = {3,7,2,5,2} num_shift = 2 -> vect = {3,7,2,5,2}.
inline void UINT64_vect_up_shift(vector<uint64_t>& vect, uint64_t num_shift)
{
    vector<uint64_t> tmp(num_shift, UINT64_MAX);
    uint64_t pos;
    for (pos = 0; pos < vect.size(); pos++) {
        if (pos < num_shift)
            tmp[pos] = vect[pos];
        else
            vect[pos - num_shift] = vect[pos];
    }
    uint64_t refill_pos = vect.size() - num_shift;
    for (pos = 0; pos < num_shift; pos++) vect[refill_pos + pos] = tmp[pos];
}

//  Input: the endpoints of two CONSECUTIVE edges u = <u0,u1>, v=<v0,v1>
// Output: the common endpoint.
inline uint32_t consecEdges_common_endpt(uint32_t u0, uint32_t u1, uint32_t v0, uint32_t v1)
{
    if (u0 == v0 || u0 == v1) return u0;
    if (u1 == v0 || u1 == v1) return u1;
    // If goes here, something wrong.
    printf("\n[BSP.cpp]consecEdges_common_endpt: ERROR no match.\n");
    return UINT32_MAX;
}

//  Input: the endpoints of an edge u = <u0,u1>,
//         the common endpoint between an edge v CONSECUTIVE to u,
//         and u itself: w.
// Output: the endpoint of u different from v_comm.
inline uint32_t other_edge_endpt(uint32_t u0, uint32_t u1, uint32_t w)
{
    if (w == u0) return u1;
    if (w == u1) return u0;
    // If goes here, something wrong.
    printf("\n[BSP.cpp]other_edge_endpt_ind: ERROR no match.\n");
    return UINT32_MAX;
}

//----------------
// BSPedge Methods
//----------------

BSPedge BSPedge::split(uint32_t new_point)
{
    BSPedge e;
    if (meshVertices[2] == UINT32_MAX) {
        e.meshVertices[0] = meshVertices[0];
        e.meshVertices[1] = meshVertices[1];
        e.meshVertices[2] = UINT32_MAX;
        // A mesh-edge keeps UINT32_MAX in positions 2..5 (see BSPedge ctor and
        // the meshVertices invariant). The old code left [3..5] uninitialized,
        // which edges_ShareCommonPlanes() later read -> nondeterministic result.
        e.meshVertices[3] = UINT32_MAX;
        e.meshVertices[4] = UINT32_MAX;
        e.meshVertices[5] = UINT32_MAX;
    } else {
        e.meshVertices[0] = meshVertices[0];
        e.meshVertices[1] = meshVertices[1];
        e.meshVertices[2] = meshVertices[2];
        e.meshVertices[3] = meshVertices[3];
        e.meshVertices[4] = meshVertices[4];
        e.meshVertices[5] = meshVertices[5];
    }
    e.vertices[0] = vertices[0];
    e.vertices[1] = new_point;
    vertices[0] = new_point;
    e.conn_face_0 = conn_face_0;
    return e;
}


//----------------
// BSPface Methods
//----------------
inline void BSPface::exchange_conn_cell(uint64_t cell, uint64_t newCell)
{
#ifdef DEBUG_BSP
    if (conn_cells[0] != cell && conn_cells[1] != cell)
        printf(
            "\n[BSP.h]BSPface::exchange_conn_cell: ERROR no macth for cell "
            "#%llu with this face.\n",
            cell);
#endif

    if (conn_cells[0] == cell)
        conn_cells[0] = newCell;
    else
        conn_cells[1] = newCell;
}

inline void BSPface::removeEdge(uint64_t edge_face_ind)
{
    // Special case: edge is the last element of edges.
    if (edge_face_ind == edges.size() - 1) edges.pop_back();
    // General case: erase it, but later.. at the moment mark it
    else
        edges[edge_face_ind] = UINT64_MAX;
}

//----------------
// BSPcell Methods
//----------------

inline void BSPcell::removeFace(uint64_t face_pos)
{
    if (face_pos != faces.size() - 1) faces[face_pos] = faces[faces.size() - 1];
    faces.pop_back();
}

//--------------------
// BSPcomplex Methods
//--------------------

//  Input: the index of a BSPedge: edge,
//         the index of a BSPface: face.
// Output: nothing.
// Note. BSPedge is attached to BSPface and viceversa.
// Note. this function is used ONLY during the conversion from Delaunay mesh to
//       BSPcomplex.
inline void BSPcomplex::assigne_edge_to_face(uint64_t edge, uint64_t face)
{
    faces[face].edges.push_back(edge);
    // edges[edge].conn_faces.push_back(face);
    edges[edge].conn_face_0 = face;
}

//  Input: the index of 2 cells w.r.t. the vector cells: c1, c2,
// Output: the index of the fece, w.r.t. the vector faces, that lies between
//         the two input-cells.
// Note. It is assumed that both the input-cells are bounded by the input-face.
uint64_t BSPcomplex::faceSharedWithCell(uint64_t c1, uint64_t c2)
{
    BSPcell& cell1 = cells[c1];
    for (uint64_t i = 0; i < cell1.faces.size(); i++)
        if (faces[cell1.faces[i]].conn_cells[0] == c2 || faces[cell1.faces[i]].conn_cells[1] == c2)
            return cell1.faces[i];

#ifdef DEBUG_BSP
    printf(
        "\n[BSP.cpp]BSPcomplex::faceSharedWithCell: ERROR no common face "
        "between cell #%llu and cell #%llu.\n",
        c1,
        c2);
#endif
    return UINT64_MAX; // Never reached.
}

//  Input: a BSPcell: cell.
// Output: returns the number of edges of the BSPcell.
uint64_t BSPcomplex::count_cellEdges(const BSPcell& cell)
{
    uint64_t num_halfedges = 0;
    for (const uint64_t f : cell.faces) num_halfedges += faces[f].edges.size();
    return num_halfedges / 2;
}

//  Input: a BSPcell: cell.
// Output: returns the number of edges of the BSPcell.
uint32_t BSPcomplex::count_cellVertices(const BSPcell& cell, uint64_t* num_cellEdges)
{
    if ((*num_cellEdges) == UINT64_MAX) (*num_cellEdges) = count_cellEdges(cell);
    // Euler formula: num_cellVrts = num_cellEdges + 2 - num_cellFaces
    return (uint32_t)((*num_cellEdges) + 2 - cell.faces.size());
}

//  Input: a BSPcell: cell,
//         vector of type edge index: cell_edges.
// Output: by using cell_edges returns the indices of
//         the edges (w.r.t. vector edges) of the BSPcell.
void BSPcomplex::list_cellEdges(BSPcell& cell, vector<uint64_t>& cell_edges)
{
    uint64_t edge_ind, ce_ind = 0;
    for (uint64_t f = 0; f < cell.faces.size(); f++) {
        BSPface& face = faces[cell.faces[f]];
        for (uint64_t e = 0; e < face.edges.size(); e++) {
            edge_ind = face.edges[e];
            if (edge_visit[edge_ind] == 0) {
                cell_edges[ce_ind++] = edge_ind;
                edge_visit[edge_ind] = 1;
            }
        }
    }

    // Reset edge_visit
    for (uint64_t e = 0; e < cell_edges.size(); e++) edge_visit[cell_edges[e]] = 0;
}

//  Input: a BSPcell: cell,
//         the number of edges of the BSPcell: num_cellEdges,
//         vector of type vertex index: cell_vrts.
// Output: by using cell_vrts returns the indices of
//         the vertices (w.r.t. vector vertices) of the BSPcell.
void BSPcomplex::list_cellVertices(
    BSPcell& cell,
    uint64_t num_cellEdges,
    vector<uint32_t>& cell_vrts)
{
    vector<uint64_t> cell_edges(num_cellEdges, UINT64_MAX);
    list_cellEdges(cell, cell_edges);

    uint32_t v, cv_ind = 0;
    for (uint64_t e = 0; e < cell_edges.size(); e++) {
        BSPedge& edge = edges[cell_edges[e]];
        v = edge.vertices[0];
        if (vrts_visit[v] == 0) {
            cell_vrts[cv_ind++] = v;
            vrts_visit[v] = 1;
        }
        v = edge.vertices[1];
        if (vrts_visit[v] == 0) {
            cell_vrts[cv_ind++] = v;
            vrts_visit[v] = 1;
        }
    }

    // Reset vrts_visit
    for (uint32_t u = 0; u < cell_vrts.size(); u++) vrts_visit[cell_vrts[u]] = 0;
}

//  Input: a BSPface: face,
//         vector of type vertex index: face_vrts.
// Output: by using face_vrts returns the indices of
//         the vertices (w.r.t. vector vertices) of the BSPface.
void BSPcomplex::list_faceVertices(BSPface& face, vector<uint32_t>& face_vrts)
{
    uint32_t fv_ind = 0;

    // Add both endpoints of first edge.
    BSPedge& edge0 = edges[face.edges[0]];
    uint32_t e0 = edge0.vertices[0];
    uint32_t e1 = edge0.vertices[1];
    face_vrts[fv_ind++] = e0;
    face_vrts[fv_ind++] = e1;

    // Find the common endpoint between first and second edge.
    BSPedge& edge1 = edges[face.edges[1]];
    uint32_t link_vrt = consecEdges_common_endpt(e0, e1, edge1.vertices[0], edge1.vertices[1]);
    if (link_vrt == e0) {
        face_vrts[0] = e1;
        face_vrts[1] = e0;
    }

    // Walk the boundary and add the endpoint != link_vrt.
    for (uint64_t e = 1; e < face.edges.size() - 1; e++) {
        BSPedge& edge = edges[face.edges[e]];

        if (link_vrt == edge.vertices[0])
            link_vrt = edge.vertices[1];
        else
            link_vrt = edge.vertices[0];

        face_vrts[fv_ind++] = link_vrt;
    }
}

//  Input: a BSPcell: cell,
//         vector of edges indices type: cell_edges,
//         vector of vertices indices type: cell_vrts.
// Output: by using cell_edges returns the indices (w.r.t. vector edges)
//         of cell's edges,
//         by using cell_vrts returns the indices (w.r.t. vector vertices)
//         of cell's vertices,
void BSPcomplex::fill_cell_locDS(
    BSPcell& cell,
    vector<uint64_t>& cell_edges,
    vector<uint32_t>& cell_vrts)
{
    uint64_t edge_ind, ce_ind = 0, cv_ind = 0;
    uint32_t e0, e1;
    for (uint64_t f = 0; f < cell.faces.size(); f++) {
        BSPface& face = faces[cell.faces[f]];
        for (uint64_t e = 0; e < face.edges.size(); e++) {
            edge_ind = face.edges[e];
            BSPedge& edge = edges[edge_ind];

            if (edge_visit[edge_ind] == 0) {
                // Fill cell_edges.
                cell_edges[ce_ind++] = edge_ind;
                edge_visit[edge_ind] = 1;

                e0 = edge.vertices[0];
                e1 = edge.vertices[1];

                // Fill cell_vrts.
                if (vrts_visit[e0] == 0) {
                    cell_vrts[cv_ind++] = e0;
                    vrts_visit[e0] = 1;
                }
                if (vrts_visit[e1] == 0) {
                    cell_vrts[cv_ind++] = e1;
                    vrts_visit[e1] = 1;
                }
            }
        }
    }

    // Reset edge_visit and vrts_visit.
    for (uint64_t e = 0; e < cell_edges.size(); e++) edge_visit[cell_edges[e]] = 0;
    for (uint64_t v = 0; v < cell_vrts.size(); v++) vrts_visit[cell_vrts[v]] = 0;
}

//  Input: a BSPface: face,
//         indices of two edge endpoints (w.r.t. vector vertices): u, v,
// Output: if exists returns the index of a BSPedge with endpoints u and v
//         that belong to faces[face], otherwise return UINT64_MAX.
inline uint64_t BSPcomplex::find_face_edge(const BSPface& face, uint32_t v, uint32_t u)
{
#ifdef DEBUG_BSP
    if (face.edges.size() == 0) printf("\n[BSP.cpp]find_face_edge: ERROR face have no edges.\n");
#endif

    uint64_t edge_ind;
    for (uint64_t e = 0; e < face.edges.size(); e++) {
        edge_ind = face.edges[e];
        BSPedge& edge = edges[edge_ind];
        if (SAME_EDGE_ENDPTS(u, v, edge.vertices[0], edge.vertices[1])) return edge_ind;
    }

#ifdef DEBUG_BSP
    printf(
        "[BSP.cpp]BSPcomplex::find_face_edge: "
        "WARNING no match for edge <%u,%u> with this face\n",
        u,
        v);
#endif

    return UINT64_MAX; // Never reached if <u,v> belong to the BSPcell.
}

//
//
uint64_t BSPcomplex::count_cellFaces_inc_cellVrt(const BSPcell& cell, uint32_t v)
{
    uint64_t k = 0;
    for (const uint64_t fi : cell.faces)
        for (const uint64_t ei : faces[fi].edges)
            if (edges[ei].vertices[0] == v || edges[ei].vertices[1] == v) k++;
    return k / 2;
}

//
//
void BSPcomplex::cell_VFrelation(const BSPcell& cell, uint32_t v, vector<uint64_t>& v_incFaces_ind)
{
    uint64_t k = 0;
    for (const uint64_t fi : cell.faces)
        for (const uint64_t ei : faces[fi].edges)
            if (edges[ei].vertices[0] == v || edges[ei].vertices[1] == v) {
                v_incFaces_ind[k++] = fi;
                break;
            }
}

//
//
void BSPcomplex::COMPL_cell_VFrelation(
    const BSPcell& cell,
    uint32_t v,
    vector<uint64_t>& v_NOT_incFaces_ind)
{
    uint64_t k = 0;
    for (const uint64_t fi : cell.faces) {
        bool has_edge = false;
        for (const uint64_t ei : faces[fi].edges)
            if (edges[ei].vertices[0] == v || edges[ei].vertices[1] == v) {
                has_edge = true;
                break;
            }
        if (!has_edge) v_NOT_incFaces_ind[k++] = fi;
    }
}

//
//
bool BSPcomplex::is_virtual(uint32_t constr_ind)
{
    return (constr_ind >= first_virtual_constraint);
}

// -- Geometric Predicates ---------------------------------------------------

// Returns TRUE if v either is a vertex of the plane p0, p1, p2, or it was built
// by intersecting such a plane with other simplexes.
bool isVertexBuiltFromPlane(
    const genericPoint* v,
    const explicitPoint3D* p0,
    const explicitPoint3D* p1,
    const explicitPoint3D* p2)
{
    if (v->isExplicit3D()) {
        return ((v == p0) || (v == p1) || (v == p2));
    } else if (v->isLPI()) {
        const implicitPoint3D_LPI& lpi = v->toLPI();
        return (
            (((p0 == &lpi.P()) && (p1 == &lpi.Q())) || ((p1 == &lpi.P()) && (p0 == &lpi.Q()))) ||
            (((p1 == &lpi.P()) && (p2 == &lpi.Q())) || ((p2 == &lpi.P()) && (p1 == &lpi.Q()))) ||
            (((p2 == &lpi.P()) && (p0 == &lpi.Q())) || ((p0 == &lpi.P()) && (p2 == &lpi.Q()))) ||
            ((p0 == &lpi.R()) && (p1 == &lpi.S()) && (p2 == &lpi.T())));
    } else {
        // TPI
        const implicitPoint3D_TPI& tpi = v->toTPI();
        return (
            ((p0 == &tpi.V1()) && (p1 == &tpi.V2()) && (p2 == &tpi.V3())) ||
            ((p0 == &tpi.W1()) && (p1 == &tpi.W2()) && (p2 == &tpi.W3())) ||
            ((p0 == &tpi.U1()) && (p1 == &tpi.U2()) && (p2 == &tpi.U3())));
    }
}


//  Input: vector of vertices indices (w.r.t. vector vertices): vrts_inds,
//         3 indices of vertices that define a plane: plane_pt0, plane_pt1,
//         plane_pt2.
// Output: by using the global vector vrts_orBin returns the orientations of
//         each point of vrts_inds w.r.t. the plane for
//         {plane_pt0, plane_pt1, plane_pt2}.
void BSPcomplex::vrts_orient_wrtPlane(
    const vector<uint32_t>& vrts_inds,
    uint32_t plane_pt0,
    uint32_t plane_pt1,
    uint32_t plane_pt2,
    uint32_t count)
{
    const explicitPoint3D& p0 = vertices[plane_pt0]->toExplicit3D();
    const explicitPoint3D& p1 = vertices[plane_pt1]->toExplicit3D();
    const explicitPoint3D& p2 = vertices[plane_pt2]->toExplicit3D();

    for (uint32_t v = 0; v < vrts_inds.size(); v++) {
        genericPoint* vrt = vertices[vrts_inds[v]];
        if (isVertexBuiltFromPlane(vrt, &p0, &p1, &p2))
            vrts_orBin[vrts_inds[v]] = 0;
        else {
            vrts_orBin[vrts_inds[v]] = genericPoint::orient3D(*vrt, p0, p2, p1);
        }
    }
}


//  Input:
// Output:
inline void BSPcomplex::count_vrt_orBin(
    const vector<uint32_t>& inds,
    uint32_t* pos,
    uint32_t* neg,
    uint32_t* zero)
{
    (*pos) = 0;
    (*neg) = 0;
    (*zero) = 0;

    for (uint32_t i = 0; i < inds.size(); i++)
        if (vrts_orBin[inds[i]] == 0)
            (*zero)++;
        else if (vrts_orBin[inds[i]] == 1)
            (*pos)++;
        else if (vrts_orBin[inds[i]] == -1)
            (*neg)++;
#ifdef DEBUG_BSP
        else
            printf(
                "[BSP.cpp]BSPcomplex::count_vrt_orBin: ERROR "
                "wrong access to vrt_orBin[%u]\n",
                inds[i]);
#endif
}

//  Input: a BSPedge: edge,
//         vector of the inices of the vertices of the BSPcell to which the
//         BSPedge belong to: cell_vrts.
// Output: true if the constraint intersects the edge interior,
//         false otherwise.
inline bool BSPcomplex::constraint_innerIntersects_edge(
    const BSPedge& e,
    const vector<uint32_t>& cell_vrts)
{
    return OPPOSITE_SIGNE(vrts_orBin[e.vertices[0]], vrts_orBin[e.vertices[1]]);
}

//  Input: vector of the inices of the vertices of the BSPface: face_vrts.
// Output: true if the constraint intersects the face interior,
//         false otherwise.
inline bool BSPcomplex::constraint_innerIntersects_face(const vector<uint32_t>& face_vrts)
{
    // Face vertices disposition w.r.t. constraint-plane.
    uint32_t vrtsOVER, vrtsUNDER, vrtsON;
    count_vrt_orBin(face_vrts, &vrtsOVER, &vrtsUNDER, &vrtsON);

    // Check if faces[face_ind] has to be splitted
    // (at least 2 vertices with non-zero opposite cell_vrts_orient)
    return (vrtsUNDER > 0 && vrtsOVER > 0);
}

//
//
bool BSPcomplex::coplanar_constraint_innerIntersects_face(
    const std::vector<uint64_t>& fedges,
    const uint32_t tri[3],
    int xyz)
{
    // ASSUMPTION: No face vertices are in the interior of the constraint

    // Process: intersection must be bound by at least three unaligned points on the boundary of tri
    // These points must be on either two different tri_edges,
    // or one of a vertex and two on the opposite edge, or on three vertices

    int mask = 0;

    const BSPedge& edge0 = edges[fedges.back()];
    const BSPedge& edge1 = edges[fedges[0]];
    uint32_t vid_0 = consecEdges_common_endpt(
        edge0.vertices[0],
        edge0.vertices[1],
        edge1.vertices[0],
        edge1.vertices[1]);

    // 1) Cerca i tre vertici sul bordo della faccia
    // Per ogni tri_vertex
    //  cerca vertici coincidenti in faccia: se trovi attiva campo in mask e passa al vertice
    //  successivo se campo in mask non è gia attivo, cerca edge della faccia che lo contengono:
    //     se trovi attiva campo in mask e passa al vertice successivo
    // Se mask == 7 torna true
    for (int i = 0; i < 3; i++) {
        uint32_t vid = vid_0;
        for (uint64_t e = 0; e < fedges.size(); e++) {
            const BSPedge& edge = edges[fedges[e]];
            if (vid == edge.vertices[0])
                vid = edge.vertices[1];
            else
                vid = edge.vertices[0];
            if (tri[i] == vid) {
                mask |= (1 << i);
                break;
            }
        }
    }
    if (mask == 7) return true;

    for (int i = 0; i < 3; i++)
        if (!(mask & (1 << i))) {
            for (uint64_t e = 0; e < fedges.size(); e++) {
                const BSPedge& edge = edges[fedges[e]];
                const genericPoint* ev1 = vertices[edge.vertices[0]];
                const genericPoint* ev2 = vertices[edge.vertices[1]];
                if (genericPoint::pointInInnerSegment(*vertices[tri[i]], *ev1, *ev2, xyz)) {
                    mask |= (1 << i);
                    break;
                }
            }
        }
    if (mask == 7) return true;

    // 2) Cerca vertici della faccia su edge del constraint
    //  Per ogni tri_edge
    //   cerca vertici in faccia che stanno innerEdge: se trovi attiva campi in mask e edge
    //   successivo se campi in mask non sono gia attivi, cerca edge della faccia che lo
    //   innerCrossano: se trovi attiva campi in mask e passa a edge successivo
    // Se mask == 7 torna true altrimenti false
    for (int i = 0; i < 3; i++) {
        uint32_t ti0 = (i + 1) % 3;
        uint32_t ti1 = (i + 2) % 3;
        if ((mask & (1 << ti0)) && (mask & (1 << ti1))) continue;
        uint32_t vid = vid_0;
        for (uint64_t e = 0; e < fedges.size(); e++) {
            const BSPedge& edge = edges[fedges[e]];
            if (vid == edge.vertices[0])
                vid = edge.vertices[1];
            else
                vid = edge.vertices[0];
            const genericPoint* ev1 = vertices[tri[ti0]];
            const genericPoint* ev2 = vertices[tri[ti1]];
            if (genericPoint::pointInInnerSegment(*vertices[vid], *ev1, *ev2, xyz)) {
                mask |= (1 << ti0);
                mask |= (1 << ti1);
                break;
            }
        }
    }
    if (mask == 7) return true;

    for (int i = 0; i < 3; i++) {
        uint32_t ti0 = (i + 1) % 3;
        uint32_t ti1 = (i + 2) % 3;
        if ((mask & (1 << ti0)) && (mask & (1 << ti1))) continue;
        for (uint64_t e = 0; e < fedges.size(); e++) {
            const BSPedge& edge = edges[fedges[e]];
            const genericPoint* ev1 = vertices[edge.vertices[0]];
            const genericPoint* ev2 = vertices[edge.vertices[1]];
            const genericPoint* fv1 = vertices[tri[ti0]];
            const genericPoint* fv2 = vertices[tri[ti1]];
            if (genericPoint::innerSegmentsCross(*ev1, *ev2, *fv1, *fv2, xyz)) return true;
        }
    }

    return false;
}

// Returns 1 if point is on the boundary of the triangle, 2 if it is in the interior, 0 otherwise
int localizedPointInTriangle(
    const genericPoint& P,
    const genericPoint& A,
    const genericPoint& B,
    const genericPoint& C,
    int xyz)
{
    int o1, o2, o3;
    if (xyz == 2) {
        o1 = genericPoint::orient2Dxy(P, A, B);
        o2 = genericPoint::orient2Dxy(P, B, C);
        o3 = genericPoint::orient2Dxy(P, C, A);
    } else if (xyz == 0) {
        o1 = genericPoint::orient2Dyz(P, A, B);
        o2 = genericPoint::orient2Dyz(P, B, C);
        o3 = genericPoint::orient2Dyz(P, C, A);
    } else {
        o1 = genericPoint::orient2Dzx(P, A, B);
        o2 = genericPoint::orient2Dzx(P, B, C);
        o3 = genericPoint::orient2Dzx(P, C, A);
    }
    return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0)) +
           ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
}

//-Upload Delaunay triangolation----------------------

//  Input: pointer to mesh,
//         vector of tetrahedra index type: new_order.
// Output: by using new_order returns new cells indexing (the same as Delaunay
//         tetrahedra indexing, but without ghost-tets).
uint64_t BSPcomplex::removing_ghost_tets(const TetMesh* mesh, vector<uint64_t>& new_order)
{
    // new_order have as many element as mesh->tet_num,
    // new_order[i-th tet] =
    //    i - (num of ghost_tet between 0 and i) IF i-th tet is non-ghost
    //    UINT64_MAX                             IF i-th tet is ghost
    uint64_t ghost_tet_count = 0;
    for (uint64_t tet_ind = 0; tet_ind < mesh->tet_num; tet_ind++) {
        if (IS_GHOST_TET(tet_ind)) {
            new_order[tet_ind] = UINT64_MAX;
            ghost_tet_count++;
        } else
            new_order[tet_ind] = tet_ind - ghost_tet_count;
    }
    return mesh->tet_num - ghost_tet_count;
}

//  Input: pointer to mesh,
//         the 2 endpoints of a tetrahedron edge: e0, e1,
//         the index of the tetrahedron (tet) to which edge belongs: tet_ind,
//         new cells indexing (the same as tetrahedra indexing, but without
//         ghost-tets): new_order.
// Output: the index of the edge (w.r.t. the vector edges) of the
//         edge <endpt0, endpt1>.
// Note. tetrahedra are added in crescent index order.
uint64_t BSPcomplex::add_tetEdge(
    const TetMesh* mesh,
    uint32_t e0,
    uint32_t e1,
    uint64_t tet_ind,
    const vector<uint64_t>& new_order)
{
    // Tetrahedra incident in <endpt0,endpt1>
    uint32_t edge_ends[2] = {e0, e1};
    uint64_t num_incTet;
    uint64_t* incTet = mesh->ETrelation(edge_ends, tet_ind, &num_incTet);
    // Note. ETrelation can return ghost-tet.

    uint64_t min = UINT64_MAX;
    for (uint32_t i = 0; i < num_incTet; i++)
        if (!IS_GHOST_TET(incTet[i]) && incTet[i] < min) min = incTet[i];

    free(incTet);

    if (min == tet_ind) {
        // None of the tetrahedra incident in <e0,e1> has been visited,
        // furthermore the edge <e0,e1> do not belongs to a face of tet_ind
        // (as a consequence of the constructor BSPcomplex).
        // A new BSPedge has to be created.
        edges.push_back(BSPedge(e0, e1, e0, e1));
        return edges.size() - 1;
    }

    uint64_t edge_ind;
    for (uint64_t f = 0; f < 4; f++) { // Note. also with f<3 sholud work.
        BSPface& face = faces[cells[new_order[min]].faces[f]];
        for (uint64_t e = 0; e < 3; e++) {
            edge_ind = face.edges[e];
            BSPedge& edge = edges[edge_ind];
            if (SAME_EDGE_ENDPTS(e0, e1, edge.vertices[0], edge.vertices[1])) break;
        }
        BSPedge& edge = edges[edge_ind];
        if (SAME_EDGE_ENDPTS(e0, e1, edge.vertices[0], edge.vertices[1])) break;
    }
    return edge_ind;
}

//  Input: three face (triangle) vertices: v0, v1, v2,
//         index of a BSPcell (w.r.t. vector cells) owning the face: cell_ind,
//         index of a BSPcell (w.r.t. vector cells) faced at the previous one
//         throught the face, if it does not exists (i.e. the previous cell is
//         a complex boundary cell) then it is UINT64_MAX: adjCell_ind.
// Output: returns the index (w.r.t. the vector faces) of the new BSPface.
inline uint64_t BSPcomplex::add_tetFace(
    uint32_t v0,
    uint32_t v1,
    uint32_t v2,
    uint64_t cell_ind,
    uint64_t adjCell_ind)
{
    // Note. adjCell==UINT64_MAX -> convex-hull face.
    faces.push_back(BSPface(v0, v1, v2, cell_ind, adjCell_ind));
    uint64_t face_ind = faces.size() - 1;
    cells[cell_ind].faces.push_back(face_ind);
    if (!IS_GHOST_CELL(adjCell_ind)) cells[adjCell_ind].faces.push_back(face_ind);

    return face_ind;
}

//  Input: the index of a Delaunay mesh tetrahedron: tet_ind,
//         the index of a Delaunay mesh tetrahedron adjacent to tet: adjTet_ind,
//         the index of the BSPcell corresponding to the adjacent tetrahedron,
//         (in the case it is a ghost-tet it is UINT64_MAX): adjCell_ind.
// Output: returns true if the face between the two input tetrahedra has not
//         been turned into a BSPface yet, false otherwise.
inline bool BSPcomplex::tet_face_isNew(uint64_t tet_ind, uint64_t adjTet_ind, uint64_t adjCell_ind)
{
    // The face is the common one between tetrahedra indexed as tet and adjTet.
    // Assuming that the tetrahedron indexed as tet has not been visited yet.
    // Check if:
    // - adjTet has not been visited yet -> new face.
    // - adjTet has been already visited -> face already exists.(not new)
    // - adjTet is ghost (i.e. adjCell is UINT64_MAX by new_order) -> new face.
    return (adjTet_ind > tet_ind || IS_GHOST_CELL(adjCell_ind));
}

//
//
inline void BSPcomplex::fill_face_colour(
    uint64_t tet_ind,
    uint64_t face_ind,
    const uint32_t** map_fi,
    const uint32_t* num_map_fi)
{
    if (num_map_fi[tet_ind] == 0)
        faces[face_ind].colour = WHITE;
    else {
        // Count non-virtual constraints.
        uint32_t n = 0;
        for (uint32_t cc = 0; cc < num_map_fi[tet_ind]; cc++)
            if (!is_virtual(map_fi[tet_ind][cc])) n++;

        if (n == 0)
            faces[face_ind].colour = WHITE;
        else {
            faces[face_ind].colour = GREY;
            faces[face_ind].coplanar_constraints.resize(n);
            uint32_t pos = 0;
            for (uint32_t cc = 0; cc < num_map_fi[tet_ind]; cc++)
                if (!is_virtual(map_fi[tet_ind][cc]))
                    faces[face_ind].coplanar_constraints[pos++] = map_fi[tet_ind][cc];
        }
    }
}


// bool faceHasCorrectOrientation(BSPcomplex* cpx, uint64_t f_id)
//{
//     const BSPface& f = cpx->faces[f_id];
//     const uint64_t c_id = f.conn_cells[0];
//     BSPcell& c = cpx->cells[c_id];
//     const uint64_t e0_id = f.edges[0];
//     const uint64_t e1_id = f.edges[1];
//     const uint64_t e2_id = f.edges[2];
//     const BSPedge& e0 = cpx->edges[e0_id];
//     const BSPedge& e1 = cpx->edges[e1_id];
//     const BSPedge& e2 = cpx->edges[e2_id];
//     const uint32_t v0_id = cpx->getFaceVertex(f, 0);
//     const uint32_t v1_id = cpx->getFaceVertex(f, 1);
//     genericPoint* v0 = cpx->vertices[v0_id];
//     genericPoint* v1 = cpx->vertices[v1_id];
//     genericPoint* v2;
//
//     size_t i;
//     for (i = 2; i < f.edges.size(); i++)
//     {
//         const uint32_t v2_id = cpx->getFaceVertex(f, i);
//         v2 = cpx->vertices[v2_id];
//         if (genericPoint::misaligned(*v0, *v1, *v2)) break;
//     }
//     if (i == f.edges.size()) ip_error("Degenerate face\n");
//
//     uint64_t num_cellEdges = UINT64_MAX;
//     uint32_t num_cellVrts = cpx->count_cellVertices(c, &num_cellEdges);
//     vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
//     cpx->list_cellVertices(c, num_cellEdges, cell_vrts);
//     for (uint32_t vi : cell_vrts) if (!cpx->faceHasVertex(f, vi))
//     {
//         genericPoint* ov = cpx->vertices[vi];
//
//         int ori = genericPoint::orient3D(*ov, *v0, *v1, *v2);
//         if (ori == 0) continue;
//         return (ori < 0);
//     }
//     ip_error("Degenerate cell\n");
// }
//
//// Returns the index of the v_ind'th vertex in f.
//// Vertex 'i' is the common vertex between edge 'i' and edge 'i+1'%num_edges
// uint32_t BSPcomplex::getFaceVertex(const BSPface& f, uint32_t v_ind)
//{
//     const BSPedge& e0 = edges[f.edges[v_ind]];
//     const BSPedge& e1 = edges[f.edges[(v_ind + 1) % (f.edges.size())]];
//     return consecEdges_common_endpt(e0.vertices[0], e0.vertices[1], e1.vertices[0],
//     e1.vertices[1]);
// }
//
// bool BSPcomplex::faceHasVertex(const BSPface& f, uint32_t v_ind)
//{
//     for (uint64_t eid : f.edges)
//     {
//         const BSPedge& e = edges[eid];
//         if (e.vertices[0] == v_ind || e.vertices[1] == v_ind) return true;
//     }
//     return false;
// }

//
// Fills the data scruture with the information of the Delauany mesh.
BSPcomplex::BSPcomplex(
    const TetMesh* mesh,
    const constraints_t* _constraints,
    const uint32_t** map,
    const uint32_t* num_map,
    const uint32_t** map_f0,
    const uint32_t* num_map_f0,
    const uint32_t** map_f1,
    const uint32_t* num_map_f1,
    const uint32_t** map_f2,
    const uint32_t* num_map_f2,
    const uint32_t** map_f3,
    const uint32_t* num_map_f3)
{
    // Uploading the vertices of the mesh
    vertices.resize(mesh->num_vertices);
    for (uint32_t vrt = 0; vrt < mesh->num_vertices; vrt++)
        vertices[vrt] = new explicitPoint3D(
            mesh->vertices[vrt].coord[0],
            mesh->vertices[vrt].coord[1],
            mesh->vertices[vrt].coord[2]);

    // Initialize vrts_orBin:
    // since orient3D can be -1, 0 or 1 all elements are set to 2.
    vrts_orBin.resize(mesh->num_vertices, 2);

    // Uploading the constraints (the last num_virtual_triangles constraints are virtual.)
    first_virtual_constraint = _constraints->num_triangles - _constraints->num_virtual_triangles;
    // Carve the real-constraint range [0, first_virtual_constraint) into surface / cap /
    // edge / point sub-ranges (see the BSP.h layout comment). Point-forcing triangles sit
    // just before the virtual ones, edge-forcing just before those; when no holes were
    // filled (num_input_triangles unset) every real constraint is input surface.
    num_edge_triangles = _constraints->num_edge_triangles;
    num_point_triangles = _constraints->num_point_triangles;
    first_point_constraint = first_virtual_constraint - num_point_triangles;
    first_edge_constraint = first_point_constraint - num_edge_triangles;
    first_fake_constraint = (_constraints->num_input_triangles == UINT32_MAX)
        ? first_virtual_constraint
        : _constraints->num_input_triangles;
    constraints_vrts.resize(3 * _constraints->num_triangles);
    constraint_group.resize(_constraints->num_triangles);
    for (uint32_t i = 0; i < _constraints->num_triangles; i++) {
        constraints_vrts[3 * i] = _constraints->tri_vertices[3 * i];
        constraints_vrts[3 * i + 1] = _constraints->tri_vertices[3 * i + 1];
        constraints_vrts[3 * i + 2] = _constraints->tri_vertices[3 * i + 2];
        constraint_group[i] = _constraints->constr_group[i];
    }
    // Copy the input-file triangle index of each real (non-virtual) constraint.
    if (_constraints->tri_original_index) {
        constraint_original_index.resize(first_virtual_constraint);
        for (uint32_t i = 0; i < first_virtual_constraint; i++)
            constraint_original_index[i] = _constraints->tri_original_index[i];
    }

    // Establish new tetrahedtra-(cell) indexing: only non-ghost cell are indexed.
    vector<uint64_t> new_order(mesh->tet_num, UINT64_MAX);
    uint64_t cell_num = removing_ghost_tets(mesh, new_order);

    // Creating as many empty cells as the number of non-ghost tet_
    cells.resize(cell_num);
    edges.reserve(cell_num + mesh->num_vertices);
    faces.reserve(cell_num * 2);


    // Loading the cells, creating faces and eadges of the BSP:
    // cells -> the non-ghost tetrahedra in the mesh,
    // faces -> the faces of the non-ghost tetrahedra in the mesh,
    // edges -> the edges of the non-ghost tetrahedra in the mesh.
    for (uint64_t tet_ind = 0; tet_ind < mesh->tet_num; tet_ind++) {
        // Here each BSPcell is a non-ghost tetrahedron of the mesh:
        // consider a tetrahedron (tet) whose index is tet_ind.
        uint64_t cell_ind = new_order[tet_ind];
        if (IS_GHOST_CELL(cell_ind)) continue; // Skip ghost-tet.

        // Create BSPcells from tetrahedra by following increasing indexing:
        // all non-ghost tetrahedra which have index lower than tet_ind
        // have been already turned into BSP cells.

        // Constraints improperly intersecated by tet.
        if (num_map[tet_ind] > 0) {
            cells[cell_ind].constraints.resize(num_map[tet_ind]);
            for (uint32_t i = 0; i < num_map[tet_ind]; i++)
                cells[cell_ind].constraints[i] = map[tet_ind][i];
        }

        // Adding BSPface and BSPedges to create a BSPcell conformed to tetrahedron.
        uint32_t v[4]; // Indices of tet vertices.
        v[0] = mesh->tet_node[4 * tet_ind];
        v[1] = mesh->tet_node[4 * tet_ind + 1];
        v[2] = mesh->tet_node[4 * tet_ind + 2];
        v[3] = mesh->tet_node[4 * tet_ind + 3];

        uint64_t face_ind, adjCell_ind, adjTet_ind;
        uint64_t tet_edge[6];
        // --- face <v0,v1,v2> -----------------------------
        adjTet_ind = mesh->tet_neigh[4 * tet_ind + 3] >> 2;
        adjCell_ind = new_order[adjTet_ind];
        if (tet_face_isNew(tet_ind, adjTet_ind, adjCell_ind)) {
            face_ind = add_tetFace(v[0], v[1], v[2], cell_ind, adjCell_ind);
            // At most three edges may have to be created <v0,v1>, <v1,v2>, <v2,v0>.
            tet_edge[0] = add_tetEdge(mesh, v[0], v[1], tet_ind, new_order);
            assigne_edge_to_face(tet_edge[0], face_ind);
            tet_edge[2] = add_tetEdge(mesh, v[2], v[0], tet_ind, new_order);
            assigne_edge_to_face(tet_edge[2], face_ind);
            tet_edge[1] = add_tetEdge(mesh, v[1], v[2], tet_ind, new_order);
            assigne_edge_to_face(tet_edge[1], face_ind);
            // Color and coplanar-constraints
            fill_face_colour(tet_ind, face_ind, map_f3, num_map_f3);
        } else {
            face_ind = faceSharedWithCell(cell_ind, adjCell_ind);
            BSPface& face = faces[face_ind];
            tet_edge[0] = find_face_edge(face, v[0], v[1]);
            tet_edge[1] = find_face_edge(face, v[1], v[2]);
            tet_edge[2] = find_face_edge(face, v[2], v[0]);
        }
        // --- face <v3,v0,v1> -----------------------------
        adjTet_ind = mesh->tet_neigh[4 * tet_ind + 2] >> 2;
        adjCell_ind = new_order[adjTet_ind];
        if (tet_face_isNew(tet_ind, adjTet_ind, adjCell_ind)) {
            face_ind = add_tetFace(v[3], v[0], v[1], cell_ind, adjCell_ind);
            // At most two edges may have to be created <v3,v0>, <v1,v3>.
            tet_edge[4] = add_tetEdge(mesh, v[1], v[3], tet_ind, new_order);
            assigne_edge_to_face(tet_edge[4], face_ind);
            tet_edge[3] = add_tetEdge(mesh, v[3], v[0], tet_ind, new_order);
            assigne_edge_to_face(tet_edge[3], face_ind);
            // <v0,v1> is tet_edge[0].
            assigne_edge_to_face(tet_edge[0], face_ind);
            // Color and coplanar-constraints
            fill_face_colour(tet_ind, face_ind, map_f2, num_map_f2);
        } else {
            face_ind = faceSharedWithCell(cell_ind, adjCell_ind);
            BSPface& face = faces[face_ind];
            tet_edge[3] = find_face_edge(face, v[3], v[0]);
            tet_edge[4] = find_face_edge(face, v[1], v[3]);
        }
        // --- face <v2,v3,v0> -----------------------------
        adjTet_ind = mesh->tet_neigh[4 * tet_ind + 1] >> 2;
        adjCell_ind = new_order[adjTet_ind];
        if (tet_face_isNew(tet_ind, adjTet_ind, adjCell_ind)) {
            face_ind = add_tetFace(v[2], v[3], v[0], cell_ind, adjCell_ind);
            // At most one edges may have to be created <v2,v3>.
            tet_edge[5] = add_tetEdge(mesh, v[2], v[3], tet_ind, new_order);
            assigne_edge_to_face(tet_edge[5], face_ind);
            // <v3,v0> is tet_edge[3], <v0,v2> is tet_edge[2].
            assigne_edge_to_face(tet_edge[2], face_ind);
            assigne_edge_to_face(tet_edge[3], face_ind);
            // Color and coplanar-constraints
            fill_face_colour(tet_ind, face_ind, map_f1, num_map_f1);
        } else {
            face_ind = faceSharedWithCell(cell_ind, adjCell_ind);
            tet_edge[5] = find_face_edge(faces[face_ind], v[2], v[3]);
        }
        // --- face <v1,v2,v3> -----------------------------
        adjTet_ind = mesh->tet_neigh[4 * tet_ind] >> 2;
        adjCell_ind = new_order[adjTet_ind];
        if (tet_face_isNew(tet_ind, adjTet_ind, adjCell_ind)) {
            face_ind = add_tetFace(v[1], v[2], v[3], cell_ind, adjCell_ind);
            // No edges have to be created.
            // <v1,v2> is tet_edge[1], <v2,v3> is tet_edge[5], <v3,v1> is tet_edge[4].
            assigne_edge_to_face(tet_edge[1], face_ind);
            assigne_edge_to_face(tet_edge[5], face_ind);
            assigne_edge_to_face(tet_edge[4], face_ind);
            // Color and coplanar-constraints
            fill_face_colour(tet_ind, face_ind, map_f0, num_map_f0);
        }
    }

    // Initialize visit-flag vectors: all the values are set to upper-limit.
    vrts_visit.resize(vertices.size(), 0);
    edge_visit.resize(edges.size(), 0);


    // Verify that conn_cell[0] is below the face for every face
    // for (size_t fid = 0; fid < faces.size(); fid++)
    //    if (!faceHasCorrectOrientation(this, fid))
    //        ip_error("Wrong orientation\n");
}

//-BSP subdivision----------------------

//  Input: index of a BSPedge (w.r.t. vector face.edges): edge_face_ind,
//         index of two BSPfaces (w.r.t. vector faces): face_ind, newFace_ind.
// Output: nothing.
// Note. edge is removed from faces[face_ind] and assigned to faces[newFace_ind].
inline void BSPcomplex::move_edge(uint64_t edge_face_ind, uint64_t face_ind, uint64_t newFace_ind)
{
    uint64_t edge_ind = faces[face_ind].edges[edge_face_ind];
    edges[edge_ind].conn_face_0 = newFace_ind;
    faces[newFace_ind].edges.push_back(edge_ind);
    faces[face_ind].removeEdge(edge_face_ind);
}

//  Input: index of a BSPface (w.r.t. vector cells[cell].faces): face_cell_ind,
//         index of two BSPcells (w.r.t. vector cells): cell_ind, newCell_ind.
// Output: nothing.
// Note. face is removed from cells[cell_ind] and assigned to cells[newCell_ind].
inline void BSPcomplex::move_face(uint64_t face_cell_ind, uint64_t cell_ind, uint64_t newCell_ind)
{
    uint64_t face_ind = cells[cell_ind].faces[face_cell_ind];
    faces[face_ind].exchange_conn_cell(cell_ind, newCell_ind);
    cells[newCell_ind].faces.push_back(face_ind);
    cells[cell_ind].removeFace(face_cell_ind);
}

//  Input: index of a constraint (w.r.t. vector cells[cell].constraints):
//         constr_cell_ind,
//         index of a BSPcells (w.r.t. vector cells): cell_ind.
// Output: nothing.
inline void BSPcomplex::remove_constraint(uint32_t constr_cell_ind, uint64_t cell_ind)
{
    BSPcell& cell = cells[cell_ind];
    size_t last_ind = cell.constraints.size() - 1;
    if (constr_cell_ind != last_ind) cell.constraints[constr_cell_ind] = cell.constraints[last_ind];

    cell.constraints.pop_back();
}

// This function is used to identify the edges that have been marked as to
// remove (with UINT64_MAX) from a certain faces[face].edge vector
inline bool remove_edge_from_face(uint64_t edge_face_ind)
{
    return edge_face_ind == UINT64_MAX;
}

//  Input: index of the splitted BSPface (w.r.t. vector faces) that is going
//         to be turned into the downSubface: face_ind,
//         index of the upSubface (w.r.t. vector faces): newFace_ind.
// Output: nothing.
void BSPcomplex::edgesPartition(uint64_t face_ind, uint64_t newFace_ind)
{
    // Partition of the edges between upSubface and downSubface.
    BSPface& face = faces[face_ind];

    // Search the first edge having first vertex ==0 and second vertex <0
    // (in face order). This will be the first in upFace.
    uint64_t e;
    for (e = 0; e < face.edges.size(); e++) {
        const uint64_t edge_ind = face.edges[e];
        const uint64_t next_edge_ind = face.edges[((e + 1) == face.edges.size()) ? (0) : (e + 1)];
        const BSPedge& edge = edges[edge_ind];
        const BSPedge& nedge = edges[next_edge_ind];
        const uint32_t comm_vert =
            (edge.vertices[0] == nedge.vertices[0] || edge.vertices[0] == nedge.vertices[1]) ? 0
                                                                                             : 1;
        if (vrts_orBin[edge.vertices[comm_vert]] < 0 && vrts_orBin[edge.vertices[!comm_vert]] == 0)
            break;
    }
    if (e == face.edges.size()) ip_error("pippo1\n");

    std::rotate(face.edges.begin(), face.edges.begin() + e, face.edges.end());

    // Search the last edge belonging to upFace
    for (e = 1; e < face.edges.size(); e++) {
        const uint64_t edge_ind = face.edges[e];
        const BSPedge& edge = edges[edge_ind];
        if (vrts_orBin[edge.vertices[0]] == 0 || vrts_orBin[edge.vertices[1]] == 0) break;
    }
    if (e == face.edges.size()) ip_error("pippo2\n");
    e++;
    // Move tail edges to new face
    faces[newFace_ind].edges.assign(face.edges.begin() + e, face.edges.end());
    face.edges.resize(e);
    for (uint64_t ei : faces[newFace_ind].edges) edges[ei].conn_face_0 = newFace_ind;

    // std::move(face.edges.begin() + e, face.edges.end(), faces[newFace_ind].edges.begin());

    // for(uint64_t e=0; e<face.edges.size(); e++){
    //   const uint64_t edge_ind = face.edges[e];
    //   BSPedge& edge = edges[edge_ind];

    //  if(vrts_orBin[ edge.vertices[0] ]>0 || vrts_orBin[ edge.vertices[1] ] >0 )
    //      move_edge(e, face_ind, newFace_ind);
    //}

    //// Remove all edges of face.edges marked as UINT64_MAX
    // face.edges.erase(
    //   std::remove_if(face.edges.begin(), face.edges.end(), remove_edge_from_face),
    //   face.edges.end() );
}

//  Input: index of the splitted BSPcell (w.r.t. vector cells) that is going
//         to be turned into the downSubcell: cell_ind,
//         index of the upSubcell (w.r.t. vector cells): newCell_ind,
//         vector of the indices of the vertices (w.r.t. vector vertices)
//         of the BSPcell to which the BSPface belong to: cell_vrts.
// Output: nothing.
void BSPcomplex::facesPartition(
    uint64_t cell_ind,
    uint64_t newCell_ind,
    const vector<uint32_t>& cell_vrts)
{
    // Faces whose indices (w.r.t. vector faces) are listed in
    // cells[cell_ind].faces have to be partitioned between
    // upSubcell (i.e. cells[newCell]) and downSubcell (i.e. cells[cell]).
    BSPcell& cell = cells[cell_ind];
    uint64_t num_faces = cell.faces.size();
    uint64_t face_ind;
    for (uint64_t f = 0; f < num_faces; f++) {
        face_ind = cell.faces[f];
        BSPface& face = faces[face_ind];

        vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
        list_faceVertices(face, face_vrts);

        // Face vertices disposition w.r.t. constraint-plane.
        uint32_t vrtsOVER, vrtsUNDER, vrtsON;
        count_vrt_orBin(face_vrts, &vrtsOVER, &vrtsUNDER, &vrtsON);

#ifdef DEBUG_BSP_DEEP
        print_BSPface_edges(edges, face_ind, face.edges);
        print_vrt_orBin(vrts_orBin, face_vrts);
#endif

        // It is impossible that a face has:
        // - all vertices ON the constraint-plane,
        // - two (or more) vertices on opposite sides w.r.t. the constraint-plane.

#ifdef DEBUG_BSP
        if (vrtsOVER == 0 && vrtsUNDER == 0)
            printf(
                "\n[BSP.cpp]BSPcomplex::facesPartition: ERROR face have "
                "all vertices on the constraint plane.\n");
        if (vrtsOVER > 0 && vrtsUNDER > 0)
            printf(
                "\n[BSP.cpp]BSPcomplex::facesPartition: ERROR face intersects "
                "the constraint plane.\n");
#endif

        // IF one of the face vertices is OVER the constraint-plane (vrtsOVER>0),
        // the face is assigned to the upSubcell (i.e. cells[newCell]).
        if (vrtsOVER > 0) {
            move_face(f, cell_ind, newCell_ind);
            f--;
            num_faces--;
        }
        // OTHERWISE
        // faces[face_ind] goes to down-subcell (i.e. remain to cells[cell_ind])
        // indeed, all its vertices have non-positive cell_vrts_orient.

#ifdef DEBUG_BSP_DEEP
        if (vrtsOVER > 0)
            printf("\n\tface #%llu goes to up-subcell (cell #%llu).\n", face_ind, newCell_ind);
        else
            printf("\n\tface #%llu goes to down-subcell (cell #%llu).\n", face_ind, cell_ind);
#endif
    }
}

//  Input: index of the constraint that has divided the original BSPcell
//         cells[down_cell] in to the current sub-cells cells[up_cell] and
//         the sub-cells cells[down_cell]: ref_constr,
//         index of the down sub-cell (w.r.t. vector cells): down_cell_ind,
//         index of the up sub-cell (w.r.t. vector cells): up_cell_ind,
//         vector of the indices of the vertices (w.r.t. vector vertices)
//         of the down sub-cell (before splitting): cell_vrts,
// Output: nothing.
void BSPcomplex::constraintsPartition(
    uint32_t ref_constr,
    uint64_t down_cell_ind,
    uint64_t up_cell_ind,
    const vector<uint32_t>& cell_vrts)
{
    // At this point down sub-cell has all the constraints,
    // while up sub-cell has none.
    BSPcell& down_cell = cells[down_cell_ind];
    BSPcell& up_cell = cells[up_cell_ind];

#ifdef DEBUG_BSP_DEEP
    printf("\tCurrent constraint is: ");
    print_constraint(constraints_vrts, ref_constr);
    vector<uint32_t> vrts_to_print;
    for (uint32_t v = 0; v < cell_vrts.size(); v++)
        if (vrts_orBin[cell_vrts[v]] >= 0) vrts_to_print.push_back(cell_vrts[v]);
    print_BSPcell_vrts(vrts_to_print, up_cell_ind);
    vrts_to_print.clear();
    for (uint32_t v = 0; v < cell_vrts.size(); v++)
        if (vrts_orBin[cell_vrts[v]] <= 0) vrts_to_print.push_back(cell_vrts[v]);
    print_BSPcell_vrts(vrts_to_print, down_cell_ind);
    printf("\tConstraints intersecting (original)cell #%llu are:\n", down_cell_ind);
    for (uint32_t c = 0; c < down_cell.constraints.size(); c++) {
        printf("\t");
        print_constraint(constraints_vrts, down_cell.constraints[c]);
    }
    printf("\n");
#endif

    // 3 vertices of the ref_constraint seen as triangle (k0,k1,k2).
    uint32_t ref_constr_ID = 3 * ref_constr;
    uint32_t k0 = constraints_vrts[ref_constr_ID];
    uint32_t k1 = constraints_vrts[ref_constr_ID + 1];
    uint32_t k2 = constraints_vrts[ref_constr_ID + 2];

    uint64_t num_constr = down_cell.constraints.size();
    uint32_t constr, constr_ID;
    vector<uint32_t> constr_vrts(3, UINT32_MAX);
    for (uint32_t c = 0; c < num_constr; c++) {
        constr = down_cell.constraints[c];
        constr_ID = 3 * constr;
        constr_vrts[0] = constraints_vrts[constr_ID];
        constr_vrts[1] = constraints_vrts[constr_ID + 1];
        constr_vrts[2] = constraints_vrts[constr_ID + 2];

        // commFace_vrts disposition w.r.t. constr vertices.
        vrts_orient_wrtPlane(constr_vrts, k0, k1, k2, 2);
        uint32_t vrtsOVER, vrtsUNDER, vrtsON;
        count_vrt_orBin(constr_vrts, &vrtsOVER, &vrtsUNDER, &vrtsON);

#ifdef DEBUG_BSP_DEEP
        printf(
            "\n\t orient3d of constraint %u vertices w.r.t. plane for ",
            down_cell.constraints[c]);
        print_constraint(constraints_vrts, ref_constr);
        printf("\n");
        print_vrt_orBin(vrts_orBin, constr_vrts);
        printf("\n");
#endif

        // If constr and commFace_vrts define the same plane, remove constr since
        // the cut will not produce further cell-split.
        if (vrtsOVER == 0 && vrtsUNDER == 0) {
            remove_constraint(c, down_cell_ind);
            c--;
            num_constr--;

#ifdef DEBUG_BSP_DEEP
            printf(
                "\tConsraints %u and %u define the same plane, constraint %u have "
                "been removed.\n",
                ref_constr,
                constr,
                constr);
            if (vrtsON == 0)
                printf(
                    "\n[BSP.cpp]BSPcomplex::constraintsPartition: ERROR vertices "
                    "of consraints %u have an undefined position wrt constraint %u.\n",
                    constr,
                    ref_constr);
#endif

            continue; // jump to next constraint.
        }

        const bool up = (vrtsOVER > 0);
        const bool down = (vrtsUNDER > 0);

#ifdef DEBUG_BSP
        if (!up && !down)
            printf(
                "[BSP.cpp]BSPcomplex::constraintsPartition: WARNING constraint "
                "#%u will be removed from down sub-cell #%llu, there are no "
                "intersection with both sub-cells.\n",
                constr,
                down_cell_ind);
#endif

        if (up) up_cell.constraints.push_back(constr);
        if (!down) {
            remove_constraint(c, down_cell_ind);
            c--;
            num_constr--;
        }
        // OTHERWISE
        // the constraint indexed as constr goes to down-subcell (i.e. remain
        // to cells[down_cell]).

#ifdef DEBUG_BSP_DEEP
        if (up && down)
            printf(
                "\tconstraint #%u goes to both subcell (cell "
                "#%llu and #%llu).\n",
                constr,
                up_cell_ind,
                down_cell_ind);

        if (up && !down)
            printf(
                "\tconstraint #%u goes to up-subcell (cell "
                "#%llu).\n",
                constr,
                up_cell_ind);
        if (!up && down)
            printf(
                "\tconstraint #%u goes to down-subcell (cell "
                "#%llu).\n",
                constr,
                down_cell_ind);
#endif
    }
}

//
//
void BSPcomplex::add_edgeToOrdFaceEdges(BSPface& face, uint64_t newEdge_ind)
{
    BSPedge& newEdge = edges[newEdge_ind];
    uint64_t edge_ind, num_faceEdges = face.edges.size();
    uint32_t n0, n1, e0, e1;
    n0 = newEdge.vertices[0];
    n1 = newEdge.vertices[1];

    for (uint64_t e = 0; e < num_faceEdges; e++) {
        edge_ind = face.edges[e];
        e0 = edges[edge_ind].vertices[0];
        e1 = edges[edge_ind].vertices[1];

        if (CONSECUTIVE_EDGES(n0, n1, e0, e1)) {
            // Special case: edge is the first vector element
            if (e == 0) {
                BSPedge& cons = edges[face.edges.back()];
                if (CONSECUTIVE_EDGES(n0, n1, cons.vertices[0], cons.vertices[1]))
                    face.edges.push_back(newEdge_ind);
                else { // cons is faces[face].edges[1]
                    face.edges.push_back(edge_ind);
                    face.edges[0] = newEdge_ind;
                }

                return;
            }

            // Special case: edge is the last vector element
            if (e == num_faceEdges - 1) {
                BSPedge& cons = edges[face.edges[0]];
                if (CONSECUTIVE_EDGES(n0, n1, cons.vertices[0], cons.vertices[1]))
                    face.edges.push_back(newEdge_ind);
                else { // cons is faces[face].edges.size() -2
                    face.edges.push_back(edge_ind);
                    face.edges[num_faceEdges - 1] = newEdge_ind;
                }

                return;
            }

            // General case
            BSPedge& cons = edges[face.edges[e + 1]];
            if (CONSECUTIVE_EDGES(n0, n1, cons.vertices[0], cons.vertices[1]))
                face.edges.insert(face.edges.begin() + e + 1, newEdge_ind);
            else // cons = edges[ faces[face].edges[e-1] ]
                face.edges.insert(face.edges.begin() + e, newEdge_ind);

            return;
        }
    }
}

//  Input: index of the constraint splitting the BSPface: constr,
//         index of the splitted BSPface (w.r.t. vector faces) that is going
//         to be turned into the downSubface: face_ind,
//         index of the upSubface (w.r.t. vector cells): newFace_ind,
//         vector of the indices of the 2 vertices (w.r.t. vector vertices)
//         of the original BSPface faces[face] that lie
//         on the constraint-plane: endpts.
// Output: nothing.
void BSPcomplex::add_commonEdge(
    uint32_t constr,
    uint64_t face_ind,
    uint64_t newFace_ind,
    const uint32_t* endpts)
{
    BSPface& face = faces[face_ind];
    BSPface& newFace = faces[newFace_ind];
    // The new faces "face" and "newFace", originated by splitting the original
    // "face" with the constraint-plane (constr), have a common edge.
    // The endpoint of the common edge are those that have vrt_orBin=0, i.e.
    // endpts[0], endpts[1].

    // New edge originated by the face-constraint intersection:
    // it is defined by the intersection of two planes (face and constraint)
    uint32_t constr_ID = 3 * constr;
    edges.push_back(BSPedge(
        endpts[0],
        endpts[1],
        face.meshVertices[0],
        face.meshVertices[1],
        face.meshVertices[2],
        constraints_vrts[constr_ID],
        constraints_vrts[constr_ID + 1],
        constraints_vrts[constr_ID + 2]));
    uint64_t commEdge_ind = edges.size() - 1;
    BSPedge& commEdge = edges[commEdge_ind];
    // commEdge.conn_faces.push_back(newFace_ind);
    // commEdge.conn_faces.push_back(face_ind);
    commEdge.conn_face_0 = face_ind;

    // Add an element to global vector edge_visit.
    edge_visit.push_back(0);

    // Add the common edge to face and newFace.
    // add_edgeToOrdFaceEdges(face, commEdge_ind);
    // add_edgeToOrdFaceEdges(newFace, commEdge_ind);
    face.edges.push_back(commEdge_ind);
    newFace.edges.push_back(commEdge_ind);
}

//  Input: an empty BSPface, created by intersecating a BSPcell with a
//         constraint: face,
//         the indices of all BSPedges (w.r.t. vector edges) that bound the
//         BSPface: edges_ind.
// Output: nothing.
void BSPcomplex::add_edges_toCommFaceEdges(BSPface& face, const vector<uint64_t>& edges_ind)
{
    // Find face vertices: since the face boundary is closed,
    //                     there are as many vertices as are the edges.
    vector<uint32_t> face_vrts(edges_ind.size(), UINT32_MAX);
    uint64_t edge_ind;
    uint32_t e0, e1, fv_ind = 0;
    for (uint64_t e = 0; e < edges_ind.size(); e++) {
        BSPedge& edge = edges[edges_ind[e]];
        e0 = edge.vertices[0];
        if (vrts_visit[e0] == 0) {
            face_vrts[fv_ind++] = e0;
            vrts_visit[e0] = 1;
        }
        e1 = edge.vertices[1];
        if (vrts_visit[e1] == 0) {
            face_vrts[fv_ind++] = e1;
            vrts_visit[e1] = 1;
        }
    }

    // Relate each face vertex with its two incident edges.
    vector<uint64_t> rel_VE(2 * face_vrts.size(), UINT64_MAX);

    // Set vrts_visit of face_vrts indices in order to memory positions.
    for (uint32_t u = 0; u < face_vrts.size(); u++) vrts_visit[face_vrts[u]] = UINT32_MAX;

    // Fill rel_VE
    uint32_t pos = 0;
    for (uint64_t e = 0; e < edges_ind.size(); e++) {
        edge_ind = edges_ind[e];
        BSPedge& edge = edges[edge_ind];
        e0 = edge.vertices[0];
        e1 = edge.vertices[1];
        if (vrts_visit[e0] == UINT32_MAX) {
            rel_VE[2 * pos] = edge_ind;
            vrts_visit[e0] = pos++;
        } else
            rel_VE[2 * vrts_visit[e0] + 1] = edge_ind;

        if (vrts_visit[e1] == UINT32_MAX) {
            rel_VE[2 * pos] = edge_ind;
            vrts_visit[e1] = pos++;
        } else
            rel_VE[2 * vrts_visit[e1] + 1] = edge_ind;
    }

    // Fill faces[face_ind].edges
    uint32_t num_ins_vrts = 0, next_vrt = face_vrts[0];
    edge_ind = rel_VE[2 * vrts_visit[next_vrt]];

    face.edges.resize(edges_ind.size());
    while (num_ins_vrts < face_vrts.size()) {
        if (edge_visit[edge_ind] == 0) {
            edge_visit[edge_ind] = 1;
            face.edges[num_ins_vrts++] = edge_ind;

            e0 = edges[edge_ind].vertices[0];
            e1 = edges[edge_ind].vertices[1];
            if (next_vrt == e0)
                next_vrt = e1;
            else
                next_vrt = e0;
            edge_ind = rel_VE[2 * vrts_visit[next_vrt]];
        } else {
            if (edge_ind == rel_VE[2 * vrts_visit[next_vrt]])
                edge_ind = rel_VE[2 * vrts_visit[next_vrt] + 1];
            else
                edge_ind = rel_VE[2 * vrts_visit[next_vrt]];
        }
    }

    // Reset vrts_visit and edge_visit
    for (uint32_t u = 0; u < face_vrts.size(); u++) vrts_visit[face_vrts[u]] = 0;
    for (uint64_t u = 0; u < edges_ind.size(); u++) edge_visit[edges_ind[u]] = 0;
}

//  Input: index of the constraint splitting the BSPcell: constr,
//         index of the splitted BSPcell (w.r.t. vector cells) that is going
//         to be turned into the downSubcell: cell_ind,
//         index of the upSubcell (w.r.t. vector cells): newCell_ind,
//         vector of the indices of the vertices (w.r.t. vector vertices)
//         of the BSPcell to which the BSPface belong to: cell_vrts.
// Output: nothing.
void BSPcomplex::add_commonFace(
    uint32_t constr,
    uint64_t cell_ind,
    uint64_t newCell_ind,
    const vector<uint32_t>& cell_vrts,
    const vector<uint64_t>& cell_edges)
{
    // Common face between up-subcell and down-subcell: the edge of that face
    // are those of cells[cell_ind] that have vrts_orBin = 0.
    uint32_t constr_ID = 3 * constr;
    COLOUR_T colour = GREY;
    if (is_virtual(constr)) colour = WHITE;
    faces.push_back(BSPface(
        constraints_vrts[constr_ID],
        constraints_vrts[constr_ID + 1],
        constraints_vrts[constr_ID + 2],
        cell_ind,
        newCell_ind,
        colour));
    uint64_t face_ind = faces.size() - 1;

    uint64_t edge_ind, num_commFace_edges = 0;
    // Count edges whose endpoints are both on the constraint-plane.
    for (uint64_t e = 0; e < cell_edges.size(); e++) {
        BSPedge& edge = edges[cell_edges[e]];
        if (vrts_orBin[edge.vertices[0]] == 0 && vrts_orBin[edge.vertices[1]] == 0)
            num_commFace_edges++;
    }
    // Fill a vector with those edges.
    vector<uint64_t> commFace_edges(num_commFace_edges, UINT64_MAX);
    num_commFace_edges = 0;
    for (uint64_t e = 0; e < cell_edges.size(); e++) {
        edge_ind = cell_edges[e];
        BSPedge& edge = edges[edge_ind];
        if (vrts_orBin[edge.vertices[0]] == 0 && vrts_orBin[edge.vertices[1]] == 0)
            commFace_edges[num_commFace_edges++] = edge_ind;
    }

    add_edges_toCommFaceEdges(faces[face_ind], commFace_edges);

    // for (uint64_t e = 0; e < commFace_edges.size(); e++)
    //     edges[commFace_edges[e]].conn_faces.push_back(face_ind);
    for (uint64_t e = 0; e < commFace_edges.size(); e++)
        edges[commFace_edges[e]].conn_face_0 = face_ind;

    cells[cell_ind].faces.push_back(face_ind);
    cells[newCell_ind].faces.push_back(face_ind);

    fixCommonFaceOrientation(face_ind);
}

inline bool edges_ShareCommonPlanes(const BSPedge& a, const BSPedge& b)
{
    const uint32_t* mva = a.meshVertices;
    const uint32_t* mvb = b.meshVertices;
    return (
        mva[0] == mvb[0] && mva[1] == mvb[1] && mva[2] == mvb[2] && mva[3] == mvb[3] &&
        mva[4] == mvb[4] && mva[5] == mvb[5]);
}

void BSPcomplex::fixCommonFaceOrientation(uint64_t cf_id)
{
    BSPface& f = faces[cf_id];
    const uint32_t* mv = f.meshVertices;
    double mvc[9]; // Coords of the original input triangle
    vertices[mv[0]]->getApproxXYZCoordinates(mvc[0], mvc[1], mvc[2]);
    vertices[mv[1]]->getApproxXYZCoordinates(mvc[3], mvc[4], mvc[5]);
    vertices[mv[2]]->getApproxXYZCoordinates(mvc[6], mvc[7], mvc[8]);
    int xyz = genericPoint::maxComponentInTriangleNormal(
        mvc[0],
        mvc[1],
        mvc[2],
        mvc[3],
        mvc[4],
        mvc[5],
        mvc[6],
        mvc[7],
        mvc[8]);

    int ori0;
    if (xyz == 2)
        ori0 = orient2d(mvc[0], mvc[1], mvc[3], mvc[4], mvc[6], mvc[7]);
    else if (xyz == 0)
        ori0 = orient2d(mvc[1], mvc[2], mvc[4], mvc[5], mvc[7], mvc[8]);
    else
        ori0 = orient2d(mvc[2], mvc[0], mvc[5], mvc[3], mvc[8], mvc[6]);

    const BSPedge& edge0 = edges[f.edges.back()];
    const BSPedge& edge1 = edges[f.edges[0]];
    uint32_t vid = consecEdges_common_endpt(
        edge0.vertices[0],
        edge0.vertices[1],
        edge1.vertices[0],
        edge1.vertices[1]);

    genericPoint* v0 = vertices[vid];
    if (vid == edge1.vertices[0])
        vid = edge1.vertices[1];
    else
        vid = edge1.vertices[0];
    genericPoint* v1 = vertices[vid];
    for (uint64_t e = 1; e < f.edges.size(); e++) {
        const BSPedge& edge = edges[f.edges[e]];
        if (vid == edge.vertices[0])
            vid = edge.vertices[1];
        else
            vid = edge.vertices[0];
        if (edges_ShareCommonPlanes(edge1, edge)) continue;
        genericPoint* v2 = vertices[vid];
        int ori = ori0 * genericPoint::orient2D(*v0, *v1, *v2, xyz);
        if (ori < 0) return;
        if (ori > 0) {
            if (f.conn_cells[1] == UINT64_MAX) {
                printf("Mmh... this should not happen\n");
                return;
            }
            std::swap(f.conn_cells[0], f.conn_cells[1]);
            return;
        }
    }
    ip_error("Degenerate face\n");
}

//  Input: a BSPedge: edge,
//         index of a constraint (w.r.t. vector constraints_vrt): constr.
// Output: returns the index of the new vertex (w.r.t. vector vertices).
// Note: a LPI (Line Plane Intersection) vertex can be generated only as
//       intesection beween a constraint-plane and an edge that is part of the
//       Delaunay triangulation (or a pice of Delaunay edge).
uint32_t BSPcomplex::add_LPIvrt(const BSPedge& edge, uint32_t constr)
{
    uint32_t constr_ID = 3 * constr;
    genericPoint* c0 = vertices[constraints_vrts[constr_ID]];
    genericPoint* c1 = vertices[constraints_vrts[constr_ID + 1]];
    genericPoint* c2 = vertices[constraints_vrts[constr_ID + 2]];
    genericPoint* e0 = vertices[edge.meshVertices[0]];
    genericPoint* e1 = vertices[edge.meshVertices[1]];

#ifdef DEBUG_BSP
    if (!(e0->isExplicit3D()) || !(e1->isExplicit3D()) || !(c0->isExplicit3D()) ||
        !(c1->isExplicit3D()) || !(c2->isExplicit3D()))
        printf(
            "\n[BSP.cpp]BSPcomplex::add_LPIvrt: ERROR explicitPoint3D are "
            "expected (LPI).\n");
#endif

    vertices.push_back(new implicitPoint3D_LPI(
        e0->toExplicit3D(),
        e1->toExplicit3D(),
        c0->toExplicit3D(),
        c1->toExplicit3D(),
        c2->toExplicit3D()));

    // Add new element to global vectors vrts_orBin and vrts_visit.
    vrts_orBin.push_back(2);
    vrts_visit.push_back(0);

    return (uint32_t)(vertices.size() - 1);
}

// Return TRUE if the intersection of the sets {p,q,r} and {t,u,v} has at least two elements.
// Fill 'e' with the intersecting elements.
inline bool twoEqualVertices(
    uint32_t p,
    uint32_t q,
    uint32_t r,
    uint32_t s,
    uint32_t t,
    uint32_t u,
    uint32_t* e)
{
    int i = 0;
    if (p == s || p == t || p == u) e[i++] = p;
    if (q == s || q == t || q == u) e[i++] = q;
    if (r == s || r == t || r == u) e[i++] = r;
    return (i >= 2);
}

//  Input: a BSPedge: edge,
//         index of a constraint (w.r.t. vector constraints_vrt): constr.
// Output: returns the index of the new vertex (w.r.t. vector vertices).
// Note: a TPI (Triple Plane Intersection) vertex can be generated only as
//       intesection beween a constraint-plane and an edge that is defined as
//       the intersection between 3 constraint-planes (i.e. not included in
//       the original Delaunay triangulation).
uint32_t BSPcomplex::add_TPIvrt(const BSPedge& edge, uint32_t constr)
{
    const uint32_t constr_ID = 3 * constr;

    const uint32_t ic0 = constraints_vrts[constr_ID];
    const uint32_t ic1 = constraints_vrts[constr_ID + 1];
    const uint32_t ic2 = constraints_vrts[constr_ID + 2];
    const uint32_t ie0 = edge.meshVertices[0];
    const uint32_t ie1 = edge.meshVertices[1];
    const uint32_t ie2 = edge.meshVertices[2];
    const uint32_t ie3 = edge.meshVertices[3];
    const uint32_t ie4 = edge.meshVertices[4];
    const uint32_t ie5 = edge.meshVertices[5];

    // If two of the three triangles share two vertices -> create an LPI
    uint32_t comm[3];
    if (twoEqualVertices(ic0, ic1, ic2, ie0, ie1, ie2, comm))
        vertices.push_back(new implicitPoint3D_LPI(
            vertices[comm[0]]->toExplicit3D(),
            vertices[comm[1]]->toExplicit3D(),
            vertices[ie3]->toExplicit3D(),
            vertices[ie4]->toExplicit3D(),
            vertices[ie5]->toExplicit3D()));
    else if (twoEqualVertices(ic0, ic1, ic2, ie3, ie4, ie5, comm))
        vertices.push_back(new implicitPoint3D_LPI(
            vertices[comm[0]]->toExplicit3D(),
            vertices[comm[1]]->toExplicit3D(),
            vertices[ie0]->toExplicit3D(),
            vertices[ie1]->toExplicit3D(),
            vertices[ie2]->toExplicit3D()));
    else if (twoEqualVertices(ie3, ie4, ie5, ie0, ie1, ie2, comm))
        vertices.push_back(new implicitPoint3D_LPI(
            vertices[comm[0]]->toExplicit3D(),
            vertices[comm[1]]->toExplicit3D(),
            vertices[ic0]->toExplicit3D(),
            vertices[ic1]->toExplicit3D(),
            vertices[ic2]->toExplicit3D()));
    else {
        vertices.push_back(new implicitPoint3D_TPI(
            vertices[ie0]->toExplicit3D(),
            vertices[ie1]->toExplicit3D(),
            vertices[ie2]->toExplicit3D(),
            vertices[ie3]->toExplicit3D(),
            vertices[ie4]->toExplicit3D(),
            vertices[ie5]->toExplicit3D(),
            vertices[ic0]->toExplicit3D(),
            vertices[ic1]->toExplicit3D(),
            vertices[ic2]->toExplicit3D()));
    }

    // Add new element to global vectors vrts_orBin and vrts_visit.
    vrts_orBin.push_back(2);
    vrts_visit.push_back(0);

    return (uint32_t)(vertices.size() - 1);
}


// 2 funzioni
// 1 -

// Given a cell 'c', one of its faces 'f0', and one of f0 edges 'e0'
// return the face in 'c' that shares 'e0' with 'f0'
uint64_t BSPcomplex::getOppositeEdgeFace(const uint64_t e0, const uint64_t f0, const uint64_t c)
{
    const std::vector<uint64_t>& cfaces = cells[c].faces;
    for (uint64_t fid : cfaces)
        if (fid != f0) {
            const std::vector<uint64_t>& fedges = faces[fid].edges;
            for (uint64_t e : fedges)
                if (e == e0) return fid;
        }
    return UINT64_MAX;
}

static inline uint64_t oppositeCellId(const uint64_t c_id, const BSPface& f)
{
    if (f.conn_cells[0] == c_id)
        return f.conn_cells[1];
    else
        return f.conn_cells[0];
}

void BSPcomplex::makeEFrelation(const uint64_t e_id, std::vector<uint64_t>& ef)
{
    const BSPedge& e = edges[e_id];
    uint64_t f = e.conn_face_0;
    uint64_t c = faces[f].conn_cells[0];
    ef.push_back(f);

    for (;;) {
        f = getOppositeEdgeFace(e_id, f, c);
        if (f == e.conn_face_0)
            return;
        else
            ef.push_back(f);
        c = oppositeCellId(c, faces[f]);
        if (c == UINT64_MAX) break;
    }

    f = e.conn_face_0;
    if ((c = faces[f].conn_cells[1]) == UINT64_MAX) return;

    for (;;) {
        f = getOppositeEdgeFace(e_id, f, c);
        ef.push_back(f);
        c = oppositeCellId(c, faces[f]);
        if (c == UINT64_MAX) return;
    }
}
//  Input: a BSPedge: edge,
//         index of a constraint that intersects the edge interior: constr.
// Output: the index of the sub-edge (w.r.t. vector edges) originated by the
//         intersection between the edge and the constraint, whose endpoints
//         have non-negative orient3D w.r.t. the constraint plane.
void BSPcomplex::splitEdge(uint64_t edge_id, uint32_t constr)
{
    BSPedge& edge = edges[edge_id];

    std::vector<uint64_t> ef;
    makeEFrelation(edge_id, ef);

    // Add the point in which the edge intersects the constraint-plane
    // to the vector vertices.
    // edge is a Delauany mesh edge -> Line Plane Intersection.
    // edge is a 2 constraints intersection -> Triple Plane Intersection.
    uint32_t new_point;
    if (edge.meshVertices[2] == UINT32_MAX)
        new_point = add_LPIvrt(edge, constr);
    else
        new_point = add_TPIvrt(edge, constr);

    // Split the edge <e0,e1> -> <e0,new_point> + <new_point,e1>
    // edge <- <new_point,e1>
    // new_edge = edges[edges.size()-1] <- <e0,new_point>.
    edges.push_back(edge.split(new_point));
    uint64_t new_edge_ind = edges.size() - 1;
    // Note. new edge is created with the same conn_faces of the old edge.

    // Add new element to edge_visit
    edge_visit.push_back(0);

    // Add new edge to all conn_faces of old edge.
    // for(uint64_t f=0; f< edges[edge_id].conn_faces.size(); f++)
    //  add_edgeToOrdFaceEdges(faces[edges[edge_id].conn_faces[f] ], new_edge_ind);

    for (uint64_t f : ef) add_edgeToOrdFaceEdges(faces[f], new_edge_ind);
}

//  Input: index of a BSPface: face_ind,
//         index of a constraint that intersects the face interior: constr,
//         index of the BSPcell to which the BSPface belongs to: cell_ind,
//         a vector with the indices of the face vertices: face_vrts.
// Output: nothing.
void BSPcomplex::splitFace(
    uint64_t face_ind,
    uint32_t constr,
    uint64_t cell_ind,
    const vector<uint32_t>& face_vrts)
{
#ifdef DEBUG_BSP_DEEP
    printf("\n\tDividing face #%llu with constraint %u.\n", face_ind, constr);
#endif

    BSPface& face = faces[face_ind];

    // The face faces[face] is divided in two subfaces:
    // the up-subface and the down-subface.
    // The edges of faces[face] have to be partitioned between them.
    // up-subface inherits edges whose vertices have non-negative vrt_orBin.
    // down-subface inherits edges whose vertices have non-positive vrt_orBin.

    // Face faces[face], after splitting, becomes the down-subface.
    // New BSPface is the up-subface of the original faces[face].
    faces.push_back(BSPface(
        face.meshVertices[0],
        face.meshVertices[1],
        face.meshVertices[2],
        face.conn_cells[0],
        face.conn_cells[1],
        face.colour,
        face.coplanar_constraints));
    uint64_t newFace_ind = faces.size() - 1;
    // Note. the edge of the new face (i.e. up-subface) will be assigned later.

#ifdef DEBUG_BSP_DEEP
    printf(
        "\tSplit face %llu -> up-subface (face #%llu) + down-subface "
        "(face #%llu)\n",
        face_ind,
        newFace_ind,
        face_ind);
#endif

    // Add the new face to its adjacent cell (the same of faces[face]).
    // If it is a convex-hull face (i.e. conn_cells[1] = UINT64_MAX),
    // there is no adjacent cell.

    cells[faces[face_ind].conn_cells[0]].faces.push_back(newFace_ind);
    if (!IS_GHOST_CELL(faces[face_ind].conn_cells[1]))
        cells[faces[face_ind].conn_cells[1]].faces.push_back(newFace_ind);

    // Find face vertices that have vrts_orBin=0.
    uint32_t zero_vrts[2];
    uint32_t pos = 0;
    for (uint32_t v = 0; v < face_vrts.size(); v++)
        if (vrts_orBin[face_vrts[v]] == 0) zero_vrts[pos++] = face_vrts[v];

#ifdef DEBUG_BSP
    if (pos != 2)
        printf(
            "\n[BSP.cpp]BSPcomplex::splitFace: ERROR there must be exactly "
            "2 (not %u) face-vertices on the constraint plane.\n",
            pos);
#endif

    // Partition of the edges between up-subface and down-subface.
    edgesPartition(face_ind, newFace_ind);

    // The new faces faces[face_ind] and faces[newFace_ind] have a common edge.
    add_commonEdge(constr, face_ind, newFace_ind, zero_vrts);


    // if (!faceHasCorrectOrientation(this, face_ind)) ip_error("Wrong f1\n");
    // if (!faceHasCorrectOrientation(this, newFace_ind)) ip_error("Wrong f2\n");
}

//
//
void BSPcomplex::find_coplanar_constraints(
    uint64_t cell_ind,
    uint32_t constr,
    vector<uint32_t>& coplanar_c)
{
    BSPcell& cell = cells[cell_ind];
    uint32_t constr_ID = 3 * constr;
    uint32_t c0 = constraints_vrts[constr_ID];
    uint32_t c1 = constraints_vrts[constr_ID + 1];
    uint32_t c2 = constraints_vrts[constr_ID + 2];

    // Count coplanar constraints.
    uint32_t num_coplanar = 0;
    vector<uint32_t> k_vrts(3, UINT32_MAX);
    for (uint32_t k = 0; k < cell.constraints.size(); k++) {
        if (is_virtual(cell.constraints[k])) continue;
        uint32_t kID = 3 * cell.constraints[k];
        k_vrts[0] = constraints_vrts[kID];
        k_vrts[1] = constraints_vrts[kID + 1];
        k_vrts[2] = constraints_vrts[kID + 2];
        vrts_orient_wrtPlane(k_vrts, c0, c1, c2, 0);
        if (vrts_orBin[k_vrts[0]] == 0 && vrts_orBin[k_vrts[1]] == 0 &&
            vrts_orBin[k_vrts[2]] == 0) {
#ifdef DEBUG_BSP_DEEP
            printf(
                "cell #%llu: constraints %u and %u are aligned.\n",
                cell_ind,
                constr_ID / 3,
                kID / 3);
#endif

            num_coplanar++;
        }
    }
    if (!is_virtual(constr)) num_coplanar++;

    // Save coplanar constraints in a new vector and remove them from cell's list.
    if (num_coplanar > 0) {
        coplanar_c.resize(num_coplanar, UINT32_MAX);
        uint32_t pos = 0;
        for (uint32_t k = 0; k < cell.constraints.size(); k++) {
            if (is_virtual(cell.constraints[k])) continue;
            uint32_t kID = 3 * cell.constraints[k];
            k_vrts[0] = constraints_vrts[kID];
            k_vrts[1] = constraints_vrts[kID + 1];
            k_vrts[2] = constraints_vrts[kID + 2];

            if (vrts_orBin[k_vrts[0]] == 0 && vrts_orBin[k_vrts[1]] == 0 &&
                vrts_orBin[k_vrts[2]] == 0) {
                coplanar_c[pos++] = cell.constraints[k];
                cell.constraints[k] = cell.constraints.back();
                cell.constraints.pop_back();
                k--;
            }
        }
        if (!is_virtual(constr)) coplanar_c[pos++] = constr;
    }
}

//  Input: index of a BSPcell: cell_ind.
// Output: nothing.
// Note. it is assumed that the BScell is (only) convex,
//       strictly convexity is not guaranteed.
void BSPcomplex::splitCell(uint64_t cell_ind)
{
    BSPcell& cell = cells[cell_ind];
    auto& cts = cell.constraints;

    // Distinguish between two mutually exclusive cases:
    // CASE. NO SPLIT: only cell boundary elements (face or edge) lie on the
    //                 constraint-plane, while the other elements belong to
    //                 the same half-space w.r.t. the constraint-plane.
    // CASE. SPLIT-INTERIOR: constraint-plane pass through the cell interior.
    // Note. it is not possible that only a vertex lies on the constraint-plane.

    // Create a local richer data structure to avoid multilple extracions of cell
    // edges and vertices.
    //
    // This is built once and then reused across the NO SPLIT constraints below. Such a
    // constraint leaves the cell untouched -- no edge, face or vertex of it changes -- so
    // the lists stay valid and rebuilding them would reproduce exactly what we already
    // have. On inputs where most constraints handed to a cell do not actually cut it, that
    // rebuild dominates the run time.
    uint64_t num_cellEdges = count_cellEdges(cell);
    uint64_t num_cellVrts = count_cellVertices(cell, &num_cellEdges);
    vector<uint64_t> cell_edges(num_cellEdges, UINT64_MAX);
    vector<uint32_t> cell_vrts(num_cellVrts, UINT32_MAX);
    fill_cell_locDS(cell, cell_edges, cell_vrts);

    uint32_t constr = UINT32_MAX, c0 = 0, c1 = 0, c2 = 0;
    uint32_t vrtsON, vrtsOVER, vrtsUNDER;
    vector<uint32_t> coplanar_constr;
    bool splits = false;

    // Take constraints from the back, exactly as before; the only difference is that a
    // constraint which turns out not to cut the cell is now skipped here instead of
    // returning to the caller, which would call us straight back to rebuild the same
    // lists. The order in which constraints are examined, and everything done to each one,
    // are unchanged.
    while (!cts.empty()) {
        // Extract the last contraint that intersect the cell and remove it from
        // the list. The cell will be splitted by that constraint.
        constr = cts.back();
        cts.pop_back();

        // Virtual constraints should be used only when non-virtual constraints are over
        if (is_virtual(constr)) {
            for (size_t i = 0; i < cts.size(); i++)
                if (!is_virtual(cts[i])) {
                    std::swap(cts[i], constr);
                    break;
                }
        }

        const uint32_t constr_ID = 3 * constr;
        c0 = constraints_vrts[constr_ID];
        c1 = constraints_vrts[constr_ID + 1];
        c2 = constraints_vrts[constr_ID + 2];

        // Search for coplanar constraints.
        coplanar_constr.clear();
        find_coplanar_constraints(cell_ind, constr, coplanar_constr);

        // Compute the orientation of cell vertices w.r.t. the constraint plane.
        vrts_orient_wrtPlane(cell_vrts, c0, c1, c2, 1);

        // Analysis of cell vertices disposition w.r.t. constraint-plane.
        count_vrt_orBin(cell_vrts, &vrtsOVER, &vrtsUNDER, &vrtsON);

        // CASE. NO SPLIT:
        // at least two cell_vrts_or are 0, the other (!=0) have the same signe.
        // The cell is unchanged, so keep the local data structure and try the next one.
        if (vrtsUNDER == 0 || vrtsOVER == 0) continue;

        splits = true;
        break;
    }
    if (!splits) return;

    // (else) CASE. SPLIT-INTERIOR:
    // at least two cell_vrts_or have opposite signe.

    // Split cell edges whose endpoints have opposite cell_vrts_or signe.

    for (uint64_t e = 0; e < cell_edges.size(); e++) {
        BSPedge& edge = edges[cell_edges[e]];
        if (constraint_innerIntersects_edge(edge, cell_vrts)) {
            splitEdge(cell_edges[e], constr); // Here the edge is splitted.
            // Add the new vertex to cell_vrts, compute its orient3D w.r.t. the
            // constraint-plane and add it to cell_vrts_or.
            uint32_t new_vrt = (uint32_t)(vertices.size() - 1);
            cell_vrts.push_back(new_vrt);
            cell_edges.push_back(edges.size() - 1);
            vrts_orBin[new_vrt] = 0;

#ifdef DEBUG_BSP_DEEP
            vector<uint32_t> vrt_to_print;
            vrt_to_print.push_back(new_vrt);
            printf("\n");
            print_vrt_orBin(vrts_orBin, vrt_to_print);
#endif
        }
    }
    // Now all the BSPcell edges are over or under the constraint.

    // Split cell-faces that have at lesat two vertices with opposite
    // cell_vrts_or signe.

    uint64_t num_faces = cell.faces.size();
    for (uint64_t f = 0; f < num_faces; f++) {
        uint64_t face_ind = cell.faces[f];
        BSPface& face = faces[face_ind];
        vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
        list_faceVertices(face, face_vrts);

        if (constraint_innerIntersects_face(face_vrts)) {
            // Here the face is splitted.
            // The intersection between face and constraint-plane is a new edge.
            splitFace(face_ind, constr, cell_ind, face_vrts);

            // Add new edge (the face-slpitting one) to cell_edges
            cell_edges.push_back(edges.size() - 1);
        }
    }

    // The cell is divided in two subcells: the up-subcell and the down-subcell.
    // The up-subcell have vertices with non-negative cell_vrts_or.
    // The down-subcell have vertices with non-positive cell_vrts_or.
    //
    // The cell faces are partitioned between the two subcells.
    // A cell face intesecting the constraint-plane is splitted in sub-faces.
    //
    // A new face is created: the common face between up-subcell and down-subcell.
    //
    // Conventionally the down-subcell replace cell in the vector cells,
    // while up-subcell is appended to the vector cells.

    // up-subcell
    cells.push_back(BSPcell());
    uint64_t newCell_ind = cells.size() - 1;

    facesPartition(cell_ind, newCell_ind, cell_vrts);
    // Add common face between up-subcell and down-subcell.
    add_commonFace(constr, cell_ind, newCell_ind, cell_vrts, cell_edges);

    // If there are coplanar constraints (to the one used for splitting) those
    // constraints have to be aded to the common face.
    uint32_t num_coplanar_constr = (uint32_t)coplanar_constr.size();
    faces.back().coplanar_constraints.resize(1 + num_coplanar_constr);
    faces.back().coplanar_constraints[0] = constr;
    if (num_coplanar_constr > 0)
        for (uint32_t cc = 0; cc < num_coplanar_constr; cc++)
            faces.back().coplanar_constraints[1 + cc] = coplanar_constr[cc];

    // Constraints that have to be partitioned between up-subcell
    // and down-subcell are: cells[cell].constarints.

    constraintsPartition(constr, cell_ind, newCell_ind, cell_vrts);
}

//--Complex tetrahedralization-------------------------

//  Input: the index of a BSPface w.r.t. vector faces: face_ind.
// Output: nothing.
// Note. The last 2 edges of the vector face[face_ind].edges are used to
//       detach a triangualr face from the face faces[face_ind].
//       This 2 last edges are replaced in the vector face[face_ind].edges by
//       one new edge: the one that closes the detached triangualr face.
void BSPcomplex::triangle_detach(uint64_t face_ind)
{
    BSPface& face = faces[face_ind];
    uint64_t num_face_edges = face.edges.size();

    // A triangle will be created by using:
    // - edges[face.edges.size()-1] and edges[face.edges.size()-2],
    // - introducing a new_edge to close the triangle.
    uint64_t s_01_ind = face.edges[num_face_edges - 1];
    //  BSPedge& s_01 = edges[s_01_ind];
    uint64_t s_12_ind = face.edges[num_face_edges - 2];
    //  BSPedge& s_12 = edges[s_12_ind];
    uint32_t t1 = consecEdges_common_endpt(
        edges[s_01_ind].vertices[0],
        edges[s_01_ind].vertices[1],
        edges[s_12_ind].vertices[0],
        edges[s_12_ind].vertices[1]);
    uint32_t t0 = other_edge_endpt(edges[s_01_ind].vertices[0], edges[s_01_ind].vertices[1], t1);
    uint32_t t2 = other_edge_endpt(edges[s_12_ind].vertices[0], edges[s_12_ind].vertices[1], t1);

    // Connect t2 and t0 with a new edge.
    edges.push_back(BSPedge());
    uint64_t s_20_ind = edges.size() - 1;
    BSPedge& s_20 = edges.back();
    s_20.vertices[0] = t2;
    s_20.vertices[1] = t0;
    // meshVertices -> all UINT32_MAX.
    s_20.meshVertices[0] = UINT32_MAX;
    s_20.meshVertices[1] = UINT32_MAX;
    s_20.meshVertices[2] = UINT32_MAX;
    s_20.meshVertices[3] = UINT32_MAX;
    s_20.meshVertices[4] = UINT32_MAX;
    s_20.meshVertices[5] = UINT32_MAX;
    // s_20.conn_faces.push_back(face_ind);

    // The other conn_faces element is the new triangle that is going to be
    // created below.

    // Add an element to edge_visit
    edge_visit.push_back(0);

    // Remove triangle <t0,t1,t2> from current BSPface and save it as a new one.
    face.edges.pop_back();
    face.edges[face.edges.size() - 1] = s_20_ind;
    uint64_t i = 0;
    // REMOVE_ELEM_VECT(face_ind, edges[s_01_ind].conn_faces);
    // REMOVE_ELEM_VECT(face_ind, edges[s_12_ind].conn_faces);

    // New face-triangle is created.
    faces.push_back(BSPface(
        face.meshVertices[0],
        face.meshVertices[1],
        face.meshVertices[2],
        face.conn_cells[0],
        face.conn_cells[1],
        face.colour));
    uint64_t new_face_ind = faces.size() - 1;
    BSPface& new_face = faces.back();
    new_face.edges.push_back(s_20_ind);
    new_face.edges.push_back(s_12_ind);
    new_face.edges.push_back(s_01_ind);
    // The sub-triangle lies on the same plane as its parent, so it inherits the same
    // coplanar input constraints (needed by face-provenance tracking). Empty for
    // WHITE (interior) faces, so this is free there.
    new_face.coplanar_constraints = faces[face_ind].coplanar_constraints;

    cells[faces[face_ind].conn_cells[0]].faces.push_back(new_face_ind);
    if (!IS_GHOST_CELL(faces[face_ind].conn_cells[1]))
        cells[faces[face_ind].conn_cells[1]].faces.push_back(new_face_ind);

    // edges[s_01_ind].conn_faces.push_back(new_face_ind);
    // edges[s_12_ind].conn_faces.push_back(new_face_ind);
    // s_20.conn_faces.push_back(new_face_ind);
    edges[s_01_ind].conn_face_0 = new_face_ind;
    edges[s_12_ind].conn_face_0 = new_face_ind;
    s_20.conn_face_0 = new_face_ind;

#ifdef DEBUG_BSP_DEEP
    printf("detached triangle %llu -> <%u,%u,%u>\n", new_face_ind, t0, t1, t2);
    triangular_BSPface_isDegenerate(faces, edges, vertices, new_face_ind);
#endif
}

//  Input:
// Output:
// Note. assumes that aligned face-edges are sub-edges of the same original edge.
bool BSPcomplex::aligned_face_edges(uint64_t fe0, uint64_t fe1, const BSPface& face)
{
    BSPedge& edge0 = edges[face.edges[fe0]];
    BSPedge& edge1 = edges[face.edges[fe1]];
    if (edge0.meshVertices[0] == UINT32_MAX) return false;
    if (edge0.meshVertices[0] != edge1.meshVertices[0]) return false;
    if (edge0.meshVertices[1] != edge1.meshVertices[1]) return false;
    if (edge0.meshVertices[2] != edge1.meshVertices[2]) return false;
    if (edge0.meshVertices[2] == UINT32_MAX) return true;
    if (edge0.meshVertices[3] != edge1.meshVertices[3]) return false;
    if (edge0.meshVertices[4] != edge1.meshVertices[4]) return false;
    if (edge0.meshVertices[5] != edge1.meshVertices[5]) return false;
    return true;
}

//
//
//  Input: the index of a BSPface: face_ind.
// Output: nothing.
// Do the unordered vertex triples {t[0],t[1],t[2]} and {x,y,z} describe the same
// triangle? (All three of t's vertices must be one of x,y,z.)
static inline bool
same_triangle(const std::array<uint32_t, 3>& t, uint32_t x, uint32_t y, uint32_t z)
{
    return (t[0] == x || t[0] == y || t[0] == z) && (t[1] == x || t[1] == y || t[1] == z) &&
        (t[2] == x || t[2] == y || t[2] == z);
}

// Triangulate the CONVEX polygon whose boundary vertices are `poly` (in boundary
// order, global vertex indices) into positive-area triangles, returned as vertex-
// index triples. `is_flat[i]` marks poly[i] as a "flat" vertex -- a Steiner point
// a split left in the middle of a straight boundary edge (as opposed to a "corner",
// where the boundary actually turns); the caller classifies these cheaply (see
// triangulateFace).
//
// BSP faces are convex but carry these flat vertices, so a naive fan/ear-clipping
// produces (or strands) zero-area triangles. This method is provably all-positive
// and, crucially, needs NO geometric predicate of its own -- it is pure integer
// bookkeeping:
//   1. Fan the corners from corners[0]: triangles (corners[0], corners[j],
//      corners[j+1]). Three corners of a convex polygon are never collinear, so
//      every fan triangle is positive.
//   2. Insert the flats. The flats between two consecutive corners c_a, c_b form a
//      straight boundary chain that is an edge of exactly one fan triangle. Walk
//      that run in boundary order; each flat f splits the triangle currently
//      carrying the edge (near, c_b) -- found by integer vertex membership -- as
//      (A,B,W) -> (A,f,W) + (f,B,W), where W is that triangle's apex (necessarily
//      off the chain's line, since the triangle is non-degenerate) and f lies
//      strictly between near and c_b on the line, so both pieces are positive.
//      `near` then advances to f. Reading W from the found triangle keeps this
//      correct even where two chains meet the same fan triangle (at corners[1] and
//      corners[m-1]): whichever chain splits first, the other simply finds the new,
//      smaller triangle carrying its edge.
// Developed and stress-tested in tests/triangulation_playground.py.
std::vector<std::array<uint32_t, 3>>
BSPcomplex::triangulateConvexFace(
    const std::vector<uint32_t>& poly, const std::vector<char>& is_flat)
{
    const uint32_t n = (uint32_t)poly.size();
    std::vector<std::array<uint32_t, 3>> tris;
    tris.reserve(n - 2);

    // Corner indices, in boundary order.
    tri_corner_list.clear();
    for (uint32_t i = 0; i < n; i++)
        if (!is_flat[i]) tri_corner_list.push_back(i);
    const size_t m = tri_corner_list.size();
    assert(m >= 3 && "triangulateConvexFace: fewer than 3 corners (degenerate face)");

    // 1) Fan the corners.
    for (size_t j = 1; j + 1 < m; j++)
        tris.push_back({poly[tri_corner_list[0]], poly[tri_corner_list[j]],
                        poly[tri_corner_list[j + 1]]});

    // 2) Insert the flats, run by run (chain between two consecutive corners).
    for (size_t k = 0; k < m; k++) {
        const uint32_t ca = tri_corner_list[k];
        const uint32_t cb = tri_corner_list[(k + 1) % m];
        if ((ca + 1) % n == cb) continue; // corners adjacent -> no flats on this chain
        const uint32_t cb_v = poly[cb];
        uint32_t near_v = poly[ca];
        for (uint32_t i = (ca + 1) % n; i != cb; i = (i + 1) % n) {
            const uint32_t f = poly[i];
            // Find the current triangle carrying edge (near_v, cb_v) by integer match.
            size_t ti = tris.size();
            int p = -1;
            for (size_t t = 0; t < tris.size() && p < 0; t++)
                for (int e = 0; e < 3; e++) {
                    const uint32_t a = tris[t][e], b = tris[t][(e + 1) % 3];
                    if ((a == near_v && b == cb_v) || (a == cb_v && b == near_v)) {
                        ti = t;
                        p = e;
                        break;
                    }
                }
            assert(p >= 0 && "triangulateConvexFace: flat-run edge not found");
            const uint32_t A = tris[ti][p];           // near_v or cb_v
            const uint32_t B = tris[ti][(p + 1) % 3]; // the other of the pair
            const uint32_t W = tris[ti][(p + 2) % 3]; // apex, off the chain line
            // Split edge A->B at f, preserving winding: (A,f,W) + (f,B,W).
            tris[ti] = {A, f, W};
            tris.push_back({f, B, W});
            near_v = f;
        }
    }
    return tris;
}

void BSPcomplex::triangulateFace(uint64_t face_ind)
{
    uint64_t num_face_edges = faces[face_ind].edges.size();
    if (num_face_edges <= 3) return;
    const uint32_t n = (uint32_t)num_face_edges;

    // Boundary vertices (in order) and the face's dominant projection plane.
    std::vector<uint32_t> poly(num_face_edges, UINT32_MAX);
    list_faceVertices(faces[face_ind], poly);
    const int n_max = face_dominant_normal_component(faces[face_ind]);

    // Classify each boundary vertex corner/flat. poly[i] is the shared endpoint of
    // the boundary edges at positions (i-1) and i, so it is a flat (a Steiner point
    // on a straight edge) when those two edges are aligned -- sub-edges of the same
    // original edge, which is how every Steiner point arises. aligned_face_edges is a
    // pure integer meshVertices comparison, so the common flats cost no predicate. We
    // only fall back to the exact orient2D when the edges are NOT aligned; that test
    // is non-degenerate (hence cheap) for a genuine corner and pays the expensive
    // exact path only for the rare collinear-but-not-aligned vertex. This is what
    // keeps triangulation off orient2D's exact-arithmetic fallback (its dominant cost).
    tri_is_flat.assign(n, 0);
    for (uint32_t i = 0; i < n; i++)
        if (aligned_face_edges((i + n - 1) % n, i, faces[face_ind]) ||
            genericPoint::orient2D(
                *vertices[poly[(i + n - 1) % n]], *vertices[poly[i]],
                *vertices[poly[(i + 1) % n]], n_max) == 0)
            tri_is_flat[i] = 1;

    // Robust, positive-area triangulation of the convex face (pure integer work).
    std::vector<std::array<uint32_t, 3>> tris = triangulateConvexFace(poly, tri_is_flat);

#ifndef NDEBUG
    // Assert the triangulation is a valid, all-positive cover of the face: exactly
    // n-2 triangles, every vertex used, and every triangle non-degenerate with the
    // same orientation. Catches a non-convex/degenerate face or a triangulation bug
    // right here rather than as a mysterious zero-volume tet downstream.
    assert(tris.size() == (size_t)num_face_edges - 2 && "triangulateFace: wrong triangle count");
    {
        const int s0 = genericPoint::orient2D(
            *vertices[tris[0][0]], *vertices[tris[0][1]], *vertices[tris[0][2]], n_max);
        std::set<uint32_t> used;
        for (const auto& t : tris) {
            const int s =
                genericPoint::orient2D(*vertices[t[0]], *vertices[t[1]], *vertices[t[2]], n_max);
            assert(s != 0 && (s > 0) == (s0 > 0) &&
                   "triangulateFace: degenerate or inconsistently-oriented triangle");
            used.insert(t[0]);
            used.insert(t[1]);
            used.insert(t[2]);
        }
        assert(used.size() == (size_t)num_face_edges && "triangulateFace: not all vertices used");
    }
#endif

    // ---- Realize `tris` on the BSP mesh via ear-clipping ------------------
    //
    // `tris` is the abstract triangulation (vertex-index triples); we now have to
    // materialize it in the mesh data structure. The only primitive available for
    // that is triangle_detach(face_ind), which "clips an ear": it takes the LAST
    // TWO edges of face.edges -- edges[m-2] and edges[m-1], which share the corner
    // vertex t1 and span the triangle <t0,t1,t2> -- removes that triangle from the
    // face as a brand-new triangular BSPface, and closes the wound with a fresh
    // diagonal edge t2-t0, shrinking this face by exactly one edge (m -> m-1) while
    // fixing up all edge/face/cell connectivity. triangle_detach can therefore only
    // ever cut the ear that currently sits in the last-two-edges slot.
    //
    // So we cannot pick which triangle to emit; we can only clip whatever ear the
    // last two edges currently form, and rotate the edge list to change which ear
    // that is. The loop is thus a match-or-rotate cycle:
    //
    //   * Read the current last-two-edges ear <t0,t1,t2>.
    //   * If that ear is one of the triangles we still owe (same_triangle scan),
    //     clip it: drop it from `tris` and call triangle_detach. The face loses an
    //     edge; one fewer triangle remains.
    //   * Otherwise rotate the edge list by one (down-shift) so a different pair of
    //     edges becomes the last two, and try again.
    //
    // This terminates because our `tris` IS a valid triangulation of this convex
    // polygon, and every polygon triangulation is an ear decomposition: at any
    // stage the not-yet-clipped triangles tile the current sub-polygon, and such a
    // tiling always contains at least two ears (triangles with two polygon-boundary
    // edges). An ear of the sub-polygon is exactly a triangle whose two boundary
    // edges are adjacent in face.edges, so at least one owed triangle can always be
    // rotated into the last-two-edges slot -- we never get stuck, and after n-3
    // clips the face is left as the final single triangle. guard_max bounds the
    // rotations between clips (< one full turn per clip) purely as a debug backstop.
    uint64_t guard = 0;
    const uint64_t guard_max = 2 * num_face_edges * num_face_edges + 16;
    while (num_face_edges > 3) {
        const uint64_t el = faces[face_ind].edges[num_face_edges - 1];
        const uint64_t ep = faces[face_ind].edges[num_face_edges - 2];
        const uint32_t t1 = consecEdges_common_endpt(
            edges[el].vertices[0], edges[el].vertices[1], edges[ep].vertices[0],
            edges[ep].vertices[1]);
        const uint32_t t0 = other_edge_endpt(edges[el].vertices[0], edges[el].vertices[1], t1);
        const uint32_t t2 = other_edge_endpt(edges[ep].vertices[0], edges[ep].vertices[1], t1);

        // Clip this ear iff it is one of our triangles; otherwise rotate to bring
        // another ear into the last-two-edges position.
        size_t match = tris.size();
        for (size_t i = 0; i < tris.size(); i++)
            if (same_triangle(tris[i], t0, t1, t2)) {
                match = i;
                break;
            }

        if (match < tris.size()) {
            tris[match] = tris.back();
            tris.pop_back();
            triangle_detach(face_ind);
            num_face_edges--;
        } else {
            UINT64_vect_down_shift(faces[face_ind].edges, 1);
        }

        if (++guard > guard_max) {
            assert(false && "triangulateFace: ear-decomposition did not terminate");
            break;
        }
    }
    assert(tris.size() == 1 && "triangulateFace: triangulation not fully realized");
}

// Append to `vertices` a point strictly interior to `cell` (its faces are already
// triangulated), to be used as the apex that fans the cell into tetrahedra.
//
// Fast path: the double average of the cell vertices. For a convex cell the exact
// vertex centroid is strictly interior, but its DOUBLE approximation can round
// exactly onto a face plane (e.g. a shallow "pyramid" cell whose apex is barely off
// an oblique base) -- then the fan tet on that face is degenerate (orient3D == 0) and
// cannot be repaired by the winding-flip below. So we check every fan tet with the
// exact orient3D predicate; if none is degenerate, keep the cheap explicit (double)
// point. (A merely wrong-signed tet is not a problem -- the flip fixes it.)
//
// Otherwise fall back to an EXACT barycenter: implicitPoint3D_TBC of four
// non-coplanar cell vertices. The centroid of four non-coplanar points lies in the
// open interior of their tetrahedron, which is contained in the open interior of the
// convex cell, so it is strictly interior to every face; and because TBC is an
// implicit (rational) point, orient3D evaluates it exactly, never rounding onto a
// plane.
genericPoint* BSPcomplex::computeBaricenterPoint(const vector<uint32_t>& vrts, const BSPcell& cell)
{
    // Exact interior barycenter: the barycenter of four affinely-independent (non-coplanar)
    // cell vertices. It is strictly inside the tetrahedron of those four vertices, which (the
    // cell being convex) lies inside the cell -- so it is strictly interior and always a valid
    // star-center (an indirect/implicit point). Built greedily: add a vertex only if
    // independent of the ones already chosen (2nd distinct, 3rd not collinear, 4th not
    // coplanar). A non-degenerate (positive-volume) cell always has such a quadruple; quad[]
    // never holds an out-of-range index, so the predicates only see already-selected points.
    const uint32_t n = (uint32_t)vrts.size();
    uint32_t quad[4];
    uint32_t nq = 0;
    for (uint32_t i = 0; i < n && nq < 4; i++) {
        const genericPoint& t = *vertices[vrts[i]];
        bool independent;
        if (nq < 2)
            independent = true; // 1st vertex, and any distinct 2nd vertex
        else if (nq == 2)
            independent = genericPoint::orient2Dxy(*vertices[quad[0]], *vertices[quad[1]], t) != 0 ||
                genericPoint::orient2Dyz(*vertices[quad[0]], *vertices[quad[1]], t) != 0 ||
                genericPoint::orient2Dzx(*vertices[quad[0]], *vertices[quad[1]], t) != 0;
        else // nq == 3
            independent = genericPoint::orient3D(
                              *vertices[quad[0]], *vertices[quad[1]], *vertices[quad[2]], t) != 0;
        if (independent) quad[nq++] = vrts[i];
    }
    assert(nq == 4 && "computeBaricenter: cell has no 4 non-coplanar vertices (flat cell)");
    implicitPoint3D_TBC* exact = new implicitPoint3D_TBC(
        *vertices[quad[0]], *vertices[quad[1]], *vertices[quad[2]], *vertices[quad[3]]);

#ifndef VOLUMEREMESHER_BARY_ALWAYS_EXACT
    // Cheap path: the approximate (double) centroid of all cell vertices. Accept it only if it
    // is STRICTLY INTERIOR -- on the same side as the exact interior barycenter for every cell
    // face (and off every face plane). This keeps the common case on a light explicit point;
    // only thin/near-flat cells, where the double centroid can round to just outside (its fan
    // tets would then overlap the neighbour), fall back to the exact point. Using `exact` as
    // the reference is one orient3D per face -- no scan over the cell vertices.
    double cx, cy, cz;
    double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    uint32_t np = 0;
    for (const uint32_t v : vrts)
        if (vertices[v]->getApproxXYZCoordinates(cx, cy, cz)) {
            sum_x += cx;
            sum_y += cy;
            sum_z += cz;
            np++;
        }
    explicitPoint3D* cand = new explicitPoint3D(sum_x / np, sum_y / np, sum_z / np);
    bool interior = true;
    for (uint64_t fi : cell.faces) {
        vector<uint32_t> fv(3, UINT32_MAX);
        list_faceVertices(faces[fi], fv);
        const genericPoint& fa = *vertices[fv[0]];
        const genericPoint& fb = *vertices[fv[1]];
        const genericPoint& fc = *vertices[fv[2]];
        const int sc = genericPoint::orient3D(fa, fb, fc, *cand);
        const int sr = genericPoint::orient3D(fa, fb, fc, *exact);
        if (sc == 0 || (sc < 0) != (sr < 0)) {
            interior = false;
            break;
        }
    }
    if (interior) {
        delete exact;
        return cand;
    }
    delete cand;
#endif

    // Exact barycenter (always with VOLUMEREMESHER_BARY_ALWAYS_EXACT; otherwise only for the
    // thin cells whose cheap centroid was rejected above).
    return exact;
}

// Thread-safe list_cellVertices: same result and same order as list_cellVertices, but using
// local dedup containers instead of the shared edge_visit/vrts_visit scratch, so it is safe to
// call concurrently for different cells.
void BSPcomplex::list_cellVertices_ts(const BSPcell& cell, vector<uint32_t>& cell_vrts) const
{
    std::unordered_set<uint64_t> edge_seen;
    std::vector<uint64_t> cell_edges;
    for (uint64_t f : cell.faces)
        for (uint64_t e : faces[f].edges)
            if (edge_seen.insert(e).second) cell_edges.push_back(e);

    cell_vrts.clear();
    std::unordered_set<uint32_t> vrt_seen;
    for (uint64_t e : cell_edges) {
        const uint32_t a = edges[e].vertices[0], b = edges[e].vertices[1];
        if (vrt_seen.insert(a).second) cell_vrts.push_back(a);
        if (vrt_seen.insert(b).second) cell_vrts.push_back(b);
    }
}

//
//
inline uint64_t BSPcomplex::triFace_oppEdge(const BSPface& face, uint32_t v)
{
    uint64_t edge_ind = face.edges[0];
    uint32_t e0 = edges[edge_ind].vertices[0];
    uint32_t e1 = edges[edge_ind].vertices[1];
    if (e0 == v || e1 == v) {
        edge_ind = face.edges[1];
        e0 = edges[edge_ind].vertices[0];
        e1 = edges[edge_ind].vertices[1];
        if (e0 == v || e1 == v) {
            edge_ind = face.edges[2];
            e0 = edges[edge_ind].vertices[0];
            e1 = edges[edge_ind].vertices[1];
        }
    }
    return edge_ind;
}

//
//
uint64_t
BSPcomplex::triFace_shareEdge(const BSPcell& cell, uint64_t face_ind, uint64_t vOppEdge_ind)
{
    for (uint64_t f = 0; f < cell.faces.size(); f++) {
        uint64_t adj_face_ind = cell.faces[f];
        if ((faces[adj_face_ind].edges[0] == vOppEdge_ind ||
             faces[adj_face_ind].edges[1] == vOppEdge_ind ||
             faces[adj_face_ind].edges[2] == vOppEdge_ind) &&
            face_ind != adj_face_ind)
            return adj_face_ind;
    }

    // never reached
    printf("[BSP.cpp]BSPcomplex::triFace_shareEdge: ERROR return UINT64_MAX.\n");
    return UINT64_MAX;
}

//
//
bool BSPcomplex::cell_is_tetrahedrizable_from_v(const BSPcell& cell, uint32_t v)
{
    uint64_t num_incFaces = count_cellFaces_inc_cellVrt(cell, v);
    vector<uint64_t> v_incFaces(num_incFaces, UINT64_MAX);
    cell_VFrelation(cell, v, v_incFaces);

    // bool return_zero = false;

    for (uint64_t f = 0; f < num_incFaces; f++) {
        // return_zero = false;
        BSPface& face = faces[v_incFaces[f]];
        uint64_t vOppEdge_ind = triFace_oppEdge(face, v);
        uint64_t faceShareEdge_ind = triFace_shareEdge(cell, v_incFaces[f], vOppEdge_ind);
        BSPface& oppFace = faces[faceShareEdge_ind];
        if (oppFace.meshVertices[0] == face.meshVertices[0] &&
            oppFace.meshVertices[1] == face.meshVertices[1] &&
            oppFace.meshVertices[2] == face.meshVertices[2]) // return_zero = true;
            return false;
    }

    return true;
}


//
//
// Run fn over blocks of [0,n) on std::thread::hardware_concurrency() threads, dynamically
// scheduled (atomic block-fetch) so uneven per-item cost stays balanced. Falls back to a
// serial call for small n or a single hardware thread. Compiling with VOLUMEREMESHER_SERIAL_TET
// (CMake: -DVOLUMEREMESHER_PARALLEL_TETRAHEDRALIZATION=OFF) forces the serial path everywhere.
//
// CAREFUL with exact arithmetic in here: bigrational, bigfloat and expansion all allocate
// from thread-local pools, so such a value must be born, used and destroyed inside ONE call
// of fn -- let only indices and other plain data escape. Computing them in parallel and
// reading them after the join is a use-after-free that compiles cleanly and asserts nothing.
// See the THREADING section of include/VolumeRemesher/numerics.h.
template <class F>
static void parallel_blocks(uint64_t n, F&& fn)
{
#ifdef VOLUMEREMESHER_SERIAL_TET
    fn(uint64_t(0), n);
#else
    unsigned nthreads = std::thread::hardware_concurrency();
    if (nthreads == 0) nthreads = 1;
    if (nthreads == 1 || n < 512) {
        fn(uint64_t(0), n);
        return;
    }
    std::atomic<uint64_t> next{0};
    const uint64_t block = 64;
    auto runner = [&]() {
        uint64_t i;
        while ((i = next.fetch_add(block)) < n) fn(i, std::min(n, i + block));
    };
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    for (unsigned t = 0; t < nthreads; t++) pool.emplace_back(runner);
    for (std::thread& t : pool) t.join();
#endif
}

void BSPcomplex::makeTetrahedra(bool verbose, bool keep_all_cells)
{
    uint64_t tet_num = 0; // total number of tetrahedra in which the cell will
                          // be decomposed.
    std::vector<uint32_t> decomposition_type(cells.size(), 0);
    // Possible decoposition techniques are:
    // 0 - cell is a tetrahedron -> no decomposition.
    // 1 - cell can be decomposed in tetrahedra by connecting a vertex with
    //     the other vertices that do not belong to its link.
    // 2 - cell is decomposed by connecting its baricenter with all the cell's
    //     vertices.
    std::vector<uint32_t> decomposition_vrt(cells.size(), UINT32_MAX);
    // decomposition_vrt is:
    // - UINT32_MAX for a tetrhedron,
    // - a cell vertex for decomposition type 1,
    // - the cell baricenter for decomposition type 2.

    // Per-cell decomposition decision. The heavy part -- the barycenter interior test in
    // computeBaricenterPoint -- lives here, and every cell is independent, so run it in
    // parallel: each thread writes only its own cells' slots and, for a barycenter (type-2)
    // cell, creates the star-center point WITHOUT adding it to `vertices`. The points are then
    // appended serially in cell order below, so vertex indices -- and the whole output --
    // match a serial run exactly. (list_cellVertices_ts avoids the shared visit scratch.)
    std::vector<uint64_t> tet_count(cells.size(), 0);
    std::vector<genericPoint*> bary_point(cells.size(), nullptr);

    parallel_blocks(cells.size(), [&](uint64_t lo, uint64_t hi) {
        vector<uint32_t> cell_vrts;
        for (uint64_t cell_i = lo; cell_i < hi; cell_i++) {
            BSPcell& cell = cells[cell_i];
            if (!keep_all_cells && cell.place != INTERNAL_A) continue;

            // If cell has more than 4 faces -> chose between types 1 and 2
            if (cell.faces.size() > 4) {
                list_cellVertices_ts(cell, cell_vrts);

                // Check if cell is tetrahedralizable from a vertex.
                bool needs_barycenter = true;
                for (uint32_t v : cell_vrts)
                    if (cell_is_tetrahedrizable_from_v(cell, v)) {
                        decomposition_type[cell_i] = 1;
                        decomposition_vrt[cell_i] = v;
                        tet_count[cell_i] = cell.faces.size() - count_cellFaces_inc_cellVrt(cell, v);
                        needs_barycenter = false;
                        break;
                    }

                if (needs_barycenter) { // Cell needs its baricenter
                    decomposition_type[cell_i] = 2;
                    bary_point[cell_i] = computeBaricenterPoint(cell_vrts, cell);
                    tet_count[cell_i] = cell.faces.size();
                }
            } else
                tet_count[cell_i] = 1; // The cell is a tet (type 0).
        }
    });

    // Serial merge: append the barycenter points in cell order (so their indices are identical
    // to a serial run) and total the tet count.
    for (uint64_t cell_i = 0; cell_i < cells.size(); cell_i++) {
        if (decomposition_type[cell_i] == 2) {
            vertices.push_back(bary_point[cell_i]);
            vrts_visit.push_back(0);
            decomposition_vrt[cell_i] = (uint32_t)vertices.size() - 1;
        }
        tet_num += tet_count[cell_i];
    }

    final_tets.reserve(tet_num * 4);
    final_tets_parent_cell.clear();
    final_tets_parent_cell.reserve(tet_num);

    // Make tets
    for (uint64_t cell_i = 0; cell_i < cells.size(); cell_i++) {
        BSPcell& cell = cells[cell_i];
        if (!keep_all_cells && cell.place != INTERNAL_A) continue;
        if (decomposition_type[cell_i] == 0) { // Simple tet
            vector<uint32_t> cell_vrts(4, UINT32_MAX);
            list_cellVertices(cells[cell_i], 6, cell_vrts);
            final_tets.insert(final_tets.end(), cell_vrts.begin(), cell_vrts.end());
            final_tets_parent_cell.push_back((uint32_t)cell_i);
        } else if (decomposition_type[cell_i] == 1) { // Tetrahedralizable from vertex
            uint32_t v = decomposition_vrt[cell_i];
            uint64_t num_incFaces = count_cellFaces_inc_cellVrt(cells[cell_i], v);
            uint64_t num_NOT_incFaces = cells[cell_i].faces.size() - num_incFaces;
            vector<uint64_t> v_NOT_incFaces(num_NOT_incFaces, UINT64_MAX);
            COMPL_cell_VFrelation(cells[cell_i], v, v_NOT_incFaces);
            for (uint64_t face_i : v_NOT_incFaces) {
                // Simple triangle
                vector<uint32_t> face_vrts(3, UINT32_MAX);
                list_faceVertices(faces[face_i], face_vrts);
                final_tets.insert(final_tets.end(), face_vrts.begin(), face_vrts.end());
                final_tets.push_back(v);
                final_tets_parent_cell.push_back((uint32_t)cell_i);
            }
        } else { // Uses cell barycenter
            for (uint64_t face_i : cells[cell_i].faces) {
                // Simple triangle
                vector<uint32_t> face_vrts(3, UINT32_MAX);
                list_faceVertices(faces[face_i], face_vrts);
                final_tets.insert(final_tets.end(), face_vrts.begin(), face_vrts.end());
                final_tets.push_back(decomposition_vrt[cell_i]);
                final_tets_parent_cell.push_back((uint32_t)cell_i);
            }
        }
    }

    // list_cellVertices / list_faceVertices do not fix a winding, so the tets come
    // out with arbitrary orientation. Flip each negatively-oriented one (swap two
    // vertices) so every tetrahedron has strictly positive volume. The sign is
    // decided with the exact orient3D predicate (no floating point). A genuinely
    // degenerate tet (orient3D == 0) cannot be fixed by a swap and is left as-is;
    // the triangulation above is meant to prevent those, and the debug-build assertion
    // below flags any that slip through rather than silently emitting them.
    // Each tet is independent (orient3D reads only, the swap touches only this tet's slots),
    // so flip in parallel.
    parallel_blocks(final_tets.size() / 4, [&](uint64_t lo, uint64_t hi) {
        for (uint64_t k = lo; k < hi; k++) {
            const uint32_t t = (uint32_t)(k * 4);
            if (genericPoint::orient3D(
                    *vertices[final_tets[t]],
                    *vertices[final_tets[t + 1]],
                    *vertices[final_tets[t + 2]],
                    *vertices[final_tets[t + 3]]) < 0)
                std::swap(final_tets[t + 2], final_tets[t + 3]);
        }
    });

#ifndef NDEBUG
    // Debug-only post-condition: EVERY emitted tet has strictly positive volume.
    //
    // The flip above trusts orient3D, so it cannot detect the one failure mode that
    // matters -- the predicate returning a sign that contradicts the exact coordinates of
    // its own vertices, which makes the flip turn a good tet into an inverted one. That is
    // not hypothetical: an implicit point whose homogeneous denominator has an undecided
    // sign made orient3D wrong on six tets of a 10.5M-tet arrangement, and it surfaced far
    // downstream as "a face appears twice in the tet list" rather than here.
    //
    // Re-asking orient3D would be circular, so this checks against exact rational
    // coordinates instead -- the same ground truth the caller would use. Debug-only: it
    // costs more than the tetrahedralization itself.
    //
    // The rationals stay inside one parallel task (thread-local pools -- see the THREADING
    // section of include/VolumeRemesher/numerics.h); only a flag comes out. The per-call
    // cache pays off because consecutive tets are one cell's fan and repeat its vertices.
    {
        std::atomic<uint64_t> bad{0};
        parallel_blocks(final_tets.size() / 4, [&](uint64_t lo, uint64_t hi) {
            std::unordered_map<uint32_t, std::array<bigrational, 3>> cache;
            const auto co = [&](uint32_t v) -> const std::array<bigrational, 3>& {
                auto it = cache.find(v);
                if (it != cache.end()) return it->second;
                std::array<bigrational, 3> c;
                vertices[v]->getExactXYZCoordinates(c[0], c[1], c[2]);
                return cache.emplace(v, std::move(c)).first->second;
            };
            for (uint64_t k = lo; k < hi; k++) {
                const uint32_t* t = &final_tets[4 * k];
                const std::array<bigrational, 3>&a = co(t[0]), &b = co(t[1]);
                const std::array<bigrational, 3>&c = co(t[2]), &d = co(t[3]);
                bigrational e1[3], e2[3], e3[3];
                for (int j = 0; j < 3; j++) {
                    e1[j] = b[j] - a[j];
                    e2[j] = c[j] - a[j];
                    e3[j] = d[j] - a[j];
                }
                const bigrational vol = (e1[1] * e2[2] - e1[2] * e2[1]) * e3[0] +
                    (e1[2] * e2[0] - e1[0] * e2[2]) * e3[1] +
                    (e1[0] * e2[1] - e1[1] * e2[0]) * e3[2];
                if (vol.sgn() <= 0) {
                    if (bad++ < 10)
                        printf(
                            "makeTetrahedra: tet %llu = [%u, %u, %u, %u] has %s volume\n",
                            (unsigned long long)k, t[0], t[1], t[2], t[3],
                            vol.sgn() < 0 ? "NEGATIVE" : "ZERO");
                }
            }
        });
        assert(bad == 0 && "makeTetrahedra emitted a non-positive tet");
    }
#endif

    // Tag output tet faces that lie on the input surface with their input triangles.
    // Done after the winding-flip so local-face indices match the emitted vertex order.
    trackFaceProvenance();
    // Tag output tet edges/vertices coming from inserted edges/points.
    trackEdgePointProvenance();

    if (verbose) printf("Tetrahedra: %lu\n", final_tets.size() / 4);
}

// Do the two coplanar triangles f0f1f2 and t0t1t2 overlap with positive area?
// Separating-axis test in the face's dominant plane (n_max) using the exact orient2D
// predicate: the triangles' interiors are disjoint iff some edge of one has the whole
// other triangle on its closed outer side. Edge/vertex-only contact => no overlap.
static bool coplanar_tris_overlap(
    const genericPoint& f0,
    const genericPoint& f1,
    const genericPoint& f2,
    const genericPoint& t0,
    const genericPoint& t1,
    const genericPoint& t2,
    int n_max)
{
    const genericPoint* F[3] = {&f0, &f1, &f2};
    const genericPoint* T[3] = {&t0, &t1, &t2};
    // For each triangle, test each of its edges as a separating axis against the other.
    for (int pass = 0; pass < 2; pass++) {
        const genericPoint** P = pass ? T : F; // edges from P
        const genericPoint** Q = pass ? F : T; // tested against Q
        for (int e = 0; e < 3; e++) {
            const genericPoint& a = *P[e];
            const genericPoint& b = *P[(e + 1) % 3];
            const genericPoint& c = *P[(e + 2) % 3]; // interior side of edge (a,b)
            const int in_side = genericPoint::orient2D(a, b, c, n_max);
            if (in_side == 0) continue; // degenerate edge, not a valid axis
            // Separating iff no Q vertex is strictly on P's interior side of (a,b).
            bool separating = true;
            for (int q = 0; q < 3 && separating; q++) {
                const int s = genericPoint::orient2D(a, b, *Q[q], n_max);
                if (s != 0 && (s > 0) == (in_side > 0)) separating = false;
            }
            if (separating) return false;
        }
    }
    return true;
}

// Partition the real input constraints into maximal groups that are transitively
// edge-adjacent AND coplanar. Two triangles sharing an edge are in the same group
// iff their four vertices are coplanar (exact orient3D). A flat region becomes one
// group; a triangle with no coplanar neighbour is its own singleton group.
void BSPcomplex::computeCoplanarGroups()
{
    // Group only genuine input triangles [0, first_fake_constraint). Hole-cap "fake"
    // constraints and virtual constraints carry no provenance and are left ungrouped.
    const uint32_t n = first_fake_constraint;
    std::vector<uint32_t> parent(n);
    for (uint32_t i = 0; i < n; i++) parent[i] = i;
    auto find = [&](uint32_t x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    auto other_v = [&](uint32_t c, uint32_t x, uint32_t y) -> uint32_t {
        for (int e = 0; e < 3; e++) {
            const uint32_t v = constraints_vrts[3 * c + e];
            if (v != x && v != y) return v;
        }
        return UINT32_MAX;
    };

    // Map each undirected constraint edge to the constraints that use it.
    std::map<std::pair<uint32_t, uint32_t>, std::vector<uint32_t>> edge_map;
    for (uint32_t c = 0; c < n; c++)
        for (int e = 0; e < 3; e++) {
            uint32_t v0 = constraints_vrts[3 * c + e];
            uint32_t v1 = constraints_vrts[3 * c + (e + 1) % 3];
            if (v0 > v1) std::swap(v0, v1);
            edge_map[{v0, v1}].push_back(c);
        }
    // Union constraints that share an edge and are coplanar.
    for (const auto& kv : edge_map) {
        const uint32_t v0 = kv.first.first, v1 = kv.first.second;
        const std::vector<uint32_t>& cs = kv.second;
        for (size_t i = 0; i < cs.size(); i++)
            for (size_t j = i + 1; j < cs.size(); j++)
                if (genericPoint::orient3D(*vertices[v0], *vertices[v1],
                        *vertices[other_v(cs[i], v0, v1)], *vertices[other_v(cs[j], v0, v1)]) == 0) {
                    const uint32_t ra = find(cs[i]), rb = find(cs[j]);
                    if (ra != rb) parent[ra] = rb;
                }
    }

    // Compact group roots to ids 0..K-1 and count group sizes.
    constraint_coplanar_group.assign(n, 0);
    coplanar_group_size.clear();
    std::map<uint32_t, uint32_t> root_id;
    for (uint32_t c = 0; c < n; c++) {
        const uint32_t r = find(c);
        auto it = root_id.find(r);
        uint32_t g;
        if (it == root_id.end()) {
            g = (uint32_t)root_id.size();
            root_id[r] = g;
            coplanar_group_size.push_back(0);
        } else
            g = it->second;
        constraint_coplanar_group[c] = g;
        coplanar_group_size[g]++;
    }
}

void BSPcomplex::trackFaceProvenance()
{
    computeCoplanarGroups();
    // One list of output faces per coplanar group (symmetric with edge/point provenance).
    triangle_provenance.assign(coplanar_group_size.size(), {});

    // Index every BLACK (on-surface) triangular face by its sorted vertex triple, so a
    // tet face can be matched to the BSP face it came from by exact integer identity.
    std::map<std::array<uint32_t, 3>, uint64_t> black_face;
    for (uint64_t fi = 0; fi < faces.size(); fi++) {
        const BSPface& face = faces[fi];
        if (face.colour == WHITE || face.edges.size() != 3) continue;
        std::vector<uint32_t> fv(3, UINT32_MAX);
        list_faceVertices(faces[fi], fv);
        std::array<uint32_t, 3> key = {fv[0], fv[1], fv[2]};
        std::sort(key.begin(), key.end());
        black_face[key] = fi;
    }
    if (black_face.empty()) return;

    // Cache, per BLACK face, the coplanar group(s) it overlaps (positive-area SAT over
    // its coplanar constraints). Usually one group; more than one only where two
    // exactly-coplanar surfaces overlap on the same plane (all triangles overlapped by
    // one face are on that face's single plane, but may belong to distinct groups).
    std::vector<int> cache_done(faces.size(), 0);
    std::vector<std::vector<uint32_t>> cache_groups(faces.size());
    // Per group: the distinct output faces already recorded (a surface face borders two kept
    // cells, so it is met from two tets; list it once).
    std::vector<std::set<std::array<uint32_t, 3>>> seen(coplanar_group_size.size());

    for (uint32_t k = 0; 4u * k < final_tets.size(); k++) {
        const uint32_t* tet = &final_tets[4 * k];
        for (int lf = 0; lf < 4; lf++) {
            // Face lf is opposite local vertex lf.
            const uint32_t a = tet[(lf + 1) & 3], b = tet[(lf + 2) & 3], c = tet[(lf + 3) & 3];
            std::array<uint32_t, 3> key = {a, b, c};
            std::sort(key.begin(), key.end());
            auto it = black_face.find(key);
            if (it == black_face.end()) continue;
            const uint64_t fi = it->second;

            if (!cache_done[fi]) {
                cache_done[fi] = 1;
                std::vector<uint32_t>& gs = cache_groups[fi];
                const BSPface& face = faces[fi];
                const int n_max = face_dominant_normal_component(face);
                for (uint32_t constr : face.coplanar_constraints) {
                    // Skip hole-cap (fake) and virtual constraints: only genuine input
                    // triangles [0, first_fake_constraint) contribute provenance, so a
                    // face lying only on caps stays untagged and is not reported.
                    if (constr >= first_fake_constraint) continue;
                    const uint32_t cID = 3 * constr;
                    if (coplanar_tris_overlap(
                            *vertices[a], *vertices[b], *vertices[c],
                            *vertices[constraints_vrts[cID]], *vertices[constraints_vrts[cID + 1]],
                            *vertices[constraints_vrts[cID + 2]], n_max)) {
                        const uint32_t g = constraint_coplanar_group[constr];
                        if (std::find(gs.begin(), gs.end(), g) == gs.end()) gs.push_back(g);
                    }
                }
                std::sort(gs.begin(), gs.end()); // deterministic output order
            }
            // Record this output face (tet k + its three vertices) under every group it
            // tiles, deduped by geometry (a surface face is met from its two adjacent tets;
            // keep one representative tet).
            for (uint32_t g : cache_groups[fi])
                if (seen[g].insert(key).second)
                    triangle_provenance[g].push_back({k, a, b, c});
        }
    }
}

void BSPcomplex::trackEdgePointProvenance()
{
    edge_provenance.clear();
    point_provenance.clear();
    const uint32_t n_edges = num_edge_triangles / 2;
    const uint32_t n_points = num_point_triangles / 3;
    if (n_edges == 0 && n_points == 0) return;

    // Which vertices are used by the output tets, and, for each output tet edge, one tet
    // that contains it (so an edge can be reported as {tet_id, v0, v1}).
    std::vector<char> used(vertices.size(), 0);
    for (uint32_t v : final_tets) used[v] = 1;
    std::map<std::pair<uint32_t, uint32_t>, uint32_t> tet_edges; // sorted pair -> a tet id
    std::vector<uint32_t> vertex_to_tet(vertices.size(), UINT32_MAX); // a tet containing v
    static const int E[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
    for (uint32_t k = 0; 4u * k < final_tets.size(); k++) {
        const uint32_t* tet = &final_tets[4 * k];
        for (int i = 0; i < 4; i++)
            if (vertex_to_tet[tet[i]] == UINT32_MAX) vertex_to_tet[tet[i]] = k;
        for (auto& ee : E) {
            uint32_t a = tet[ee[0]], b = tet[ee[1]];
            if (a > b) std::swap(a, b);
            tet_edges.emplace(std::make_pair(a, b), k); // keep the first tet seen
        }
    }

    // Each inserted edge is pinned by its two forcing triangles (A,B,*): A and B are the
    // first two vertices of the constraint at first_edge_constraint + 2*e. The output tet
    // edges on segment AB are the consecutive pairs of the output vertices lying on AB.
    //
    // Testing every output vertex against every edge is O(n_edges * n_vertices) and far too
    // slow at scale, so bucket the used output vertices into a uniform spatial grid (on
    // APPROXIMATE coordinates) and, per edge, gather only the vertices near the segment.
    // The grid merely prunes candidates; membership is still decided by the exact
    // pointInSegment test, so the approximate bucketing does not affect the result.
    if (n_edges > 0) {
        const uint32_t NV = (uint32_t)vertices.size();
        std::vector<double> ax(3 * NV, 0.0);
        double blo[3] = {DBL_MAX, DBL_MAX, DBL_MAX}, bhi[3] = {-DBL_MAX, -DBL_MAX, -DBL_MAX};
        uint64_t nused = 0;
        for (uint32_t v = 0; v < NV; v++)
            if (used[v]) {
                vertices[v]->getApproxXYZCoordinates(ax[3 * v], ax[3 * v + 1], ax[3 * v + 2]);
                for (int k = 0; k < 3; k++) {
                    if (ax[3 * v + k] < blo[k]) blo[k] = ax[3 * v + k];
                    if (ax[3 * v + k] > bhi[k]) bhi[k] = ax[3 * v + k];
                }
                nused++;
            }
        int gridN = (int)std::cbrt((double)(nused > 1 ? nused : 1));
        if (gridN < 4) gridN = 4;
        if (gridN > 256) gridN = 256;
        double ext = 0.0;
        for (int k = 0; k < 3; k++) ext = std::max(ext, bhi[k] - blo[k]);
        const double cell = (ext > 0.0) ? ext / gridN : 1.0;
        const int64_t OFF = 1 << 20; // keep packed cell indices non-negative
        auto ckey = [&](double x, double y, double z) -> uint64_t {
            int64_t ix = (int64_t)std::floor((x - blo[0]) / cell) + OFF;
            int64_t iy = (int64_t)std::floor((y - blo[1]) / cell) + OFF;
            int64_t iz = (int64_t)std::floor((z - blo[2]) / cell) + OFF;
            return (uint64_t)(ix & 0x1FFFFF) | ((uint64_t)(iy & 0x1FFFFF) << 21) |
                   ((uint64_t)(iz & 0x1FFFFF) << 42);
        };
        std::unordered_map<uint64_t, std::vector<uint32_t>> grid;
        for (uint32_t v = 0; v < NV; v++)
            if (used[v]) grid[ckey(ax[3 * v], ax[3 * v + 1], ax[3 * v + 2])].push_back(v);

        std::vector<uint32_t> stamp(NV, 0);
        uint32_t cur = 0;
        std::vector<uint32_t> cand;
        for (uint32_t e = 0; e < n_edges; e++) {
            const uint32_t c = first_edge_constraint + 2 * e;
            const uint32_t ia = constraints_vrts[3 * c], ib = constraints_vrts[3 * c + 1];
            double A[3], B[3];
            vertices[ia]->getApproxXYZCoordinates(A[0], A[1], A[2]);
            vertices[ib]->getApproxXYZCoordinates(B[0], B[1], B[2]);
            const double seglen = std::sqrt(
                (B[0] - A[0]) * (B[0] - A[0]) + (B[1] - A[1]) * (B[1] - A[1]) +
                (B[2] - A[2]) * (B[2] - A[2]));
            const int nsamp = std::max(1, (int)(seglen / cell) + 1);
            cur++;
            cand.clear();
            for (int s = 0; s <= nsamp; s++) {
                const double t = (double)s / (double)nsamp;
                const double px = A[0] + t * (B[0] - A[0]);
                const double py = A[1] + t * (B[1] - A[1]);
                const double pz = A[2] + t * (B[2] - A[2]);
                for (int dz = -1; dz <= 1; dz++)
                    for (int dy = -1; dy <= 1; dy++)
                        for (int dx = -1; dx <= 1; dx++) {
                            auto it = grid.find(
                                ckey(px + dx * cell, py + dy * cell, pz + dz * cell));
                            if (it == grid.end()) continue;
                            for (uint32_t v : it->second)
                                if (stamp[v] != cur) {
                                    stamp[v] = cur;
                                    cand.push_back(v);
                                }
                        }
            }
            std::vector<uint32_t> on; // candidates exactly on the closed segment [A,B]
            for (uint32_t v : cand)
                if (genericPoint::pointInSegment(*vertices[v], *vertices[ia], *vertices[ib]))
                    on.push_back(v);
            // Order them along AB: for collinear points the full lexicographic order
            // (lessThan) is monotonic along the line.
            std::sort(on.begin(), on.end(), [&](uint32_t x, uint32_t y) {
                return genericPoint::lessThan(*vertices[x], *vertices[y]) < 0;
            });
            std::vector<std::array<uint32_t, 3>> segs; // {tet_id, v0, v1}
            for (size_t i = 0; i + 1 < on.size(); i++) {
                uint32_t a = on[i], b = on[i + 1];
                uint32_t lo = a < b ? a : b, hi = a < b ? b : a;
                auto it = tet_edges.find({lo, hi});
                if (it != tet_edges.end()) segs.push_back({it->second, a, b});
            }
            edge_provenance.push_back({e, segs});
        }
    }

    // Each inserted point is pinned by its three corner triangles (P,*,*): P is the first
    // vertex of the constraint at first_point_constraint + 3*p. Report the output vertex
    // equal to it (P itself if it survived, else any coincident used vertex).
    for (uint32_t p = 0; p < n_points; p++) {
        const uint32_t c = first_point_constraint + 3 * p;
        const uint32_t ip = constraints_vrts[3 * c];
        uint32_t out = UINT32_MAX;
        if (used[ip])
            out = ip;
        else
            for (uint32_t v = 0; v < vertices.size(); v++)
                if (used[v] && genericPoint::coincident(*vertices[v], *vertices[ip])) {
                    out = v;
                    break;
                }
        const uint32_t tet = (out == UINT32_MAX) ? UINT32_MAX : vertex_to_tet[out];
        point_provenance.push_back({p, tet, out});
    }
}


//-Decide colour of GREY faces--------------------------------------------------

//
//
int BSPcomplex::face_dominant_normal_component(const BSPface& face)
{
    const uint32_t* mv = face.meshVertices;
    double mvc[9]; // Coords of the original input triangle
    vertices[mv[0]]->getApproxXYZCoordinates(mvc[0], mvc[1], mvc[2]);
    vertices[mv[1]]->getApproxXYZCoordinates(mvc[3], mvc[4], mvc[5]);
    vertices[mv[2]]->getApproxXYZCoordinates(mvc[6], mvc[7], mvc[8]);
    return genericPoint::maxComponentInTriangleNormal(
        mvc[0],
        mvc[1],
        mvc[2],
        mvc[3],
        mvc[4],
        mvc[5],
        mvc[6],
        mvc[7],
        mvc[8]);
}

//
//
void BSPcomplex::get_approx_faceBaricenterCoord(const BSPface& face, double* bar)
{
    bar[0] = 0;
    bar[1] = 0;
    bar[2] = 0;
    double tp[3];
    const BSPedge& edge0 = edges[face.edges.back()];
    const BSPedge& edge1 = edges[face.edges[0]];
    uint32_t vid = consecEdges_common_endpt(
        edge0.vertices[0],
        edge0.vertices[1],
        edge1.vertices[0],
        edge1.vertices[1]);
    for (uint64_t e = 0; e < face.edges.size(); e++) {
        const BSPedge& edge = edges[face.edges[e]];
        if (vid == edge.vertices[0])
            vid = edge.vertices[1];
        else
            vid = edge.vertices[0];
        vertices[vid]->getApproxXYZCoordinates(tp[0], tp[1], tp[2]);
        bar[0] += tp[0];
        bar[1] += tp[1];
        bar[2] += tp[2];
    }
    bar[0] /= face.edges.size();
    bar[1] /= face.edges.size();
    bar[2] /= face.edges.size();
}

//
//
bool BSPcomplex::is_baricenter_inFace(
    const BSPface& face,
    const explicitPoint3D& face_center,
    int max_normComp)
{
    const BSPedge& edge0 = edges[face.edges.back()];
    const BSPedge& edge1 = edges[face.edges[0]];
    uint32_t vid = consecEdges_common_endpt(
        edge0.vertices[0],
        edge0.vertices[1],
        edge1.vertices[0],
        edge1.vertices[1]);

    int oro = 0;
    uint64_t e;
    for (e = 0; e < face.edges.size(); e++) {
        const BSPedge& edge = edges[face.edges[e]];
        const uint32_t pvid = vid;
        if (vid == edge.vertices[0])
            vid = edge.vertices[1];
        else
            vid = edge.vertices[0];

        const int ao =
            genericPoint::orient2D(face_center, *vertices[pvid], *vertices[vid], max_normComp);
        if (ao == 0) break;
        if (ao != oro) {
            if (oro)
                break;
            else
                oro = ao;
        }
    }

    if (e == face.edges.size()) return true;
    return false;
}

//
//
inline bool faceColour_matches_constrGroup(CONSTR_GROUP_T group, bool f_blackA, bool f_blackB)
{
    return (group == CONSTR_A && f_blackA) || (group == CONSTR_B && f_blackB);
}

//
//
COLOUR_T BSPcomplex::blackAB_or_white(uint64_t face_ind, bool two_input)
{
    const BSPface& face = faces[face_ind];

    // Get dominant normal component
    int xyz = face_dominant_normal_component(face);

    // Calculate approximated face barycenter
    double p[3];
    get_approx_faceBaricenterCoord(face, p);
    const explicitPoint3D face_center(p[0], p[1], p[2]);

    // Needed only for two input case.
    bool face_is_blackA = false;
    bool face_is_blackB = false;

    // Check whether the barycenter is indeed inside the face (might be not due to approximation)
    if (is_baricenter_inFace(face, face_center, xyz)) // Barycenter is inside the face: just check
                                                      // that it is in one of the constraints too
    {
        for (uint32_t c = 0; c < face.coplanar_constraints.size(); c++) {
            const uint32_t constr = face.coplanar_constraints[c];
            const CONSTR_GROUP_T c_group = constraint_group[constr];

            if (two_input &&
                faceColour_matches_constrGroup(c_group, face_is_blackA, face_is_blackB))
                continue;

            const uint32_t constr_ID = 3 * constr;
            const genericPoint* c0 = vertices[constraints_vrts[constr_ID]];
            const genericPoint* c1 = vertices[constraints_vrts[constr_ID + 1]];
            const genericPoint* c2 = vertices[constraints_vrts[constr_ID + 2]];

            if (genericPoint::pointInTriangle(face_center, *c0, *c1, *c2, xyz)) {
                if (!two_input) return BLACK_A;

                if (c_group == CONSTR_A)
                    face_is_blackA = true;
                else if (c_group == CONSTR_B)
                    face_is_blackB = true;
                if (face_is_blackA && face_is_blackB) return BLACK_AB;
            }
        }

        if (two_input) {
            if (face_is_blackA)
                return BLACK_A;
            else if (face_is_blackB)
                return BLACK_B;
        }

        return WHITE;
    } else // Barycenter is not inside the face: revert to slow version
    {
        const BSPedge& edge0 = edges[face.edges.back()];
        const BSPedge& edge1 = edges[face.edges[0]];
        uint32_t vid = consecEdges_common_endpt(
            edge0.vertices[0],
            edge0.vertices[1],
            edge1.vertices[0],
            edge1.vertices[1]);

        for (uint64_t e = 0; e < face.edges.size(); e++) {
            const BSPedge& edge = edges[face.edges[e]];
            if (vid == edge.vertices[0])
                vid = edge.vertices[1];
            else
                vid = edge.vertices[0];

            uint32_t out_from_all = 0;
            const genericPoint* face_pt = vertices[vid];
            for (uint32_t c = 0; c < face.coplanar_constraints.size(); c++) {
                const uint32_t constr = face.coplanar_constraints[c];
                const CONSTR_GROUP_T c_group = constraint_group[constr];

                if (two_input &&
                    faceColour_matches_constrGroup(c_group, face_is_blackA, face_is_blackB))
                    continue;

                const uint32_t constr_ID = 3 * constr;
                const uint32_t vid1 = constraints_vrts[constr_ID];
                const uint32_t vid2 = constraints_vrts[constr_ID + 1];
                const uint32_t vid3 = constraints_vrts[constr_ID + 2];
                const genericPoint* c0 = vertices[vid1];
                const genericPoint* c1 = vertices[vid2];
                const genericPoint* c2 = vertices[vid3];

                if (vid == vid1 || vid == vid2 || vid == vid3) break; // point on boundary.
                int lpt = localizedPointInTriangle(*face_pt, *c0, *c1, *c2, xyz);
                // lpt = 1 -> point on boundary, lpt = 2 -> point in interior, otherwise lpt = 0.

                if (lpt == 2) {
                    if (!two_input) return BLACK_A;

                    if (c_group == CONSTR_A)
                        face_is_blackA = true;
                    else if (c_group == CONSTR_B)
                        face_is_blackB = true;
                    if (face_is_blackA && face_is_blackB) return BLACK_AB;
                }

                if (lpt) break;

                out_from_all++;
            }
            if (out_from_all == face.coplanar_constraints.size()) return WHITE;
        }

        // All face vertices are on the boundary of some coplanar constraints.
        for (uint32_t c = 0; c < face.coplanar_constraints.size(); c++) {
            const uint32_t constr = face.coplanar_constraints[c];
            const CONSTR_GROUP_T c_group = constraint_group[constr];
            if (two_input &&
                faceColour_matches_constrGroup(c_group, face_is_blackA, face_is_blackB))
                continue;

            const uint32_t* constraint = constraints_vrts.data() + 3 * constr;
            if (coplanar_constraint_innerIntersects_face(face.edges, constraint, xyz)) {
                if (!two_input) return BLACK_A;

                if (c_group == CONSTR_A)
                    face_is_blackA = true;
                else if (c_group == CONSTR_B)
                    face_is_blackB = true;
                if (face_is_blackA && face_is_blackB) return BLACK_AB;
            }
        }

        if (face_is_blackA)
            return BLACK_A;
        else if (face_is_blackB)
            return BLACK_B;

        return WHITE;
    }
}


//----------------
// BSP subdivision
//----------------

//
//

void BSPcomplex::extractSkinTriMesh(
    const char* filename,
    const char bool_opcode,
    double** coords,
    uint32_t* npts,
    uint32_t** tri_idx,
    uint32_t* ntri)
{
    const uint64_t num_faces = faces.size();
    for (uint64_t face_ind = 0; face_ind < num_faces; face_ind++) triangulateFace(face_ind);

    if (bool_opcode == 'U') { // Union
        for (BSPcell& cell : cells)
            cell.place =
                (cell.place == INTERNAL_A || cell.place == INTERNAL_B || cell.place == INTERNAL_AB)
                    ? (INTERNAL_A)
                    : (EXTERNAL);
    }

    else if (bool_opcode == 'I') { // Intersection
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_AB) ? (INTERNAL_A) : (EXTERNAL);
    }

    else if (bool_opcode == 'D') { // Difference A\B
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_A && !(cell.place == INTERNAL_AB)) ? (INTERNAL_A)
                                                                                    : (EXTERNAL);
    }

    // Set "internal" depending on bool_opcode and find border faces to save
    vector<uint64_t> mark(faces.size(), 0);
    for (BSPcell& cell : cells)
        if (cell.place == INTERNAL_A)
            for (uint64_t fi = 0; fi < cell.faces.size(); fi++) mark[cell.faces[fi]]++;

    for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
    for (size_t i = 0; i < edges.size(); i++) edge_visit[i] = 0;

    uint64_t num_border_faces = 0;
    for (uint64_t f_i = 0; f_i < faces.size(); f_i++)
        if (mark[f_i] == 1) {
            num_border_faces++;
            for (uint64_t eid : faces[f_i].edges) edge_visit[eid] = 1;
        }
    for (size_t i = 0; i < edges.size(); i++)
        if (edge_visit[i]) vrts_visit[edges[i].vertices[0]] = vrts_visit[edges[i].vertices[1]] = 1;

    std::vector<uint32_t> vmap(vertices.size(), 0);
    size_t num_v = 0;
    for (size_t i = 0; i < vrts_visit.size(); i++) {
        vmap[i] = (uint32_t)num_v;
        if (vrts_visit[i]) num_v++;
    }

    *coords = (double*)malloc(sizeof(double) * 3 * num_v);
    *tri_idx = (uint32_t*)malloc(sizeof(uint32_t) * 3 * num_border_faces);
    *npts = (uint32_t)num_v;
    *ntri = (uint32_t)num_border_faces;

    // Store vertex coordinates
    for (uint32_t v = 0, i = 0; v < vertices.size(); v++)
        if (vrts_visit[v]) {
            vertices[v]
                ->getApproxXYZCoordinates((*coords)[i], (*coords)[i + 1], (*coords)[i + 2], true);
            i += 3;
        }

    // Print border faces
    for (uint32_t f_i = 0, i = 0; f_i < faces.size(); f_i++)
        if (mark[f_i] == 1) {
            BSPface& face = faces[f_i];
            vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
            list_faceVertices(face, face_vrts);

            if (cells[face.conn_cells[0]].place == INTERNAL_A)
                for (uint32_t v = (uint32_t)face_vrts.size(); v > 0; v--)
                    (*tri_idx)[i++] = vmap[face_vrts[v - 1]];
            else
                for (uint32_t v = 0; v < face_vrts.size(); v++)
                    (*tri_idx)[i++] = vmap[face_vrts[v]];
        }
}


void BSPcomplex::saveSkin(const char* filename, const char bool_opcode, bool triangulate)
{
    if (triangulate) {
        const uint64_t num_faces = faces.size();
        for (uint64_t face_ind = 0; face_ind < num_faces; face_ind++) triangulateFace(face_ind);
    }

    if (bool_opcode == 'U') { // Union
        for (BSPcell& cell : cells)
            cell.place =
                (cell.place == INTERNAL_A || cell.place == INTERNAL_B || cell.place == INTERNAL_AB)
                    ? (INTERNAL_A)
                    : (EXTERNAL);
    }

    else if (bool_opcode == 'I') { // Intersection
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_AB) ? (INTERNAL_A) : (EXTERNAL);
    }

    else if (bool_opcode == 'D') { // Difference A\B
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_A && !(cell.place == INTERNAL_AB)) ? (INTERNAL_A)
                                                                                    : (EXTERNAL);
    }

    // Set "internal" depending on bool_opcode and find border faces to save
    vector<uint64_t> mark(faces.size(), 0);
    for (BSPcell& cell : cells)
        if (cell.place == INTERNAL_A)
            for (uint64_t fi = 0; fi < cell.faces.size(); fi++) mark[cell.faces[fi]]++;

    for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
    for (size_t i = 0; i < edges.size(); i++) edge_visit[i] = 0;

    uint64_t num_border_faces = 0;
    for (uint64_t f_i = 0; f_i < faces.size(); f_i++)
        if (mark[f_i] == 1) {
            num_border_faces++;
            for (uint64_t eid : faces[f_i].edges) edge_visit[eid] = 1;
        }
    for (size_t i = 0; i < edges.size(); i++)
        if (edge_visit[i]) vrts_visit[edges[i].vertices[0]] = vrts_visit[edges[i].vertices[1]] = 1;

    std::vector<uint32_t> vmap(vertices.size(), 0);
    size_t num_v = 0;
    for (size_t i = 0; i < vrts_visit.size(); i++) {
        vmap[i] = (uint32_t)num_v;
        if (vrts_visit[i]) num_v++;
    }

    // Binary mode so line endings stay LF on every OS (byte-identical output).
    ofstream f(filename, std::ios::binary);

    if (!f) ip_error("BSPcomplex::saveSkin: cannot open the file.\n");

    f << "OFF\n";
    f << num_v << " ";
    f << num_border_faces << " ";
    f << "0\n";

    // Print vertices coordinates
    for (uint32_t v = 0; v < vertices.size(); v++)
        if (vrts_visit[v]) f << (*(vertices)[v]) << "\n";

    // Print border faces
    for (uint32_t f_i = 0; f_i < faces.size(); f_i++)
        if (mark[f_i] == 1) {
            BSPface& face = faces[f_i];
            vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
            list_faceVertices(face, face_vrts);
            f << face_vrts.size();

            if (cells[face.conn_cells[0]].place == INTERNAL_A)
                for (uint32_t v = (uint32_t)face_vrts.size(); v > 0; v--)
                    f << " " << vmap[face_vrts[v - 1]];
            else
                for (uint32_t v = 0; v < face_vrts.size(); v++) f << " " << vmap[face_vrts[v]];

            f << "\n";
        }

    f.close();
}


//
//
// Portable decimal "[-]num[/den]" of an exact rational, independent of the bignum
// backend (gmpxx's get_str vs the in-house get_dec_str both emit base-10 num/den).
static inline std::string rational_to_string(const bigrational& r)
{
#ifdef USE_GNU_GMP_CLASSES
    return r.get_str();
#else
    return r.get_dec_str();
#endif
}

void BSPcomplex::saveMesh(
    const char* filename, const char bool_opcode, bool tetrahedrize, bool export_rational)
{
    // Binary mode so line endings stay LF on every OS (byte-identical output).
    ofstream f(filename, std::ios::binary);

    if (!f) ip_error("\nBSPcomplex::[BSP.cpp]saveTetMesh: FATAL ERROR cannot open the file.\n");

    // Full double precision: the default ostream precision (6 significant figures)
    // collapses/inverts elements on fine or large-coordinate meshes.
    f.precision(17);

    const uint64_t num_faces = faces.size();

    if (tetrahedrize)
        for (uint64_t face_ind = 0; face_ind < num_faces; face_ind++) triangulateFace(face_ind);

    if (bool_opcode == 'U') { // Union
        for (BSPcell& cell : cells)
            cell.place =
                (cell.place == INTERNAL_A || cell.place == INTERNAL_B || cell.place == INTERNAL_AB)
                    ? (INTERNAL_A)
                    : (EXTERNAL);
    }

    else if (bool_opcode == 'I') { // Intersection
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_AB) ? (INTERNAL_A) : (EXTERNAL);
    }

    else if (bool_opcode == 'D') { // Difference A\B
        for (BSPcell& cell : cells)
            cell.place = (cell.place == INTERNAL_A && !(cell.place == INTERNAL_AB)) ? (INTERNAL_A)
                                                                                    : (EXTERNAL);
    }

    if (tetrahedrize) {
        makeTetrahedra();

        for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
        for (uint32_t t = 0; t < final_tets.size(); t++) vrts_visit[final_tets[t]] = 1;
        uint32_t final_numver = 0;
        for (size_t i = 0; i < vrts_visit.size(); i++)
            if (vrts_visit[i]) final_numver++;

        std::vector<uint32_t> vmap(vertices.size(), 0);
        size_t num_v = 0;
        for (size_t i = 0; i < vrts_visit.size(); i++) {
            vmap[i] = (uint32_t)num_v;
            if (vrts_visit[i]) num_v++;
        }

        f << final_numver << " vertices\n";
        f << final_tets.size() / 4 << " tets\n";

        // Print vertices coordinates.
        // Use operator<< (approximate double coordinates), NOT get_str(): get_str()
        // serializes the exact rational in a base that depends on the bignum backend
        // (decimal with gmpxx, binary without), so its bytes diverge across platforms
        // -- e.g. MSVC (no gmpxx) vs Linux/macOS. operator<< prints portable doubles,
        // matching the .msh path below and read_TET_file's %lf parsing.
        for (uint32_t v = 0; v < vertices.size(); v++)
            if (vrts_visit[v]) f << (*vertices[v]) << "\n";

        // Print tets
        for (uint32_t t = 0; t < final_tets.size(); t += 4)
            f << "4 " << vmap[final_tets[t]] << " " << vmap[final_tets[t + 1]] << " "
              << vmap[final_tets[t + 2]] << " " << vmap[final_tets[t + 3]] << "\n";

        // Sidecar 1: input-surface provenance, keyed by coplanar group (symmetric with the
        // edge/point sidecars below). One line per group listing the output tet faces tiling
        // it: "<group_id> <num_faces> <tet v0 v1 v2> <tet v0 v1 v2> ...", where tet is the
        // output tet index (0-based, matching the tet block above) and v0 v1 v2 are its
        // face's output vertex indices. (A face on two overlapping coplanar surfaces is
        // listed under both groups; group sizes are in the .groups sidecar.)
        {
            ofstream pf((std::string(filename) + ".triangleprov").c_str(), std::ios::binary);
            pf << "# group_id num_faces tet v0 v1 v2 ...\n";
            pf << triangle_provenance.size() << "\n";
            for (uint32_t g = 0; g < triangle_provenance.size(); g++) {
                pf << g << " " << triangle_provenance[g].size();
                for (const auto& f : triangle_provenance[g])
                    pf << " " << f[0] << " " << vmap[f[1]] << " " << vmap[f[2]] << " " << vmap[f[3]];
                pf << "\n";
            }
        }

        // Sidecar 2: the coplanar group of each input triangle, so the group areas can
        // be reconstructed. "<num_constraints> <num_groups>" then, per real constraint,
        // "<input_triangle_id> <coplanar_group_id>". input_triangle_id is the triangle's
        // index in the INPUT file (degenerate/collinear input triangles are dropped by
        // the mesher, so this is not simply 0..n-1 when the input has any).
        {
            ofstream gf((std::string(filename) + ".groups").c_str(), std::ios::binary);
            gf << "# input_triangle_id coplanar_group_id\n";
            gf << constraint_coplanar_group.size() << " " << coplanar_group_size.size() << "\n";
            for (uint32_t c = 0; c < constraint_coplanar_group.size(); c++) {
                const uint32_t off = constraint_original_index.empty() ? c : constraint_original_index[c];
                gf << off << " " << constraint_coplanar_group[c] << "\n";
            }
        }

        // Sidecar 3: provenance of inserted edges. One line per input edge:
        // "<edge_id> <num_out_edges> <tet v0 v1> <tet v0 v1> ...", where tet is the output
        // tet index and v0 v1 are its edge's output vertex indices (matching the .tet block).
        // The listed output tet edges tile the input segment.
        if (!edge_provenance.empty()) {
            ofstream ef((std::string(filename) + ".edgeprov").c_str(), std::ios::binary);
            ef << "# edge_id num_out_edges tet v0 v1 ...\n";
            ef << edge_provenance.size() << "\n";
            for (const EdgeProvenance& ep : edge_provenance) {
                ef << ep.edge_id << " " << ep.out_edges.size();
                for (const auto& oe : ep.out_edges)
                    ef << " " << oe[0] << " " << vmap[oe[1]] << " " << vmap[oe[2]];
                ef << "\n";
            }
        }

        // Sidecar 4: provenance of inserted points. One line per input point:
        // "<point_id> <tet> <out_vertex>", where out_vertex is the output vertex equal to the
        // point and tet is an output tet containing it (both -1 if the point did not survive).
        if (!point_provenance.empty()) {
            ofstream pf((std::string(filename) + ".pointprov").c_str(), std::ios::binary);
            pf << "# point_id tet out_vertex(-1 -1 if absent)\n";
            pf << point_provenance.size() << "\n";
            for (const PointProvenance& pp : point_provenance) {
                if (pp.out_vertex == UINT32_MAX)
                    pf << pp.point_id << " -1 -1\n";
                else
                    pf << pp.point_id << " " << pp.tet << " " << vmap[pp.out_vertex] << "\n";
            }
        }

        // Sidecar 2 (optional, for exact verification): the output vertex coordinates
        // as exact rationals, same order/indexing as the .tet vertex block. Format per
        // line: "<x> <y> <z>", each a rational "[-]num[/den]" (portable decimal).
        if (export_rational) {
            ofstream rf((std::string(filename) + ".rational").c_str(), std::ios::binary);
            rf << final_numver << "\n";
            bigrational rx, ry, rz;
            for (uint32_t v = 0; v < vertices.size(); v++)
                if (vrts_visit[v]) {
                    vertices[v]->getExactXYZCoordinates(rx, ry, rz);
                    rf << rational_to_string(rx) << " " << rational_to_string(ry) << " "
                       << rational_to_string(rz) << "\n";
                }
        }
    } else {
        size_t internal_cell_num = 0;
        std::vector<uint32_t> face_visit(faces.size(), 0);
        for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
        for (size_t i = 0; i < edge_visit.size(); i++) edge_visit[i] = 0;
        for (BSPcell& cell : cells)
            if (cell.place == INTERNAL_A) {
                internal_cell_num++;
                for (uint64_t f : cell.faces) face_visit[f] = 1;
            }
        for (size_t f = 0; f < faces.size(); f++)
            if (face_visit[f])
                for (uint64_t e : faces[f].edges) edge_visit[e] = 1;
        for (size_t e = 0; e < edges.size(); e++)
            if (edge_visit[e])
                vrts_visit[edges[e].vertices[0]] = vrts_visit[edges[e].vertices[1]] = 1;
        uint32_t final_numver = 0;
        for (size_t i = 0; i < vrts_visit.size(); i++)
            if (vrts_visit[i]) final_numver++;

        std::vector<uint32_t> vmap(vertices.size(), 0);
        size_t num_v = 0;
        for (size_t i = 0; i < vrts_visit.size(); i++) {
            vmap[i] = (uint32_t)num_v;
            if (vrts_visit[i]) num_v++;
        }

        f << final_numver << " vertices\n";
        f << internal_cell_num << " cells\n";
        // Print vertices coordinates
        for (uint32_t v = 0; v < vertices.size(); v++)
            if (vrts_visit[v]) {
                f << (*vertices[v]) << "\n";
            }

        // Print cells
        for (size_t v = 0; v < vertices.size(); v++) vrts_visit[v] = 0;
        for (size_t e = 0; e < edges.size(); e++) edge_visit[e] = 0;
        for (BSPcell& cell : cells)
            if (cell.place == INTERNAL_A) {
                uint64_t numce = count_cellEdges(cell);
                uint32_t numcv = count_cellVertices(cell, &numce);
                std::vector<uint32_t> cell_vrts(numcv);
                list_cellVertices(cell, numce, cell_vrts);
                f << cell_vrts.size() << " ";
                for (uint32_t i : cell_vrts) f << vmap[i] << " ";
                f << "\n";
            }
    }

    f.close();

#ifndef NDEBUG
    // In tetrahedralization mode, every emitted tetrahedron must have strictly
    // positive volume. makeTetrahedra() guarantees this (it drops zero-volume
    // tets and orients the rest positively); assert it with the exact orient3D
    // predicate -- no floating point -- so any regression is caught immediately.
    if (tetrahedrize)
        for (uint32_t t = 0; t < final_tets.size(); t += 4)
            assert(genericPoint::orient3D(
                       *vertices[final_tets[t]],
                       *vertices[final_tets[t + 1]],
                       *vertices[final_tets[t + 2]],
                       *vertices[final_tets[t + 3]]) > 0 &&
                   "saveMesh: generated tetrahedron does not have positive volume");
#endif

    final_tets.clear();
}


void BSPcomplex::saveBlackFaces(const char* filename, bool triangulate)
{
    // Binary mode so line endings stay LF on every OS (byte-identical output).
    ofstream f(filename, std::ios::binary);

    if (!f) ip_error("BSPcomplex::saveBlackFaces: cannot open the file.\n");

    if (triangulate) {
        const uint64_t num_faces = faces.size();
        for (uint64_t face_ind = 0; face_ind < num_faces; face_ind++) triangulateFace(face_ind);
    }

    for (size_t i = 0; i < vrts_visit.size(); i++) vrts_visit[i] = 0;
    for (size_t i = 0; i < edges.size(); i++) edge_visit[i] = 0;

    uint64_t num_border_faces = 0;
    for (uint64_t f_i = 0; f_i < faces.size(); f_i++)
        if (faces[f_i].colour != WHITE) {
            num_border_faces++;
            for (uint64_t eid : faces[f_i].edges) edge_visit[eid] = 1;
        }
    for (size_t i = 0; i < edges.size(); i++)
        if (edge_visit[i]) vrts_visit[edges[i].vertices[0]] = vrts_visit[edges[i].vertices[1]] = 1;

    std::vector<uint32_t> vmap(vertices.size(), 0);
    size_t num_v = 0;
    for (size_t i = 0; i < vrts_visit.size(); i++) {
        vmap[i] = (uint32_t)num_v;
        if (vrts_visit[i]) num_v++;
    }

    f << "OFF\n";
    f << num_v << " ";
    f << num_border_faces << " ";
    f << "0\n";

    // Print vertices coordinates
    for (uint32_t v = 0; v < vertices.size(); v++)
        if (vrts_visit[v]) f << (*(vertices)[v]) << "\n";

    // Print border faces
    for (uint32_t f_i = 0; f_i < faces.size(); f_i++)
        if (faces[f_i].colour != WHITE) {
            BSPface& face = faces[f_i];
            vector<uint32_t> face_vrts(face.edges.size(), UINT32_MAX);
            list_faceVertices(face, face_vrts);
            f << face_vrts.size();

            for (uint32_t v = 0; v < face_vrts.size(); v++) f << " " << vmap[face_vrts[v]];

            f << "\n";
        }

    f.close();
}


size_t BSPcomplex::getStructureSize() const
{
    size_t tot = 0;

    // Size for points
    for (genericPoint* p : vertices) {
        if (p->isExplicit3D())
            tot += sizeof(explicitPoint3D);
        else if (p->isLPI())
            tot += sizeof(implicitPoint3D_LPI);
        else
            tot += sizeof(implicitPoint3D_TPI);
    }

    // Size of the array of pointers
    tot += vertices.size() * sizeof(genericPoint*);

    // Size of edge objects
    tot += edges.size() * sizeof(BSPedge);

    // Size of face objects (including pointed arrays)
    for (const BSPface& f : faces) tot += f.getSize();

    // Size of cell objects (including pointed arrays)
    for (const BSPcell& c : cells) tot += c.getSize();

    // And all the other vectors use by the structure...
    tot += sizeof(uint32_t) * constraints_vrts.size();
    tot += sizeof(CONSTR_GROUP_T) * constraint_group.size();
    tot += sizeof(uint32_t) * final_tets.size();
    tot += sizeof(char) * vrts_orBin.size();
    tot += sizeof(uint32_t) * vrts_visit.size();
    tot += sizeof(uint64_t) * edge_visit.size();
    tot += sizeof(BSPcomplex);

    return tot;
}
} // namespace vol_rem