// volume_remesher._core: nanobind module over driver.h.
//
// This layer only moves data: numpy arrays in, the driver's snapshots out as read-only numpy
// views that keep their owner alive, and the writers of mesh_output.h / io2d.h for save(). The
// Python package (python/volume_remesher/) builds the user-facing API on top of it.
//
// Two process-level concerns live here:
//  - The pipeline runs on a dedicated thread with a large stack. mesh_generator links with an
//    8 MB stack; python.exe's main thread on Windows has 1 MB, and the arrangement recurses deep
//    enough to need more. The GIL is released meanwhile.
//  - ip_error throws (NFG_IP_ERROR_THROWS) and the exception reaches Python as
//    VolumeRemesherError, see the translator at the bottom.

// Platform headers first. On MSVC, windows.h must precede any header that brings `using
// namespace std` into scope (the library's do), or `byte` becomes ambiguous between std::byte
// and the Windows SDK's; lean-and-mean also keeps the COM headers that trip on it out entirely.
#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

#include <VolumeRemesher/driver.h>
#include <VolumeRemesher/embed.h>
#include <VolumeRemesher/2d/embed2d.h>
#include <VolumeRemesher/parallel.h>

#include <cfenv>
#include <cstdint>
#include <exception>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;
using namespace vol_rem;

namespace {

// ---------------------------------------------------------------------------------------------
// Run a callable on a thread with a large stack, propagating its exception.
// ---------------------------------------------------------------------------------------------
struct BigStackTask
{
    std::function<void()> fn;
    std::exception_ptr error;
    void run()
    {
        try {
            fn();
        } catch (...) {
            error = std::current_exception();
        }
    }
};

#ifdef _WIN32
unsigned __stdcall big_stack_entry(void* p)
{
    static_cast<BigStackTask*>(p)->run();
    return 0;
}
#else
void* big_stack_entry(void* p)
{
    static_cast<BigStackTask*>(p)->run();
    return nullptr;
}
#endif

void run_on_big_stack(std::function<void()> fn)
{
    // Reserved, not committed: only the pages the pipeline touches cost memory.
    constexpr size_t kStackBytes = size_t(64) << 20;
    BigStackTask task{std::move(fn), nullptr};
#ifdef _WIN32
    HANDLE h = (HANDLE)_beginthreadex(
        nullptr, (unsigned)kStackBytes, big_stack_entry, &task, STACK_SIZE_PARAM_IS_A_RESERVATION,
        nullptr);
    if (!h) throw std::runtime_error("volume_remesher: cannot start the worker thread\n");
    WaitForSingleObject(h, INFINITE);
    CloseHandle(h);
#else
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, kStackBytes);
    pthread_t th;
    const int rc = pthread_create(&th, &attr, big_stack_entry, &task);
    pthread_attr_destroy(&attr);
    if (rc != 0) throw std::runtime_error("volume_remesher: cannot start the worker thread\n");
    pthread_join(th, nullptr);
#endif
    if (task.error) std::rethrow_exception(task.error);
}

// ---------------------------------------------------------------------------------------------
// numpy <-> std::vector
// ---------------------------------------------------------------------------------------------
using F64x3 = nb::ndarray<const double, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>;
using F64x2 = nb::ndarray<const double, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>;
using U32x2 = nb::ndarray<const uint32_t, nb::shape<-1, 2>, nb::c_contig, nb::device::cpu>;
using U32x3 = nb::ndarray<const uint32_t, nb::shape<-1, 3>, nb::c_contig, nb::device::cpu>;
using U32x4 = nb::ndarray<const uint32_t, nb::shape<-1, 4>, nb::c_contig, nb::device::cpu>;

template <class T, class Arr>
std::vector<T> to_vector(const Arr& a)
{
    return std::vector<T>(a.data(), a.data() + a.size());
}

TriMeshData to_trimesh(const std::pair<F64x3, U32x3>& m)
{
    TriMeshData t;
    t.coords = to_vector<double>(m.first);
    t.tris = to_vector<uint32_t>(m.second);
    return t;
}

EdgeData to_edges(const std::pair<F64x3, U32x2>& e)
{
    EdgeData d;
    d.coords = to_vector<double>(e.first);
    d.edges = to_vector<uint32_t>(e.second);
    return d;
}

PointData to_points(const F64x3& p)
{
    PointData d;
    d.coords = to_vector<double>(p);
    return d;
}

TetMeshData to_tetmesh(const std::pair<F64x3, U32x4>& m)
{
    TetMeshData t;
    t.coords = to_vector<double>(m.first);
    t.tets = to_vector<uint32_t>(m.second);
    return t;
}

// A read-only 2-D numpy view (rows x cols) of a vector owned by `owner`. An empty vector gets a
// valid pointer to a static buffer so numpy never sees a null data pointer.
template <class T>
nb::ndarray<nb::numpy, const T> view2(nb::handle owner, const std::vector<T>& v, size_t cols)
{
    static const T empty[4] = {};
    const T* data = v.empty() ? empty : v.data();
    const size_t shape[2] = {cols ? v.size() / cols : 0, cols};
    return nb::ndarray<nb::numpy, const T>(data, 2, shape, owner);
}

template <class T, size_t N>
nb::ndarray<nb::numpy, const T> view2(nb::handle owner, const std::vector<std::array<T, N>>& v)
{
    static const T empty[N] = {};
    const T* data = v.empty() ? empty : v.data()->data();
    const size_t shape[2] = {v.size(), N};
    return nb::ndarray<nb::numpy, const T>(data, 2, shape, owner);
}

template <class T>
nb::ndarray<nb::numpy, const T> view1(nb::handle owner, const std::vector<T>& v)
{
    static const T empty[1] = {};
    const T* data = v.empty() ? empty : v.data();
    const size_t shape[1] = {v.size()};
    return nb::ndarray<nb::numpy, const T>(data, 1, shape, owner);
}

// A numpy array that owns a moved-in vector (for the readers).
template <class T>
nb::ndarray<nb::numpy, T> own2(std::vector<T>&& v, size_t cols)
{
    auto* holder = new std::vector<T>(std::move(v));
    if (holder->empty()) holder->reserve(1); // a valid data pointer for the 0-row array
    nb::capsule owner(holder, [](void* p) noexcept { delete static_cast<std::vector<T>*>(p); });
    const size_t shape[2] = {cols ? holder->size() / cols : 0, cols};
    return nb::ndarray<nb::numpy, T>(holder->data(), 2, shape, owner);
}

// ---------------------------------------------------------------------------------------------
// Result holders. The driver's snapshots plus flattened (offsets + rows) copies of the nested
// provenance lists, so Python gets a few arrays instead of millions of small objects.
// ---------------------------------------------------------------------------------------------
struct Mesh
{
    MeshOutput out;
    std::vector<uint32_t> tri_prov_offsets; // per coplanar group, into tri_prov
    std::vector<std::array<uint32_t, 4>> tri_prov; // {tet, v0, v1, v2}
    std::vector<uint32_t> edge_ids;
    std::vector<uint32_t> edge_prov_offsets; // per input edge, into edge_prov
    std::vector<std::array<uint32_t, 3>> edge_prov; // {tet, v0, v1}
    std::vector<std::array<uint32_t, 3>> point_prov; // {point_id, tet, vertex}

    void flatten()
    {
        const VolumeOutput& v = out.volume;
        tri_prov_offsets.reserve(v.triangle_provenance.size() + 1);
        tri_prov_offsets.push_back(0);
        for (const auto& g : v.triangle_provenance) {
            tri_prov.insert(tri_prov.end(), g.begin(), g.end());
            tri_prov_offsets.push_back((uint32_t)tri_prov.size());
        }
        edge_prov_offsets.reserve(v.edge_provenance.size() + 1);
        edge_prov_offsets.push_back(0);
        for (const auto& e : v.edge_provenance) {
            edge_ids.push_back(e.edge_id);
            edge_prov.insert(edge_prov.end(), e.out_edges.begin(), e.out_edges.end());
            edge_prov_offsets.push_back((uint32_t)edge_prov.size());
        }
        point_prov.reserve(v.point_provenance.size());
        for (const auto& p : v.point_provenance)
            point_prov.push_back({p.point_id, p.tet, p.out_vertex});
    }
};

struct Arrangement
{
    vr2d::ArrangementOutput out;
    std::vector<uint32_t> seg_prov_offsets; // per input segment, into seg_prov
    std::vector<std::array<uint32_t, 3>> seg_prov; // {tri, v0, v1}

    void flatten()
    {
        seg_prov_offsets.reserve(out.segment_provenance.size() + 1);
        seg_prov_offsets.push_back(0);
        for (const auto& s : out.segment_provenance) {
            seg_prov.insert(seg_prov.end(), s.begin(), s.end());
            seg_prov_offsets.push_back((uint32_t)seg_prov.size());
        }
    }
};

MeshOptions make_options(
    bool verbose,
    bool logging,
    bool triangulate,
    bool export_rational,
    bool fill_holes,
    bool keep_all_cells,
    MeshOutput::Kind output)
{
    MeshOptions o;
    o.verbose = verbose;
    o.logging = logging;
    o.triangulate = triangulate;
    o.export_rational = export_rational;
    o.fill_holes = fill_holes;
    o.keep_all_cells = keep_all_cells;
    o.output = output;
    return o;
}

char parse_opcode(const std::string& op)
{
    if (op.size() != 1 || (op[0] != '0' && op[0] != 'U' && op[0] != 'I' && op[0] != 'D'))
        throw std::invalid_argument("opcode must be one of '0', 'U', 'I', 'D'");
    return op[0];
}

// ---------------------------------------------------------------------------------------------
// Pipeline entry points
// ---------------------------------------------------------------------------------------------
std::unique_ptr<Mesh> py_make_mesh(
    std::optional<std::pair<F64x3, U32x3>> a,
    std::optional<std::string> name_a,
    const std::string& opcode,
    std::optional<std::pair<F64x3, U32x3>> b,
    std::optional<std::string> name_b,
    std::optional<std::pair<F64x3, U32x2>> edges,
    std::optional<F64x3> points,
    bool verbose,
    bool logging,
    bool triangulate,
    bool export_rational,
    bool fill_holes,
    bool keep_all_cells,
    MeshOutput::Kind output)
{
    const char op = parse_opcode(opcode);
    const MeshOptions opts = make_options(
        verbose, logging, triangulate, export_rational, fill_holes, keep_all_cells, output);
    // Copy the inputs while the GIL is held; the arrays may be temporaries.
    std::optional<TriMeshData> A, B;
    if (a) A = to_trimesh(*a);
    if (b) B = to_trimesh(*b);
    std::optional<EdgeData> E;
    if (edges) E = to_edges(*edges);
    std::optional<PointData> P;
    if (points) P = to_points(*points);

    auto result = std::make_unique<Mesh>();
    {
        nb::gil_scoped_release release;
        run_on_big_stack([&]() {
            result->out = make_mesh(
                A ? &*A : nullptr,
                name_a ? name_a->c_str() : nullptr,
                op,
                B ? &*B : nullptr,
                name_b ? name_b->c_str() : nullptr,
                E ? &*E : nullptr,
                P ? &*P : nullptr,
                opts);
            result->flatten();
        });
    }
    return result;
}

std::unique_ptr<Mesh> py_embed_mesh(
    std::pair<F64x3, U32x3> a,
    std::pair<F64x3, U32x4> tet_mesh,
    const std::string& opcode,
    std::optional<std::pair<F64x3, U32x2>> edges,
    std::optional<F64x3> points,
    bool verbose,
    bool triangulate,
    bool export_rational,
    MeshOutput::Kind output)
{
    const char op = parse_opcode(opcode);
    const MeshOptions opts =
        make_options(verbose, false, triangulate, export_rational, false, false, output);
    TriMeshData A = to_trimesh(a);
    TetMeshData T = to_tetmesh(tet_mesh);
    std::optional<EdgeData> E;
    if (edges) E = to_edges(*edges);
    std::optional<PointData> P;
    if (points) P = to_points(*points);

    auto result = std::make_unique<Mesh>();
    {
        nb::gil_scoped_release release;
        run_on_big_stack([&]() {
            result->out = embed_mesh(A, T, op, E ? &*E : nullptr, P ? &*P : nullptr, opts);
            result->flatten();
        });
    }
    return result;
}

std::unique_ptr<Arrangement> py_arrange_2d(
    F64x2 points, U32x2 segments, bool export_rational, bool verbose)
{
    std::vector<double> coords = to_vector<double>(points);
    std::vector<uint32_t> idx = to_vector<uint32_t>(segments);
    auto result = std::make_unique<Arrangement>();
    {
        nb::gil_scoped_release release;
        run_on_big_stack([&]() {
            result->out = arrange_2d(coords, idx, export_rational, verbose);
            result->flatten();
        });
    }
    return result;
}

// Portable decimal "[-]num[/den]" of an exact rational (as BSP.cpp's rational_to_string).
std::string rational_str(const bigrational& r)
{
#ifdef USE_GNU_GMP_CLASSES
    return r.get_str();
#else
    return r.get_dec_str();
#endif
}

// The outputs of embed_tri_in_poly_mesh, as plain data: exact coordinates as strings, every
// nested list as offsets + rows. Built on the worker thread, where the bigrationals live.
struct PolyhedralEmbedding
{
    std::vector<std::string> exact_coords; // 3 per vertex
    std::vector<uint32_t> facet_offsets, facet_vertices; // polygons
    std::vector<uint32_t> cell_offsets, cell_facets; // polyhedra as facet lists
    std::vector<std::array<uint32_t, 4>> tets;
    std::vector<uint32_t> tets_parent;
    std::vector<uint32_t> facets_on_input;
    std::vector<uint8_t> cells_with_faces_on_input;
    std::vector<uint32_t> tets_parent_faces_offsets, tets_parent_faces;
    std::vector<uint32_t> tri_prov_offsets;
    std::vector<std::array<uint32_t, 4>> tri_prov;
    std::vector<uint32_t> tri_group;
    std::vector<uint32_t> edge_prov_offsets;
    std::vector<std::array<uint32_t, 3>> edge_prov;
    std::vector<std::array<uint32_t, 2>> point_prov;
};

// "(n, i1..in)" sequences -> offsets + flat items.
void unpack_counted(
    const std::vector<uint32_t>& seq, std::vector<uint32_t>& offsets, std::vector<uint32_t>& flat)
{
    offsets.push_back(0);
    for (size_t k = 0; k < seq.size();) {
        const uint32_t n = seq[k++];
        for (uint32_t i = 0; i < n && k < seq.size(); i++) flat.push_back(seq[k++]);
        offsets.push_back((uint32_t)flat.size());
    }
}

template <class Row>
void flatten_nested(
    const std::vector<std::vector<Row>>& nested, std::vector<uint32_t>& offsets,
    std::vector<Row>& flat)
{
    offsets.reserve(nested.size() + 1);
    offsets.push_back(0);
    for (const auto& v : nested) {
        flat.insert(flat.end(), v.begin(), v.end());
        offsets.push_back((uint32_t)flat.size());
    }
}

// Runs inside run_on_big_stack: the bigrationals are born and die on that thread.
void run_embed_polyhedral(
    const TriMeshData& A,
    const TetMeshData& T,
    const EdgeData* E,
    const PointData* P,
    bool verbose,
    PolyhedralEmbedding& r)
{
    std::vector<bigrational> vertices;
    std::vector<uint32_t> facets, cells;
    std::vector<std::array<uint32_t, 4>> out_tets;
    std::vector<uint32_t> final_tets_parent, facets_on_input;
    std::vector<bool> cells_with_faces_on_input;
    std::vector<std::vector<uint32_t>> final_tets_parent_faces;
    std::vector<std::vector<std::array<uint32_t, 4>>> tri_prov;
    std::vector<uint32_t> tri_group;
    std::vector<std::vector<std::array<uint32_t, 3>>> edge_prov;
    std::vector<std::array<uint32_t, 2>> point_prov;
    static const std::vector<double> no_coords;
    static const std::vector<uint32_t> no_idx;

    embed_tri_in_poly_mesh(
        A.coords, A.tris, T.coords, T.tets, vertices, facets, cells, out_tets, final_tets_parent,
        facets_on_input, cells_with_faces_on_input, final_tets_parent_faces,
        E ? E->coords : no_coords, E ? E->edges : no_idx, P ? P->coords : no_coords, tri_prov,
        tri_group, edge_prov, point_prov, verbose);

    r.exact_coords.reserve(vertices.size());
    for (const bigrational& c : vertices) r.exact_coords.push_back(rational_str(c));
    unpack_counted(facets, r.facet_offsets, r.facet_vertices);
    unpack_counted(cells, r.cell_offsets, r.cell_facets);
    r.tets = std::move(out_tets);
    r.tets_parent = std::move(final_tets_parent);
    r.facets_on_input = std::move(facets_on_input);
    r.cells_with_faces_on_input.assign(
        cells_with_faces_on_input.begin(), cells_with_faces_on_input.end());
    flatten_nested(final_tets_parent_faces, r.tets_parent_faces_offsets, r.tets_parent_faces);
    flatten_nested(tri_prov, r.tri_prov_offsets, r.tri_prov);
    r.tri_group = std::move(tri_group);
    flatten_nested(edge_prov, r.edge_prov_offsets, r.edge_prov);
    r.point_prov = std::move(point_prov);
}

std::unique_ptr<PolyhedralEmbedding> py_embed_polyhedral(
    std::pair<F64x3, U32x3> a,
    std::pair<F64x3, U32x4> tet_mesh,
    std::optional<std::pair<F64x3, U32x2>> edges,
    std::optional<F64x3> points,
    bool verbose)
{
    TriMeshData A = to_trimesh(a);
    TetMeshData T = to_tetmesh(tet_mesh);
    std::optional<EdgeData> E;
    if (edges) E = to_edges(*edges);
    std::optional<PointData> P;
    if (points) P = to_points(*points);
    auto result = std::make_unique<PolyhedralEmbedding>();
    {
        nb::gil_scoped_release release;
        run_on_big_stack([&]() {
            run_embed_polyhedral(A, T, E ? &*E : nullptr, P ? &*P : nullptr, verbose, *result);
        });
    }
    return result;
}

// The outputs of embed_seg_in_tri_mesh (2D), likewise.
struct SegmentEmbedding
{
    std::vector<std::string> exact_coords; // 2 per vertex
    std::vector<std::array<uint32_t, 3>> tris;
    std::vector<uint32_t> seg_prov_offsets;
    std::vector<std::array<uint32_t, 3>> seg_prov;
    std::vector<std::array<uint32_t, 2>> point_prov;
};

std::unique_ptr<SegmentEmbedding> py_embed_segments(F64x2 points, U32x2 segments, bool verbose)
{
    std::vector<double> coords = to_vector<double>(points);
    std::vector<uint32_t> idx = to_vector<uint32_t>(segments);
    auto result = std::make_unique<SegmentEmbedding>();
    {
        nb::gil_scoped_release release;
        run_on_big_stack([&]() {
            std::vector<bigrational> vertices;
            std::vector<std::vector<std::array<uint32_t, 3>>> seg_prov;
            if (!embed_seg_in_tri_mesh(
                    coords, idx, vertices, result->tris, seg_prov, result->point_prov, verbose))
                throw std::runtime_error("Failed to build the 2D arrangement\n");
            result->exact_coords.reserve(vertices.size());
            for (const bigrational& c : vertices) result->exact_coords.push_back(rational_str(c));
            flatten_nested(seg_prov, result->seg_prov_offsets, result->seg_prov);
        });
    }
    return result;
}

} // namespace

NB_MODULE(_core, m)
{
    m.doc() = "Low-level bindings of VolumeRemesher's driver. Use the volume_remesher package.";

    // std::runtime_error is what ip_error throws; std::invalid_argument maps to ValueError
    // through nanobind's default translation.
    nb::exception<std::runtime_error>(m, "VolumeRemesherError");

    nb::enum_<MeshOutput::Kind>(m, "OutputKind")
        .value("VOLUME", MeshOutput::VOLUME)
        .value("SKIN", MeshOutput::SKIN)
        .value("BLACK_FACES", MeshOutput::BLACK_FACES);

    // --- readers (same parsing as the command line tool) ---------------------------------
    m.def(
        "read_off",
        [](const std::string& path, bool verbose) {
            TriMeshData t = read_OFF(path.c_str(), verbose);
            return std::make_pair(own2(std::move(t.coords), 3), own2(std::move(t.tris), 3));
        },
        "path"_a, "verbose"_a = false);
    m.def(
        "read_edge_off",
        [](const std::string& path, bool verbose) {
            EdgeData e = read_edge_OFF(path.c_str(), verbose);
            return std::make_pair(own2(std::move(e.coords), 3), own2(std::move(e.edges), 2));
        },
        "path"_a, "verbose"_a = false);
    m.def(
        "read_points_off",
        [](const std::string& path, bool verbose) {
            PointData p = read_points_OFF(path.c_str(), verbose);
            return own2(std::move(p.coords), 3);
        },
        "path"_a, "verbose"_a = false);
    m.def(
        "read_tet",
        [](const std::string& path, bool verbose) {
            TetMeshData t = read_TET(path.c_str(), verbose);
            return std::make_pair(own2(std::move(t.coords), 3), own2(std::move(t.tets), 4));
        },
        "path"_a, "verbose"_a = false);
    m.def(
        "read_medit",
        [](const std::string& path, bool verbose) {
            TetMeshData t = read_MEDIT(path.c_str(), verbose);
            return std::make_pair(own2(std::move(t.coords), 3), own2(std::move(t.tets), 4));
        },
        "path"_a, "verbose"_a = false);
    m.def(
        "read_obj_segments",
        [](const std::string& path) {
            std::vector<double> coords;
            std::vector<uint32_t> idx;
            if (!vr2d::read_OBJ_segments(path, coords, idx))
                throw std::runtime_error("Cannot read the input OBJ line soup\n");
            return std::make_pair(own2(std::move(coords), 2), own2(std::move(idx), 2));
        },
        "path"_a);

    // --- results -----------------------------------------------------------------------------
    nb::class_<Mesh>(m, "Mesh")
        .def_prop_ro("kind", [](const Mesh& s) { return s.out.kind; })
        // volume
        .def_prop_ro("tetrahedral", [](const Mesh& s) { return s.out.volume.tetrahedral; })
        .def_prop_ro(
            "coords",
            [](nb::handle_t<Mesh> h) {
                const Mesh& s = nb::cast<const Mesh&>(h);
                return view2(
                    h, s.out.kind == MeshOutput::VOLUME ? s.out.volume.coords : s.out.surface.coords,
                    3);
            })
        .def_prop_ro(
            "tets",
            [](nb::handle_t<Mesh> h) { return view2(h, nb::cast<const Mesh&>(h).out.volume.tets, 4); })
        .def_prop_ro(
            "cell_offsets",
            [](nb::handle_t<Mesh> h) {
                return view1(h, nb::cast<const Mesh&>(h).out.volume.cell_offsets);
            })
        .def_prop_ro(
            "cell_vertices",
            [](nb::handle_t<Mesh> h) {
                return view1(h, nb::cast<const Mesh&>(h).out.volume.cell_vertices);
            })
        .def_prop_ro(
            "exact_coords", [](const Mesh& s) { return s.out.volume.exact_coords; },
            "3 strings per vertex, each '[-]num[/den]'; empty unless exact coordinates were asked for")
        .def_prop_ro(
            "triangle_provenance_offsets",
            [](nb::handle_t<Mesh> h) { return view1(h, nb::cast<const Mesh&>(h).tri_prov_offsets); })
        .def_prop_ro(
            "triangle_provenance",
            [](nb::handle_t<Mesh> h) { return view2(h, nb::cast<const Mesh&>(h).tri_prov); })
        .def_prop_ro(
            "constraint_groups",
            [](nb::handle_t<Mesh> h) {
                return view2(h, nb::cast<const Mesh&>(h).out.volume.constraint_groups);
            })
        .def_prop_ro(
            "num_coplanar_groups", [](const Mesh& s) { return s.out.volume.num_coplanar_groups; })
        .def_prop_ro(
            "edge_ids", [](nb::handle_t<Mesh> h) { return view1(h, nb::cast<const Mesh&>(h).edge_ids); })
        .def_prop_ro(
            "edge_provenance_offsets",
            [](nb::handle_t<Mesh> h) { return view1(h, nb::cast<const Mesh&>(h).edge_prov_offsets); })
        .def_prop_ro(
            "edge_provenance",
            [](nb::handle_t<Mesh> h) { return view2(h, nb::cast<const Mesh&>(h).edge_prov); })
        .def_prop_ro(
            "point_provenance",
            [](nb::handle_t<Mesh> h) { return view2(h, nb::cast<const Mesh&>(h).point_prov); })
        // surface (skin / black faces)
        .def_prop_ro(
            "face_offsets",
            [](nb::handle_t<Mesh> h) {
                return view1(h, nb::cast<const Mesh&>(h).out.surface.face_offsets);
            })
        .def_prop_ro(
            "face_vertices",
            [](nb::handle_t<Mesh> h) {
                return view1(h, nb::cast<const Mesh&>(h).out.surface.face_vertices);
            })
        .def(
            "save",
            [](const Mesh& s, const std::string& path) {
                nb::gil_scoped_release release;
                write_mesh_output(path.c_str(), s.out);
            },
            "path"_a,
            "Write the file(s) mesh_generator would write for this result, byte for byte: the "
            "mesh at `path`, plus the .triangleprov/.groups/.edgeprov/.pointprov/.rational "
            "sidecars for a tetrahedral mesh.");

    nb::class_<Arrangement>(m, "Arrangement")
        .def_prop_ro(
            "coords",
            [](nb::handle_t<Arrangement> h) {
                return view2(h, nb::cast<const Arrangement&>(h).out.coords, 2);
            })
        .def_prop_ro(
            "triangles",
            [](nb::handle_t<Arrangement> h) {
                return view2(h, nb::cast<const Arrangement&>(h).out.tris);
            })
        .def_prop_ro("exact_coords", [](const Arrangement& s) { return s.out.exact_coords; })
        .def_prop_ro(
            "segment_provenance_offsets",
            [](nb::handle_t<Arrangement> h) {
                return view1(h, nb::cast<const Arrangement&>(h).seg_prov_offsets);
            })
        .def_prop_ro(
            "segment_provenance",
            [](nb::handle_t<Arrangement> h) {
                return view2(h, nb::cast<const Arrangement&>(h).seg_prov);
            })
        .def_prop_ro(
            "point_provenance",
            [](nb::handle_t<Arrangement> h) {
                return view2(h, nb::cast<const Arrangement&>(h).out.point_provenance);
            })
        .def(
            "save",
            [](const Arrangement& s, const std::string& base) {
                nb::gil_scoped_release release;
                if (!vr2d::write_arrangement(base, s.out))
                    throw std::runtime_error("Cannot write the output triangulation\n");
            },
            "path"_a,
            "Write the OFF triangulation at `path` plus its .segmentprov/.pointprov (and "
            ".rational when exact coordinates were asked for), byte for byte as mesh_generator -2d.");

    // --- pipeline ----------------------------------------------------------------------------
    m.def(
        "make_mesh", &py_make_mesh, "a"_a.none(), "name_a"_a.none(), "opcode"_a, "b"_a.none(),
        "name_b"_a.none(), "edges"_a.none(), "points"_a.none(), "verbose"_a, "logging"_a,
        "triangulate"_a, "export_rational"_a, "fill_holes"_a, "keep_all_cells"_a, "output"_a);
    m.def(
        "embed_mesh", &py_embed_mesh, "a"_a, "tet_mesh"_a, "opcode"_a, "edges"_a.none(),
        "points"_a.none(), "verbose"_a, "triangulate"_a, "export_rational"_a, "output"_a);
    m.def(
        "arrange_2d", &py_arrange_2d, "points"_a, "segments"_a, "export_rational"_a, "verbose"_a);

    nb::class_<PolyhedralEmbedding>(m, "PolyhedralEmbedding")
        .def_prop_ro("exact_coords", [](const PolyhedralEmbedding& s) { return s.exact_coords; })
        .def_prop_ro("facet_offsets", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).facet_offsets); })
        .def_prop_ro("facet_vertices", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).facet_vertices); })
        .def_prop_ro("cell_offsets", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).cell_offsets); })
        .def_prop_ro("cell_facets", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).cell_facets); })
        .def_prop_ro("tets", [](nb::handle_t<PolyhedralEmbedding> h) { return view2(h, nb::cast<const PolyhedralEmbedding&>(h).tets); })
        .def_prop_ro("tets_parent", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).tets_parent); })
        .def_prop_ro("facets_on_input", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).facets_on_input); })
        .def_prop_ro("cells_with_faces_on_input", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).cells_with_faces_on_input); })
        .def_prop_ro("tets_parent_faces_offsets", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).tets_parent_faces_offsets); })
        .def_prop_ro("tets_parent_faces", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).tets_parent_faces); })
        .def_prop_ro("triangle_provenance_offsets", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).tri_prov_offsets); })
        .def_prop_ro("triangle_provenance", [](nb::handle_t<PolyhedralEmbedding> h) { return view2(h, nb::cast<const PolyhedralEmbedding&>(h).tri_prov); })
        .def_prop_ro("triangle_group", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).tri_group); })
        .def_prop_ro("edge_provenance_offsets", [](nb::handle_t<PolyhedralEmbedding> h) { return view1(h, nb::cast<const PolyhedralEmbedding&>(h).edge_prov_offsets); })
        .def_prop_ro("edge_provenance", [](nb::handle_t<PolyhedralEmbedding> h) { return view2(h, nb::cast<const PolyhedralEmbedding&>(h).edge_prov); })
        .def_prop_ro("point_provenance", [](nb::handle_t<PolyhedralEmbedding> h) { return view2(h, nb::cast<const PolyhedralEmbedding&>(h).point_prov); });

    nb::class_<SegmentEmbedding>(m, "SegmentEmbedding")
        .def_prop_ro("exact_coords", [](const SegmentEmbedding& s) { return s.exact_coords; })
        .def_prop_ro("triangles", [](nb::handle_t<SegmentEmbedding> h) { return view2(h, nb::cast<const SegmentEmbedding&>(h).tris); })
        .def_prop_ro("segment_provenance_offsets", [](nb::handle_t<SegmentEmbedding> h) { return view1(h, nb::cast<const SegmentEmbedding&>(h).seg_prov_offsets); })
        .def_prop_ro("segment_provenance", [](nb::handle_t<SegmentEmbedding> h) { return view2(h, nb::cast<const SegmentEmbedding&>(h).seg_prov); })
        .def_prop_ro("point_provenance", [](nb::handle_t<SegmentEmbedding> h) { return view2(h, nb::cast<const SegmentEmbedding&>(h).point_prov); });

    m.def(
        "embed_polyhedral", &py_embed_polyhedral, "a"_a, "tet_mesh"_a, "edges"_a.none(),
        "points"_a.none(), "verbose"_a,
        "embed_tri_in_poly_mesh: cut the tet mesh by the triangles (and edges/points); exact "
        "rational output coordinates in lowest terms, polyhedral cells, tets and provenance.");
    m.def(
        "embed_segments", &py_embed_segments, "points"_a, "segments"_a, "verbose"_a,
        "embed_seg_in_tri_mesh: the 2D arrangement with exact coordinates in lowest terms.");

    // --- threads ---------------------------------------------------------------------------------
    m.def(
        "set_num_threads",
        &vol_rem::set_num_threads,
        "n"_a,
        "Threads of the parallel passes; 0 restores the default (see get_num_threads).");
    m.def(
        "get_num_threads",
        &vol_rem::num_threads,
        "Threads of the parallel passes: set_num_threads(n), else VOLUMEREMESHER_NUM_THREADS, "
        "else the CPUs this process may run on (its affinity mask on Linux).");

    // --- process state, for the tests --------------------------------------------------------
    m.def(
        "rounding_mode_is_nearest", []() { return std::fegetround() == FE_TONEAREST; },
        "True when the calling thread's floating point rounding mode is round-to-nearest. The "
        "pipeline runs on its own thread, so this must stay true for the caller.");
}
