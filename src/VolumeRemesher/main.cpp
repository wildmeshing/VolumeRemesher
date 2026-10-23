#include "driver.h"

#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>

using namespace vol_rem;
using NFG::ip_error;


/// <summary>
/// Main function
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>

static int run(int argc, char** argv);

// With NFG_IP_ERROR_THROWS (the Python build defines it for the whole library) ip_error throws
// instead of exiting; report it as the exit would have, so the tool behaves the same either way.
int main(int argc, char** argv)
{
    try {
        return run(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << e.what();
        return EXIT_FAILURE;
    }
}

static int run(int argc, char** argv)
{
    bool triangulate = false;
    bool verbose = false;
    bool logging = false;
    bool surfmesh = false;
    bool blackfaces = false;
    bool export_rational = false;
    bool fill_holes = false; // -f: opt-in cap of open boundaries (solid-output path)
    bool two_d = false; // -2d: run the 2D segment-arrangement pipeline instead
    bool keep_all_cells = false; // -a: skip in/out min-cut, keep the whole domain
    char* edge_file = NULL; // -e: extra edges to insert
    char* point_file = NULL; // -p: extra points to insert
    char* fileA_name = NULL;
    char* fileB_name = NULL;
    char bool_opcode = '0';

    if (argc < 2) {
        printf("\nUsage: mesh_generator [-v | -l | -s | -b | -t] [inputfile_A.off] [bool_opcode "
               "inputfile_B.off] [-e edges.off] [-p points.off]\n\n"
               "Defines the volume enclosed by the input OFF file(s) and saves a volume mesh to "
               "'volume.msh'. The surface (inputfile_A.off), -e edges, and -p points are all "
               "optional; provide any non-empty combination.\n\n"
               "Command line arguments:\n"
               "-v = verbose mode\n"
               "-l = logging mode (appends a line to mesh_generator.log)\n"
               "-s = save the mesh bounding surface to 'skin.off'\n"
               "-b = save the subdivided constraints to 'black_faces.off'\n"
               "-t = triangulate/tetrahedrize output\n"
               "-r = with -t, also write exact-rational vertex coords to 'volume.tet.rational'\n"
               "-2d = 2D mode: read an OBJ line soup (v/l records) and write the arrangement\n"
               "     of the input segments as a triangulation to 'triangulation.off'\n"
               "-a = keep all cells: skip the in/out min-cut and tetrahedralize the whole\n"
               "     domain, so the -t surface tracking is exact (no cells are deleted)\n"
               "-f = fill holes: cap open boundaries with ear-clipped triangles (off by\n"
               "     default; only affects the solid-output min-cut path)\n"
               "-e edges.off  = also insert these edges into the output (implies -a); OFF\n"
               "     with 'nv ne 0' then vertices then edge records '2 i j'\n"
               "-p points.off = also insert these points into the output (implies -a); OFF\n"
               "     with 'np 0 0' then vertex lines\n"
               "bool_opcode: {U, I, D}\n"
               "  U -> union (AuB),\n"
               "  I -> intersection (A^B),\n"
               "  D -> difference (A\\B)\n\n"
               "The parallel passes run on VOLUMEREMESHER_NUM_THREADS threads if set, else on\n"
               "the CPUs the process may run on.\n\n"
               "Example:\n"
               "mesh_generator ant.off U pig.off\n");
        return 0;
    }

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            // Multi-character flags must be matched in full BEFORE the single-character
            // dispatch below, which inspects only argv[i][1] and would read "-2d" as "-2".
            if (!strcmp(argv[i], "-2d")) {
                two_d = true;
                continue;
            }
            if (argv[i][1] == 'v')
                verbose = true;
            else if (argv[i][1] == 't')
                triangulate = true;
            else if (argv[i][1] == 'l')
                logging = true;
            else if (argv[i][1] == 'b')
                blackfaces = true;
            else if (argv[i][1] == 's')
                surfmesh = true;
            else if (argv[i][1] == 'r')
                export_rational = true;
            else if (argv[i][1] == 'f')
                fill_holes = true;
            else if (argv[i][1] == 'a')
                keep_all_cells = true;
            else if (argv[i][1] == 'e') {
                if (i + 1 >= argc) ip_error("-e requires an edge file argument\n");
                edge_file = argv[++i];
            } else if (argv[i][1] == 'p') {
                if (i + 1 >= argc) ip_error("-p requires a point file argument\n");
                point_file = argv[++i];
            } else
                ip_error("Unknown option\n");
        } else if (fileA_name == NULL)
            fileA_name = argv[i];
        else if (bool_opcode == '0')
            bool_opcode = argv[i][0];
        else if (fileB_name == NULL)
            fileB_name = argv[i];
        else
            ip_error("Too many args passed\n");
    }

    bool two_input = (bool_opcode != '0');

    // ---------------------------------------------------------------------------------------
    // 2D mode. A completely separate pipeline (src/VolumeRemesher/2d/): an OBJ line soup is
    // turned into the arrangement of its segments, triangulated, with per-segment provenance.
    // It shares nothing with the 3D path except the exact predicate kernel.
    // ---------------------------------------------------------------------------------------
    if (two_d) {
        if (fileA_name == NULL) ip_error("-2d needs an input .obj line soup\n");

        std::vector<double> coords;
        std::vector<uint32_t> indexes;
        if (!vol_rem::vr2d::read_OBJ_segments(fileA_name, coords, indexes))
            ip_error("Cannot read the input OBJ line soup\n");
        if (verbose)
            printf("Loaded %zu points and %zu segments from %s\n", coords.size() / 2,
                   indexes.size() / 2, fileA_name);

        const vr2d::ArrangementOutput arr = arrange_2d(coords, indexes, export_rational, verbose);

        if (!vr2d::write_arrangement("triangulation.off", arr))
            ip_error("Cannot write the output triangulation\n");

        if (verbose)
            printf("Wrote triangulation.off: %zu vertices, %zu triangles\n",
                   (size_t)arr.num_vertices(), arr.tris.size());
        return 0;
    }

    if (verbose) {
        if (fileB_name == NULL) {
            printf("\nResolve auto-intersections and/or repair.\n\n");
            printf("Loading %s\n", fileA_name);
        } else {
            printf("\nBoolean operator: ");
            if (bool_opcode == 'U')
                printf(" union.\n\n");
            else if (bool_opcode == 'I')
                printf(" intersection.\n\n");
            else if (bool_opcode == 'D')
                printf(" difference.\n\n");
            else {
                printf("INVALID\n\n");
                return 0;
            }
            printf("Loading %s and %s.\n\n", fileA_name, fileB_name);
        }
    }

    // The surface (fileA), edges (-e), and points (-p) are all optional; at least one
    // must be given. Only edges / only points / only a surface / any combination works.
    if (fileA_name == NULL && edge_file == NULL && point_file == NULL)
        ip_error("No input: provide a surface .off and/or -e edges.off and/or -p points.off\n");
    if (fileA_name == NULL && two_input)
        ip_error("A boolean operation needs a surface input\n");

    const bool file_B_is_tet =
        fileB_name != NULL && !strcmp(fileB_name + strlen(fileB_name) - 4, ".tet");
    const bool file_B_is_msh =
        fileB_name != NULL && !strcmp(fileB_name + strlen(fileB_name) - 4, ".msh");

    // A second input that is a tetrahedral mesh means "embed A into it" rather than a boolean.
    const bool embedsurf = (file_B_is_tet || file_B_is_msh);
    if (embedsurf) two_input = false;
    if (two_input && fileB_name == NULL) ip_error("A boolean operation needs a second surface\n");

    MeshOptions opts;
    opts.verbose = verbose;
    opts.logging = logging;
    opts.triangulate = triangulate;
    opts.export_rational = export_rational;
    opts.fill_holes = fill_holes;
    opts.keep_all_cells = keep_all_cells;
    // Black faces win over the skin when both are asked for, as they always did.
    opts.output = blackfaces ? MeshOutput::BLACK_FACES
                             : (surfmesh ? MeshOutput::SKIN : MeshOutput::VOLUME);

    TriMeshData A, B;
    if (fileA_name) A = read_OFF(fileA_name, verbose);
    if (two_input) B = read_OFF(fileB_name, verbose);

    // Optional extra edges/points to insert into the tetrahedralization.
    EdgeData edges;
    PointData points;
    if (edge_file) edges = read_edge_OFF(edge_file, verbose);
    if (point_file) points = read_points_OFF(point_file, verbose);

    MeshOutput out;
    if (embedsurf) {
        const TetMeshData T =
            file_B_is_tet ? read_TET(fileB_name, verbose) : read_MEDIT(fileB_name, verbose);
        out = embed_mesh(
            A, T, bool_opcode, edge_file ? &edges : nullptr, point_file ? &points : nullptr, opts);
    } else {
        out = make_mesh(
            fileA_name ? &A : nullptr,
            fileA_name,
            bool_opcode,
            two_input ? &B : nullptr,
            fileB_name,
            edge_file ? &edges : nullptr,
            point_file ? &points : nullptr,
            opts);
    }

    printf("Writing output file...\n");
    write_mesh_output(default_output_filename(opts), out);
    printf("Done.\n");

    return 0;
}
