#include "BSP.h"

using namespace vol_rem;

// Skip whitespace and '#'-comment lines (standard OFF allows them, e.g. files
// exported by meshio) so the parser lands on the next real token.
static void off_skip_comments(FILE* file)
{
    int c;
    while ((c = fgetc(file)) != EOF) {
        if (isspace(c))
            continue;
        if (c == '#') {
            while ((c = fgetc(file)) != EOF && c != '\n')
                ;
            continue;
        }
        ungetc(c, file);
        return;
    }
}

void read_OFF_file(
    const char* filename,
    double** vertices_p,
    uint32_t* npts,
    uint32_t** tri_vertices_p,
    uint32_t* ntri,
    bool verbose)
{
    FILE* file = fopen(filename, "r");
    if (file == NULL)
        ip_error("read_OFF_file: FATAL ERROR "
                 "cannot open input file.\n");

    // Check OFF mark (1st line).
    char file_ext_read[3];
    char file_ext_target[] = {'O', 'F', 'F'};
    if (fscanf(file, "%3c", file_ext_read) == 0)
        ip_error("read_OFF_file: FATAL ERROR "
                 "cannot read 1st line of input file\n");

    for (uint32_t i = 0; i < 3; i++)
        if (file_ext_read[i] != file_ext_target[i])
            ip_error("read_OFF_file: FATAL ERROR "
                     "1st line of input file is different from OFF\n");

    // Reading number of points and triangles (skipping any comment lines).
    off_skip_comments(file);
    if (fscanf(file, " %d %d %*d ", npts, ntri) == 0)
        ip_error("read_OFF_file: FATAL ERROR 2st line of "
                 "input file do not contanins point and triangles numbers.\n");

    if (verbose)
        printf(
            "file %s contains %d vertices and %d constraints (triangles)\n",
            filename,
            *npts,
            *ntri);

    // Reading points coordinates.
    *vertices_p = (double*)malloc(sizeof(double) * 3 * (*npts));
    *tri_vertices_p = (uint32_t*)malloc(sizeof(uint32_t) * 3 * (*ntri));

    for (uint32_t i = 0; i < (*npts); i++) {
        if (fscanf(
                file,
                " %lf %lf %lf ",
                (*vertices_p) + (i * 3),
                (*vertices_p) + (i * 3 + 1),
                (*vertices_p) + (i * 3 + 2)) == 0)
            ip_error("error reading input file\n");
    }

    uint32_t nv;
    for (uint32_t i = 0; i < (*ntri); i++) {
        if (fscanf(
                file,
                " %u %u %u %u ",
                &nv,
                (*tri_vertices_p) + (i * 3),
                (*tri_vertices_p) + (i * 3 + 1),
                (*tri_vertices_p) + (i * 3 + 2)) == 0)
            ip_error("error reading input file\n");
        if (nv != 3) ip_error("Non-triangular faces not supported\n");
    }
    fclose(file);
}

// Read an OFF-style edge file: "OFF", then "nv ne 0", nv vertex lines "x y z", then
// ne edge lines "2 i j" (like an OFF face record but with two vertex indices).
void read_edge_OFF(
    const char* filename,
    double** vertices_p,
    uint32_t* npts,
    uint32_t** edge_idx_p,
    uint32_t* nedges,
    bool verbose)
{
    FILE* file = fopen(filename, "r");
    if (file == NULL) ip_error("read_edge_OFF: FATAL ERROR cannot open edge file.\n");

    char ext[3];
    const char target[] = {'O', 'F', 'F'};
    if (fscanf(file, "%3c", ext) == 0) ip_error("read_edge_OFF: cannot read 1st line\n");
    for (uint32_t i = 0; i < 3; i++)
        if (ext[i] != target[i]) ip_error("read_edge_OFF: 1st line is not OFF\n");

    off_skip_comments(file);
    uint32_t dummy = 0;
    if (fscanf(file, " %u %u %u ", npts, nedges, &dummy) < 2)
        ip_error("read_edge_OFF: cannot read vertex/edge counts\n");
    if (verbose) printf("edge file %s: %u vertices, %u edges\n", filename, *npts, *nedges);

    *vertices_p = (double*)malloc(sizeof(double) * 3 * (*npts));
    *edge_idx_p = (uint32_t*)malloc(sizeof(uint32_t) * 2 * (*nedges));
    for (uint32_t i = 0; i < (*npts); i++)
        if (fscanf(
                file, " %lf %lf %lf ", (*vertices_p) + 3 * i, (*vertices_p) + 3 * i + 1,
                (*vertices_p) + 3 * i + 2) != 3)
            ip_error("read_edge_OFF: error reading vertex\n");
    uint32_t nv;
    for (uint32_t i = 0; i < (*nedges); i++) {
        if (fscanf(file, " %u %u %u ", &nv, (*edge_idx_p) + 2 * i, (*edge_idx_p) + 2 * i + 1) != 3)
            ip_error("read_edge_OFF: error reading edge\n");
        if (nv != 2) ip_error("read_edge_OFF: each edge record must list 2 vertices\n");
    }
    fclose(file);
}

// Read an OFF-style point file: "OFF", then "np 0 0", np vertex lines "x y z".
void read_points_OFF(const char* filename, double** vertices_p, uint32_t* npts, bool verbose)
{
    FILE* file = fopen(filename, "r");
    if (file == NULL) ip_error("read_points_OFF: FATAL ERROR cannot open point file.\n");

    char ext[3];
    const char target[] = {'O', 'F', 'F'};
    if (fscanf(file, "%3c", ext) == 0) ip_error("read_points_OFF: cannot read 1st line\n");
    for (uint32_t i = 0; i < 3; i++)
        if (ext[i] != target[i]) ip_error("read_points_OFF: 1st line is not OFF\n");

    off_skip_comments(file);
    uint32_t nf = 0, dummy = 0;
    if (fscanf(file, " %u %u %u ", npts, &nf, &dummy) < 1)
        ip_error("read_points_OFF: cannot read vertex count\n");
    if (verbose) printf("point file %s: %u points\n", filename, *npts);

    *vertices_p = (double*)malloc(sizeof(double) * 3 * (*npts));
    for (uint32_t i = 0; i < (*npts); i++)
        if (fscanf(
                file, " %lf %lf %lf ", (*vertices_p) + 3 * i, (*vertices_p) + 3 * i + 1,
                (*vertices_p) + 3 * i + 2) != 3)
            ip_error("read_points_OFF: error reading point\n");
    fclose(file);
}

void read_TET_file(
    const char* filename,
    double** vertices_p,
    uint32_t* npts,
    uint32_t** tet_vertices_p,
    uint32_t* ntet,
    bool verbose)
{
    FILE* file = fopen(filename, "r");

    uint32_t nv, nt;
    if (!fscanf(file, "%d vertices\n", &nv))
        ip_error("read_TET_file: FATAL ERROR\ncannot read 1st line of input file\n");
    if (!fscanf(file, "%d tets\n", &nt))
        ip_error("read_TET_file: FATAL ERROR\ncannot read 2nd line of input file\n");

    if (verbose) printf("file %s contains %u vertices and %u tetrahedra\n", filename, nv, nt);

    double* v_crd = (double*)malloc(nv * sizeof(double) * 3);

    for (uint32_t i = 0; i < nv; i++) {
        double* crd = v_crd + i * 3;
        if (fscanf(file, "%lf %lf %lf\n", crd, crd + 1, crd + 2) != 3)
            ip_error("read_TET_file: FATAL ERROR\ncannot read vertex coordinates\n");
    }

    uint32_t* t_idx = (uint32_t*)malloc(nt * sizeof(uint32_t) * 4);

    for (uint32_t i = 0; i < nt; i++) {
        uint32_t* nodes = t_idx + i * 4;
        if (fscanf(file, "4 %u %u %u %u\n", nodes, nodes + 1, nodes + 2, nodes + 3) != 4)
            ip_error("read_TET_file: FATAL ERROR\ncannot read tet indexes\n");
    }

    fclose(file);

    *vertices_p = v_crd;
    *npts = nv;
    *tet_vertices_p = t_idx;
    *ntet = nt;
}

void read_MEDIT_file(
    const char* filename,
    double** vertices_p,
    uint32_t* npts,
    uint32_t** tet_vertices_p,
    uint32_t* ntet,
    bool verbose)
{
    FILE* file = fopen(filename, "r");

    char line[1024];
    do {
        if (!fscanf(file, "%s\n", line))
            ip_error("read_MEDIT_file: Could not find Vertices keyword\n");
    } while (strcmp(line, "Vertices"));

    int dummy, nv = 0;
    if (!fscanf(file, "%d\n", &nv)) ip_error("read_MEDIT_file: Could not read num vertices\n");
    printf("Reading %d vertices\n", nv);

    double* v_crd = (double*)malloc(nv * sizeof(double) * 3);

    for (int i = 0; i < nv; i++) {
        double* crd = v_crd + i * 3;
        if (fscanf(file, "%lf %lf %lf %d\n", crd, crd + 1, crd + 2, &dummy) != 4)
            ip_error("read_MEDIT_file: Could not read vertex coords\n");
    }

    int nt;
    if (!fscanf(file, "%s\n", line) || strcmp(line, "Tetrahedra"))
        ip_error("read_MEDIT_file: Could not find Tetrahedra keyword\n");
    if (!fscanf(file, "%d\n", &nt)) ip_error("read_MEDIT_file: Could not read num tet\n");
    printf("Reading %d tetrahedra\n", nt);

    uint32_t* t_idx = (uint32_t*)malloc(nt * sizeof(uint32_t) * 4);

    for (uint32_t i = 0; i < (uint32_t)nt; i++) {
        uint32_t* nodes = t_idx + i * 4;
        if (fscanf(file, "%d %d %d %d %d\n", nodes, nodes + 1, nodes + 2, nodes + 3, &dummy) != 5)
            ip_error("read_MEDIT_file: Could not read tet indexes\n");
    }
    for (uint32_t i = 0; i < (uint32_t)nt * 4; i++) t_idx[i]--; // Decrease all indexes by one

    fclose(file);

    if (verbose) printf("file %s contains %u vertices and %u tetrahedra\n", filename, nv, nt);

    *vertices_p = v_crd;
    *npts = nv;
    *tet_vertices_p = t_idx;
    *ntet = nt;
}

#include "embed.h"
#include "2d/arrangement2d.h"
#include "2d/io2d.h"

/// <summary>
/// Main function
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>

int main(int argc, char** argv)
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
#if defined _MSC_VER && !defined NDEBUG
        fileA_name = (char*)"..\\models\\Octocat.off";
        fileB_name = (char*)"..\\models\\Octocat.bg.tet";
        bool_opcode = 'U';
        verbose = true;
#else
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
               "Example:\n"
               "mesh_generator ant.off U pig.off\n");
        return 0;
#endif
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

        vol_rem::vr2d::Arrangement2D arr;
        if (!vol_rem::vr2d::build_arrangement(coords, indexes, arr, verbose))
            ip_error("Failed to build the 2D arrangement\n");

        if (!vol_rem::vr2d::write_arrangement("triangulation.off", arr, export_rational))
            ip_error("Cannot write the output triangulation\n");

        if (verbose) {
            uint32_t nt = 0;
            for (uint32_t t = 0; t < arr.num_triangles(); t++)
                if (arr.tri_is_finite(t)) nt++;
            printf("Wrote triangulation.off: %zu vertices, %u triangles\n", arr.V.size(), nt);
        }
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

    double *coords_A = NULL, *coords_B = NULL;
    uint32_t ncoords_A = 0, ncoords_B = 0;
    uint32_t *tri_idx_A = NULL, *tri_idx_B = NULL;
    uint32_t ntriidx_A = 0, ntriidx_B = 0;

    const bool file_B_is_tet =
        fileB_name != NULL && !strcmp(fileB_name + strlen(fileB_name) - 4, ".tet");
    const bool file_B_is_msh =
        fileB_name != NULL && !strcmp(fileB_name + strlen(fileB_name) - 4, ".msh");

    bool embedsurf = (file_B_is_tet || file_B_is_msh);
    if (embedsurf) two_input = false;

    if (fileA_name)
        read_OFF_file(fileA_name, &coords_A, &ncoords_A, &tri_idx_A, &ntriidx_A, verbose);
    if (two_input)
        read_OFF_file(fileB_name, &coords_B, &ncoords_B, &tri_idx_B, &ntriidx_B, verbose);

    BSPcomplex* complex;
    if (embedsurf) {
        if (file_B_is_tet)
            read_TET_file(fileB_name, &coords_B, &ncoords_B, &tri_idx_B, &ntriidx_B, verbose);
        else
            read_MEDIT_file(fileB_name, &coords_B, &ncoords_B, &tri_idx_B, &ntriidx_B, verbose);

        // const std::vector<double> tri_vrt_coords(coords_A, coords_A + ncoords_A * 3);
        // const std::vector<uint32_t> triangle_indexes(tri_idx_A, tri_idx_A + ntriidx_A * 3);
        // const std::vector<double> tet_vrt_coords(coords_B, coords_B + ncoords_B * 3);
        // const std::vector<uint32_t> tet_indexes(tri_idx_B, tri_idx_B + ntriidx_B * 4);
        // std::vector<bigrational> vertices;
        // std::vector<uint32_t> facets;
        // std::vector<uint32_t> cells;
        // std::vector<uint32_t> facets_on_input;

        // embed_tri_in_poly_mesh(
        //     tri_vrt_coords,
        //     triangle_indexes,
        //     tet_vrt_coords,
        //     tet_indexes,
        //     vertices,
        //     facets,
        //     cells,
        //     facets_on_input,
        //     verbose
        //);

        complex = remakePolyhedralMesh(
            coords_A,
            ncoords_A,
            tri_idx_A,
            ntriidx_A,
            coords_B,
            ncoords_B,
            tri_idx_B,
            ntriidx_B,
            verbose);
        free(coords_A);
        free(tri_idx_A);
    } else {
        // Optional extra edges/points to insert into the tetrahedralization.
        extra_features_t extra;
        double *edge_verts = NULL, *point_verts = NULL;
        uint32_t *edge_idx = NULL, n_edge_verts = 0, n_edges = 0, n_points = 0;
        if (edge_file) {
            read_edge_OFF(edge_file, &edge_verts, &n_edge_verts, &edge_idx, &n_edges, verbose);
            extra.edge_verts = edge_verts;
            extra.n_edge_verts = n_edge_verts;
            extra.edge_idx = edge_idx;
            extra.n_edges = n_edges;
        }
        if (point_file) {
            read_points_OFF(point_file, &point_verts, &n_points, verbose);
            extra.point_verts = point_verts;
            extra.n_points = n_points;
        }
        complex = makePolyhedralMesh(
            fileA_name,
            coords_A,
            ncoords_A,
            tri_idx_A,
            ntriidx_A,
            fileB_name,
            coords_B,
            ncoords_B,
            tri_idx_B,
            ntriidx_B,
            bool_opcode,
            true,
            verbose,
            logging,
            fill_holes,
            keep_all_cells,
            (edge_file || point_file) ? &extra : nullptr);
        if (edge_verts) free(edge_verts);
        if (edge_idx) free(edge_idx);
        if (point_verts) free(point_verts);
    }

    printf("Writing output file...\n");
    if (blackfaces)
        complex->saveBlackFaces("black_faces.off", triangulate);
    else if (surfmesh)
        complex->saveSkin("skin.off", bool_opcode, triangulate);
    else
        complex->saveMesh(
            (triangulate) ? ("volume.tet") : ("volume.msh"),
            bool_opcode,
            triangulate,
            export_rational);
    printf("Done.\n");

    // make/remakePolyhedralMesh return an owning pointer; the caller must free it.
    delete complex;

    return 0;
}
