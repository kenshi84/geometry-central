#pragma once

#include "geometrycentral/surface/embedded_geometry_interface.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/intrinsic_geometry_interface.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/elementary_geometry.h"


// An implmementation of the Signpost datastructure from
//   > "Navigating Intrinsic Triangulations". Sharp, Soliman, and Crane. SIGGRAPH 2019


namespace geometrycentral {
namespace surface {


class SignpostIntrinsicTriangulation : public IntrinsicGeometryInterface {

public:
  // Construct an intrinsic triangulation which sits atop this input mesh. Initially, the input triangulation will
  // just be a copy of the input mesh.
  SignpostIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom);

  // Construct from serialized blob
  SignpostIntrinsicTriangulation(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& inputGeom, const std::string& serializedBlob);

  std::string toSerializedBlob() const;

  static std::unique_ptr<ManifoldSurfaceMesh> getIntrinsicMeshFromSerializedBlob(const std::string& serializedBlob);

  // ======================================================
  // ======== Core Members
  // ======================================================

  // The underlying surface on which the intrinsic triangulation has been constructed
  ManifoldSurfaceMesh& inputMesh;
  IntrinsicGeometryInterface& inputGeom;

  // The connectivity of the intrinsic triangulation
  // note that somewhat confusingly, there is a .mesh reference which points to this same mesh,
  // inherited from the geometry interface
  std::unique_ptr<ManifoldSurfaceMesh> intrinsicMesh;

  // The geometry of the intrinsic triangulation
  EdgeData<double> intrinsicEdgeLengths;            // length of each edge
  HalfedgeData<double> intrinsicHalfedgeDirections; // direction of each halfedge, in radians from [0, angleSum]
  VertexData<double> intrinsicVertexAngleSums;      // vertex cone angle sum

  // NOTE: To enable use to make efficient use of the surface tracers, this class always automatically updates the
  // halfedgeVectorsInVertex and halfedgeVectorsInFace geometry members. Could remove this requirement if we change the
  // way the tracer works.

  // Vertex locations for the intrinsic triangulation
  VertexData<SurfacePoint> vertexLocations;

  // Minor parameters
  double delaunayEPS = 1e-6; // epsilon to use when testing if an edge is Delaunay

  // Marked edges, which cannot be removed.
  // (set to an array which holds true if an edge is fixed, and should not be flipped)
  // A callback is automatically registered which will update this array as edge splits are performed, so if a marked
  // edge is split the two resulting edges will be marked.
  EdgeData<char> markedEdges;
  void setMarkedEdges(const EdgeData<char>& markedEdges);


  // ======================================================
  // ======== Queries & Accessors
  // ======================================================

  // Given a point on the input triangulation, returns the corresponding point on the intrinsic triangulation
  SurfacePoint equivalentPointOnIntrinsic(SurfacePoint pointOnInput);

  // Given a point on the intrinsic triangulation, returns the corresponding point on the input triangulation
  SurfacePoint equivalentPointOnInput(SurfacePoint pointOnIntrinsic);

  // Trace out the edges of the intrinsic triangulation along the surface of the input mesh.
  // Each path is ordered along edge.halfedge(), and includes both the start and end points
  EdgeData<std::vector<SurfacePoint>> traceEdges(bool trimEnd = true);
  std::vector<SurfacePoint> traceHalfedge(Halfedge he, bool trimEnd = true); // trace a single intrinsic halfedge

  // Given data defined on the intrinsic triangulation, samples it at the vertices of the input triangulation
  template <typename T>
  VertexData<T> sampleAtInput(const VertexData<T>& dataOnIntrinsic);

  // Returns true if the intrinsic triangulation (or edge) satisifies the intrinsic Delaunay criterion
  bool isDelaunay();
  bool isDelaunay(Edge e);

  bool isIntrinsicEdgeOriginal(Edge eIntrinsic, Edge* eInput_ptr = nullptr) const;

  bool isInputEdgePreserved(Edge eInput, Edge* eIntrinsic_ptr = nullptr) const;

  // Returns the smallest angle in the intrinsic triangulation, in degrees
  double minAngleDegrees();

  // ======================================================
  // ======== High-Level Mutators
  // ======================================================
  //
  // Call once to build a useful triangulation

  // Flips edge in the intrinsic triangulation until is satisfies teh intrinsic Delaunay criterion
  void flipToDelaunay(std::function<bool(Edge)> flippableTest = {});

  // Perform intrinsic Delaunay refinement the intrinsic triangulation until it simultaneously:
  //   - satisfies the intrinsic Delaunay criterion
  //   - has no angles smaller than `angleThreshDegrees` (values > 30 degrees may not terminate)
  //   - has no triangles larger than `circumradiusThresh`
  // Terminates no matter what after maxInsertions insertions (infinite by default)
  void delaunayRefine(double angleThreshDegrees = 25.,
                      double circumradiusThresh = std::numeric_limits<double>::infinity(),
                      size_t maxInsertions = INVALID_IND);


  // General version of intrinsic Delaunay refinement, taking a function which will be called
  // to determine if a triangle should be refined.
  // Will return only when all triangles pass this function, or maxInsertions is exceeded, so
  // be sure to chose arguments such that the function terminates.
  void delaunayRefine(const std::function<bool(Face)>& shouldRefine, size_t maxInsertions = INVALID_IND);


  // Split edges which bend along an underlying surface.
  //   - angleThreshDeg: threshold at which to split (interpret as angle between normals), 0 would mean never split
  //   - relaiveL
  void splitBentEdges(EmbeddedGeometryInterface& posGeom, double angleThreshDeg = 30, double relativeLengthEPS = 1e-4,
                      size_t maxInsertions = INVALID_IND);


  // ======================================================
  // ======== Low-Level Mutators
  // ======================================================
  //
  // Basic operations to modify the intrinsic triangulation
  // NOTE: The individual operations to not call refreshQuantities(), so you should call it if you want quantities
  // updated.

  // If the edge is not Delaunay, flip it. Returns true if flipped.
  bool flipEdgeIfNotDelaunay(Edge e);

  // If the edge can be flipped, flip it (must be combinatorially flippable and inside a convex quad). Returns true if
  // flipped.
  bool flipEdgeIfPossible(Edge e, double possibleEPS = 1e-6, bool checkOnly = false);

  // Insert a new vertex in to the intrinsic triangulation
  Vertex insertVertex(SurfacePoint newPositionOnIntrinsic);

  // Insert the circumcenter of a face in to the triangulation. Returns the newly created intrinsic vertex.
  Vertex insertCircumcenter(Face f);

  // Insert the barycenter of a face in to the triangulation. Returns the newly created intrinsic vertex.
  Vertex insertBarycenter(Face f);

  // Remove an (inserted) vertex from the triangulation.
  // Note: if something goes terribly (numerically?) wrong, will exit without removing the vertex.
  Face removeInsertedVertex(Vertex v);

  // Collapse an interior edge whose endpoints are also assumed to be interior.
  // Takes halfedge he as input in order to control which vertex to be deleted (he.vertex is kept).
  // Returns false if it fails (when he.twin.vertex (i.e. vertex to be deleted) is an original vertex,
  // or when the collapse is geometrically infeasible)
  bool collapseInteriorEdge(Halfedge he, bool checkOnly = false);

  // Reverse of collapseInteriorEdge: given {heA, heB} outgoing from a common non-boundary vertex, create
  // a new vertex at a position offset from the common vertex by traceVec, then create a pair of triangles
  // formed by the original vertex, the new vertex, and the tip of each of {heA, heB}.
  // Returns halfedge pointing from the existing vertex to the new vertex.
  // The existing vertex and the new vertex's halfedges are set to the returned halfedge and the twin of the returned halfedge, respectively.
  // The operation is deemed infeasible in some cases (e.g., traceVec is too large and exits the vertex's
  // one ring, the flattened triangle fan of faces between heA & heB is highly concave and gives rise to
  // negatively oriented face, etc) where the method throws.
  Halfedge splitVertexAlongTwoEdges(Halfedge heA, Halfedge heB, Vector2 traceVec);

  // Relocate an inserted vertex to a location inside its trimmed one-ring triangle fan. Returns true if relocated.
  // With checkOnly==true, only perform the validity check and quit without actually performing the relocation
  bool relocateInsertedVertex(Vertex v, SurfacePoint pointOnIntrinsic, bool checkOnly = false);

  // Calls ManifoldSurfaceMesh::switchHalfedgeSides, then updates halfedge directions
  Halfedge switchHalfedgeSides(Edge e);

  // Split a halfedge
  Halfedge splitEdge(Halfedge he, double tSplit);

  // ======================================================
  // ======== Callbacks
  // ======================================================
  //
  // Get called whenever mesh mutations occur. Register a callback by inserting it in to this list.
  //

  // edge E if flipped
  std::list<std::function<void(Edge)>> edgeFlipCallbackList;

  // old face F is split by new vertex V
  std::list<std::function<void(Face, Vertex)>> faceInsertionCallbackList;

  // old edge E is split to halfedge HE1,HE2 both with he.vertex() as split vertex
  std::list<std::function<void(Edge, Halfedge, Halfedge)>> edgeSplitCallbackList;


  // ======================================================
  // ======== Geometry Immediates
  // ======================================================

  // Computed on the intrinsic triangulation
  double cornerAngle(Corner c) const;
  double halfedgeCotanWeight(Halfedge he) const;
  double edgeCotanWeight(Edge e) const;
  double area(Face f) const;
  double shortestEdge(Face f) const;
  double circumradius(Face f) const;


private:
  // ======================================================
  // ======== Geometry Interface
  // ======================================================
  //
  // Satisfy the requirements of the IntrinsicGeometryInterface

  // Override the compute edge lengths method from intrinsic geometry.
  virtual void computeEdgeLengths() override;

  // Override the halfedge vectors method from intrinsic geometry
  virtual void computeHalfedgeVectorsInVertex() override;


  // ======================================================
  // ======== Helpers
  // ======================================================

  // Insertion helpers
  Vertex insertVertex_face(SurfacePoint newPositionOnIntrinsic);
  Halfedge insertVertex_edge(SurfacePoint newPositionOnIntrinsic);
  void resolveNewVertex(Vertex newV, SurfacePoint intrinsicPoint);

  // Update a signpost angle from the (counter-)clockwise neighboring angle
  void updateAngleFromCWNeighor(Halfedge he);

  // Map angle to range [0, angleSum)
  double standardizeAngle(Vertex vert, double angle) const;

  // Get vector in rescaled vertex coordinates
  Vector2 halfedgeVector(Halfedge he) const;

  // Scale factor to take Euclidean data to cone data
  double vertexAngleScaling(Vertex v) const;

  // Repopulate the member halfedgeVectorInFace
  void updateFaceBasis(Face f);

  // Is this a marked or boundary edge?
  bool isFixed(Edge e);
  bool isOnFixedEdge(Vertex v); // boundary vertex or on fixed edge

  // Isometrically lay out the vertices around a halfedge in 2D coordinates
  // he points from vertex 2 to 0; others are numbered CCW
  std::array<Vector2, 4> layoutDiamond(Halfedge he);
  std::array<Vector2, 3> vertexCoordinatesInTriangle(Face face);

  // Callback helpers
  void invokeEdgeFlipCallbacks(Edge e);
  void invokeFaceInsertionCallbacks(Face f, Vertex v);
  void invokeEdgeSplitCallbacks(Edge e, Halfedge he1, Halfedge he2);
};


} // namespace surface
} // namespace geometrycentral


#include "geometrycentral/surface/signpost_intrinsic_triangulation.ipp"

