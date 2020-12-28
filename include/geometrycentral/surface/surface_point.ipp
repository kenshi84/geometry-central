#pragma once

namespace geometrycentral {
namespace surface {


// == Constructors
inline SurfacePoint::SurfacePoint() : type(SurfacePointType::Vertex) {}
inline SurfacePoint::SurfacePoint(Vertex v) : type(SurfacePointType::Vertex), vertex(v) {}
inline SurfacePoint::SurfacePoint(Edge e, double tEdge_) : type(SurfacePointType::Edge), edge(e), tEdge(tEdge_) {}
inline SurfacePoint::SurfacePoint(Halfedge he, double tHalfedge)
    : type(SurfacePointType::Edge), edge(he.edge()), tEdge(he == he.edge().halfedge() ? tHalfedge : (1. - tHalfedge)) {}
inline SurfacePoint::SurfacePoint(Face f, Vector3 faceCoords_)
    : type(SurfacePointType::Face), face(f), faceCoords(faceCoords_) {}


// === Serialization
template <class Archive>
void SurfacePoint::serialize(Archive& ar) {
  ar(type);

  if (type == SurfacePointType::Vertex) {
    ar(vertex);
  } else if (type == SurfacePointType::Edge) {
    ar(edge, tEdge);
  } else {
    ar(face, faceCoords);
  }
}

inline SurfaceMesh* SurfacePoint::getMesh() const {
  if (type == SurfacePointType::Vertex)
    return vertex.getMesh();
  else if (type == SurfacePointType::Edge)
    return edge.getMesh();
  else
    return face.getMesh();
}

inline void SurfacePoint::setMesh(SurfaceMesh* mesh) {
  if (type == SurfacePointType::Vertex)
    vertex.setMesh(mesh);
  else if (type == SurfacePointType::Edge)
    edge.setMesh(mesh);
  else
    face.setMesh(mesh);
}

// == Methods

inline std::ostream& operator<<(std::ostream& output, const SurfacePoint& p) {
  switch (p.type) {
  case SurfacePointType::Vertex: {
    output << "[SurfacePoint: type=Vertex, vertex= " << p.vertex << "]";
    break;
  }
  case SurfacePointType::Edge: {
    output << "[SurfacePoint: type=Edge, edge= " << p.edge << " tEdge= " << p.tEdge << "]";
    break;
  }
  case SurfacePointType::Face: {
    output << "[SurfacePoint: type=Face, face= " << p.face << " faceCoords= " << p.faceCoords << "]";
    break;
  }
  }

  return output;
}


inline SurfacePoint SurfacePoint::inSomeFace() const {

  switch (type) {
  case SurfacePointType::Vertex: {

    Halfedge he = vertex.halfedge();
    Face inFace = he.face();
    Halfedge targetHe = inFace.halfedge();

    // Find the appropriate barycentric coordinate and return
    if (he == targetHe) {
      return SurfacePoint(inFace, Vector3{1., 0., 0.});
    }
    he = he.next();
    if (he == targetHe) {
      return SurfacePoint(inFace, Vector3{0., 0., 1.});
    }
    return SurfacePoint(inFace, Vector3{0., 1., 0.});

    break;
  }
  case SurfacePointType::Edge: {

    Halfedge he = edge.halfedge();
    Face inFace = he.face();
    Halfedge targetHe = inFace.halfedge();

    // Find the appropriate barycentric coordinate and return
    if (he == targetHe) {
      return SurfacePoint(inFace, Vector3{1. - tEdge, tEdge, 0.});
    }
    he = he.next();
    if (he == targetHe) {
      return SurfacePoint(inFace, Vector3{tEdge, 0., 1. - tEdge});
    }
    return SurfacePoint(inFace, Vector3{0., 1. - tEdge, tEdge});

    break;
  }
  case SurfacePointType::Face: {
    return *this;
    break;
  }
  }

  throw std::logic_error("bad switch");
  return *this;
}

inline SurfacePoint SurfacePoint::inFace(Face targetFace) const {
  GC_SAFETY_ASSERT(getMesh() && getMesh() == targetFace.getMesh(), "wrong argument");

  switch (type) {
  case SurfacePointType::Vertex: {

    Halfedge he = targetFace.halfedge();

    // Find the appropriate barycentric coordinate and return
    if (he.vertex() == vertex) {
      return SurfacePoint(targetFace, Vector3{1., 0., 0.});
    }
    he = he.next();
    if (he.vertex() == vertex) {
      return SurfacePoint(targetFace, Vector3{0., 1., 0.});
    }
    he = he.next();
    if (he.vertex() == vertex) {
      return SurfacePoint(targetFace, Vector3{0., 0., 1.});
    }

    break;
  }

  case SurfacePointType::Edge: {

    double thisT = tEdge;
    for (Halfedge targetHe : {edge.halfedge(), edge.halfedge().twin()}) {

      int i = 0;
      for (Halfedge he : targetFace.adjacentHalfedges()) {
        if (he == targetHe) {
          // Find the appropriate barycentric coordinate and return

          Vector3 bary = Vector3::zero();
          bary[i] = 1.0 - thisT;
          bary[(i + 1) % 3] = thisT;

          return SurfacePoint(targetFace, bary);
        }
        i++;
      }

      // Flip the point to be along the other halfedge
      thisT = 1. - thisT;
    }

    break;
  }

  case SurfacePointType::Face: {
    if (face == targetFace) {
      return *this;
    };
    break;
  }
  }

  throw std::logic_error("SurfacePoint " + std::to_string(*this) + " not adjacent to target face " +
                         std::to_string(targetFace));
  return *this;
}

inline SurfacePoint SurfacePoint::inEdge(Edge targetEdge) const {
  GC_SAFETY_ASSERT(getMesh() && getMesh() == targetEdge.getMesh(), "wrong argument");

  switch (type) {
  case SurfacePointType::Vertex: {

    Halfedge he = targetEdge.halfedge();

    // Find the appropriate tEdge and return
    if (he.vertex() == vertex) {
      return SurfacePoint(targetEdge, 0.);
    }
    he = he.twin();
    if (he.vertex() == vertex) {
      return SurfacePoint(targetEdge, 1.);
    }

    break;
  }

  case SurfacePointType::Edge:
    if (edge == targetEdge)
      return *this;
    break;

  case SurfacePointType::Face:
    throw std::logic_error("inEdge cannot be called on a face point");
    break;

  }

  throw std::logic_error("SurfacePoint " + std::to_string(*this) + " not adjacent to target edge " +
                         std::to_string(targetEdge));
  return *this;
}


inline Vertex SurfacePoint::nearestVertex(double* delta) const {

  switch (type) {
  case SurfacePointType::Vertex: {
    if (delta) *delta = 0.;
    return vertex;
    break;
  }
  case SurfacePointType::Edge: {
    if (tEdge < .5) {
      if (delta) *delta = tEdge;
      return edge.halfedge().vertex();
    }
    if (delta) *delta = 1. - tEdge;
    return edge.halfedge().twin().vertex();
    break;
  }
  case SurfacePointType::Face: {
    if (faceCoords.x >= faceCoords.y && faceCoords.x >= faceCoords.z) {
      if (delta) *delta = norm(Vector3{1., 0., 0.} - faceCoords);
      return face.halfedge().vertex();
    }
    if (faceCoords.y >= faceCoords.x && faceCoords.y >= faceCoords.z) {
      if (delta) *delta = norm(Vector3{0., 1., 0.} - faceCoords);
      return face.halfedge().next().vertex();
    }
    if (delta) *delta = norm(Vector3{0., 0., 1.} - faceCoords);
    return face.halfedge().next().next().vertex();
    break;
  }
  }

  throw std::logic_error("bad switch");
  return vertex;
}

inline SurfacePoint SurfacePoint::reduced() const {
  switch (type) {
  case SurfacePointType::Vertex: {
    return *this;
  }
  case SurfacePointType::Edge: {
    if (tEdge == 0. || tEdge == 1.)
      return SurfacePoint(nearestVertex());
    else
      return *this;
  }
  case SurfacePointType::Face: {
    if (faceCoords.x == 1. || faceCoords.y == 1. || faceCoords.z == 1.)
      return SurfacePoint(nearestVertex());
    else if (faceCoords.z == 0.) {
        Edge e = face.halfedge().edge();
        double tEdge = face == e.halfedge().face() ? faceCoords.y : faceCoords.x;
        return SurfacePoint(e, tEdge);
    } else if (faceCoords.x == 0.) {
        Edge e = face.halfedge().next().edge();
        double tEdge = face == e.halfedge().face() ? faceCoords.z : faceCoords.y;
        return SurfacePoint(e, tEdge);
    } else if (faceCoords.y == 0.) {
        Edge e = face.halfedge().next().next().edge();
        double tEdge = face == e.halfedge().face() ? faceCoords.x : faceCoords.z;
        return SurfacePoint(e, tEdge);
    } else {
      return *this;
    }
  }
  }

  throw std::logic_error("bad switch");
  return {};
}

template <typename T>
inline T SurfacePoint::interpolate(const VertexData<T>& data) const {

  switch (type) {
  case SurfacePointType::Vertex: {
    return data[vertex];
    break;
  }
  case SurfacePointType::Edge: {
    T valTail = data[edge.halfedge().vertex()];
    T valTip = data[edge.halfedge().twin().vertex()];
    return (1. - tEdge) * valTail + tEdge * valTip;
    break;
  }
  case SurfacePointType::Face: {
    T valA = data[face.halfedge().vertex()];
    T valB = data[face.halfedge().next().vertex()];
    T valC = data[face.halfedge().next().next().vertex()];

    return (faceCoords.x * valA) + (faceCoords.y * valB) + (faceCoords.z * valC);
    break;
  }
  }

  throw std::logic_error("bad switch");
  return data[vertex];
}

inline void SurfacePoint::validate() const {
  const double EPS = 1e-4;

  switch (type) {
  case SurfacePointType::Vertex: {
    if (vertex == Vertex()) throw std::logic_error("surface point with Type::Vertex has invalid vertex ref");
    break;
  }
  case SurfacePointType::Edge: {
    if (edge == Edge()) throw std::logic_error("surface point with Type::Edge has invalid edge ref");
    if (!std::isfinite(tEdge)) throw std::logic_error("surface point with Type::Edge has non-finite tEdge");
    if (tEdge < -EPS || tEdge > (1. + EPS))
      throw std::logic_error("surface point with Type::Edge has tEdge outside of [0,1] = " + std::to_string(tEdge));
    break;
  }
  case SurfacePointType::Face: {
    if (face == Face()) throw std::logic_error("surface point with Type::Face has invalid face ref");
    if (!isfinite(faceCoords)) throw std::logic_error("surface point with Type::Face has non-finite face coords");
    if (faceCoords.x < -EPS || faceCoords.y < -EPS || faceCoords.z < -EPS)
      throw std::logic_error("surface point with Type::Face has negative bary coord " + std::to_string(faceCoords));

    if ((faceCoords.x + faceCoords.y + faceCoords.z) > (1. + EPS))
      throw std::logic_error("surface point with Type::Face has bary coord that sum to > 1 " +
                             std::to_string(faceCoords));
    break;
  }
  }
}

inline bool SurfacePoint::operator==(const SurfacePoint& other) const {
  if (type != other.type) return false;
  switch (type) {
  case SurfacePointType::Vertex: {
    return vertex == other.vertex;
    break;
  }
  case SurfacePointType::Edge: {
    return edge == other.edge && tEdge == other.tEdge;
    break;
  }
  case SurfacePointType::Face: {
    return face == other.face && faceCoords == other.faceCoords;
    break;
  }
  }
  return false; // should never be reached
}
inline bool SurfacePoint::operator!=(const SurfacePoint& other) const { return !(*this == other); }


inline Face sharedFace(const SurfacePoint& pA, const SurfacePoint& pB) {

  switch (pA.type) {

  case SurfacePointType::Vertex:
    for (Face f : pA.vertex.adjacentFaces()) {
      if (checkAdjacent(SurfacePoint(f, Vector3::zero()), pB)) return f;
    }
    break;

  case SurfacePointType::Edge:

    if (checkAdjacent(SurfacePoint(pA.edge.halfedge().face(), Vector3::zero()), pB)) {
      return pA.edge.halfedge().face();
    }
    if (checkAdjacent(SurfacePoint(pA.edge.halfedge().twin().face(), Vector3::zero()), pB)) {
      return pA.edge.halfedge().twin().face();
    }
    break;

  case SurfacePointType::Face:
    if (checkAdjacent(pA, pB)) return pA.face;
    break;
  }

  // no shared face
  return Face();
}
} // namespace surface
} // namespace geometrycentral

namespace std {

inline std::string to_string(geometrycentral::surface::SurfacePoint p) {
  ostringstream output;
  output << p;
  return output.str();
}

template <> struct hash<geometrycentral::surface::SurfacePoint> {
  size_t operator()(const geometrycentral::surface::SurfacePoint& sp) const {
    size_t seed = 0;
    hash_combine(seed, sp.vertex);
    hash_combine(seed, sp.edge);
    hash_combine(seed, sp.face);
    hash_combine(seed, sp.tEdge);
    hash_combine(seed, sp.faceCoords);
    return seed;
  }
};

} // namespace std
