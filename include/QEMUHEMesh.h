#pragma once
#include "UHEMesh/HeMesh.h"
#include "QEMDebug.h"
#include <Eigen/Sparse>
#include "assimp_helper.h"

class QEM
{
public:
	enum SimplificationMode
	{
		Position,
		PositionNormal,
		PositionNormalUV
	};

	class V;
	class E;
	using TraitsVEP = Ubpa::HEMeshTraits_EmptyPH<V, E>;
	using P = TraitsVEP::P;
	using H = TraitsVEP::H;
	using HEMesh = Ubpa::HEMesh<TraitsVEP>;

	class V : public Ubpa::TVertex<TraitsVEP>
	{
	public:
		V(Eigen::Vector3d pos, Eigen::Vector3d normal, Eigen::Vector2d uv) :
		  pos(pos), normal(normal), uv(uv)
		{}

		V() :
		    pos(Eigen::Vector3d::Zero())
		{}
		Eigen::Vector3d pos;
		Eigen::Vector3d normal;
		Eigen::Vector2d uv;

        // TODO: add vertex associated QEM attributes here
	};

	class E : public Ubpa::TEdge<TraitsVEP>
	{
		// TODO: add edge associated QEM attributes here
	};

	using HEMesh = Ubpa::HEMesh<TraitsVEP>;

private:
	HEMesh heMesh;

public:
    inline void PrintMesh()
	{
		QEM_DEBUG("HEMesh statistics:");
		QEM_DEBUG(" - Validity: %s", heMesh.IsValid() ? "Valid" : "Invalid");
		QEM_DEBUG(" - TriMesh: %s", heMesh.IsTriMesh() ? "True" : "False");

		QEM_DEBUG(" - Vertices: %zd", heMesh.Vertices().size());
		for (auto v : heMesh.Vertices())
		{
			QEM_DEBUG("   - %zd (%zx): he=%zd pos=(%lf, %lf, %lf)",
				heMesh.Index(v), v, heMesh.Index(v->HalfEdge()), v->pos[0], v->pos[1], v->pos[2]);
		}

		QEM_DEBUG(" - Half edges: %zd", heMesh.HalfEdges().size());
		for (auto he : heMesh.HalfEdges())
		{
			QEM_DEBUG("   - %zd (%zx):", heMesh.Index(he), he);

			std::string polygonStr;
			if (he->Polygon() != nullptr)
			{
				polygonStr = std::to_string(heMesh.Index(he->Polygon()));
			}
			else
			{
				polygonStr = "(nullptr)";
			}
			
			QEM_DEBUG(
			    "     - origin=%zd, edge=%zd, polygon=%s",
			    heMesh.Index(he->Origin()),
			    heMesh.Index(he->Edge()),
			    polygonStr.c_str()
			);

			QEM_DEBUG(
			    "     - next=%zd, pair=%zd",
			    heMesh.Index(he->Next()),
			    heMesh.Index(he->Pair()));
		}

		QEM_DEBUG(" - Edges: %zd", heMesh.Edges().size());
		for (auto e : heMesh.Edges())
		{
			QEM_DEBUG("   - %zd (%zx): he=%zd", heMesh.Index(e), e, heMesh.Index(e->HalfEdge()));
		}
			
		QEM_DEBUG(" - Polygons: %zd", heMesh.Polygons().size());
		for (auto p : heMesh.Polygons())
		{
			QEM_DEBUG("   - %zd (%zx): he=%zd", heMesh.Index(p), p, heMesh.Index(p->HalfEdge()));
		}

		auto boundaries = heMesh.Boundaries();
		QEM_DEBUG(" - Bounaries: %zd", boundaries.size());
		for (auto b : heMesh.Boundaries())
		{
			std::string out = "";

			for (auto he : b->NextLoop())
			{
				out += std::to_string(heMesh.Index(he->Origin()));
				out += "-";
			}

			out += std::to_string(heMesh.Index(b->Origin()));
			QEM_DEBUG("   - %s", out.c_str());
		}
	}

public:
	inline void ImportMesh(const Mesh& mesh)
	{
		heMesh.Clear();

		for (size_t i = 0; i < mesh.V.rows(); i++)
		{
			heMesh.AddVertex(mesh.V.row(i), mesh.N.row(i), mesh.UV.row(i));
		}

		for (size_t fIndex = 0; fIndex < mesh.F.rows(); fIndex++)
		{
			std::array<HEMesh::H *, 3> heLoop;
			for (size_t i = 0; i < 3; i++)
			{
				size_t next = (i + 1) % 3;
				assert(mesh.F(fIndex, i) != mesh.F(fIndex, next));
				auto *u  = heMesh.Vertices()[mesh.F(fIndex, i)];
				auto *v  = heMesh.Vertices()[mesh.F(fIndex, next)];
				auto *he = u->HalfEdgeTo(v);
				if (!he)
				{
					he = heMesh.AddEdge(u, v)->HalfEdge();
				}
					
				heLoop[i] = he;
			}
			auto *p = heMesh.AddPolygon(heLoop);
			assert(p != nullptr);
		}
        
        // PrintMesh();
	}

	inline Mesh ExportMesh(bool delineateBounaries = false)
	{
		Mesh mesh;
		
		mesh.F.resize(heMesh.Polygons().size(), 3);
		mesh.V.resize(heMesh.Vertices().size(), 3);
		mesh.N.resize(heMesh.Vertices().size(), 3);
		mesh.UV.resize(heMesh.Vertices().size(), 2);

		for (size_t i = 0; i < heMesh.Vertices().size(); i++)
		{
			mesh.V.row(i) = heMesh.Vertices()[i]->pos;
			mesh.N.row(i) = heMesh.Vertices()[i]->normal;
			mesh.UV.row(i) = heMesh.Vertices()[i]->uv;
		}

		for (size_t i = 0; i < heMesh.Polygons().size(); i++)
		{
			auto tri = heMesh.Polygons()[i];
			auto triHeIt = tri->AdjHalfEdges().begin();

			for (int j = 0; j < 3; j++)
			{
				mesh.F(i, j) = heMesh.Index(triHeIt->Origin());
				triHeIt++;
			}

			assert(triHeIt == tri->AdjHalfEdges().end());
		}

		if (delineateBounaries)
		{
            QEM_INFO("Delineate boundary feature is used. Recommend to turn wireframe view off to avoid overlapping.");
			auto boundaries = heMesh.Boundaries();
			for (auto b : boundaries)
			{
				for (auto he : b->NextLoop())
				{
					mesh.overlayEdges.push_back(std::make_tuple(
					    he->Origin()->pos, he->Next()->Origin()->pos,
					    Eigen::RowVector3d(1.0, 0.0, 0.0)));

				}
			}
		}

		return mesh;
	}

public:
	inline bool QEMBasic(int targetFaceCount)
	{
		std::map<V *, Eigen::Matrix4d> Qs;
		for (size_t i = 0; i < heMesh.Vertices().size(); i++)
		{
			Qs[heMesh.Vertices()[i]] = Eigen::Matrix4d::Zero();
		}
		for (size_t i = 0; i < heMesh.Polygons().size(); i++)
		{
			P *triangle = heMesh.Polygons()[i];

			V *v0 = triangle->HalfEdge()->Origin();
			V *v1 = triangle->HalfEdge()->End();
			V *v2 = triangle->HalfEdge()->Next()->End();

			Eigen::Vector3d normal = ((v2->pos - v0->pos).cross(v1->pos - v0->pos)).normalized();

			double distance = -1 * v0->pos.dot(normal);

			Eigen::Vector4d p(normal.x(), normal.y(), normal.z(), distance);

			Eigen::Matrix4d kp = p * p.transpose();

			Qs[v0] += kp;
			Qs[v1] += kp;
			Qs[v2] += kp;

			std::vector<H *> h_pairs = {
			    V::HalfEdgeAlong(v1, v0),
			    V::HalfEdgeAlong(v2, v1),
			    V::HalfEdgeAlong(v0, v2)};
			for (auto h : h_pairs)
			{
				if (h->IsOnBoundary())
				{
					Eigen::Vector3d edge = h->Origin()->pos - h->End()->pos;

					Eigen::Vector3d p_normal = (edge.cross(normal)).normalized();

					double p_distance = -1 * h->Origin()->pos.dot(p_normal);

					Eigen::Vector4d p_p(p_normal.x(), p_normal.y(), p_normal.z(), p_distance);

					Eigen::Matrix4d p_kp = p_p * p_p.transpose();

					Qs[h->Origin()] += p_kp;
					Qs[h->End()] += p_kp;
				}
			}
		}

		struct pair
		{
			V              *v0;
			V              *v1;
			Eigen::Vector3d pos;
			double          error;
			Eigen::Matrix4d Q;

			bool operator()(pair &a, pair &b)
			{
				return a.error > b.error;
			}
		};

		auto createPair = [&](V *v0, V *v1) -> pair {
			pair ret;
			ret.v0 = v0;
			ret.v1 = v1;

			Eigen::Matrix4d Q0 = Qs[v0];
			Eigen::Matrix4d Q1 = Qs[v1];

			ret.Q                = Q0 + Q1;
			Eigen::Matrix4d d_Qv = ret.Q;
			d_Qv.row(3) << 0.0, 0.0, 0.0, 1.0;

			Eigen::Matrix4d d_Qv_i;
			bool            invertible;
			double          determinant;
			d_Qv.computeInverseAndDetWithCheck(d_Qv_i, determinant, invertible);
			if (invertible)
			{
				Eigen::Vector4d vv4 = d_Qv_i * Eigen::Vector4d(0, 0, 0, 1);

				ret.pos = Eigen::Vector3d(vv4.x(), vv4.y(), vv4.z()) / vv4.w();
			}
			else
			{
				ret.pos = (v0->pos + v1->pos) / 2;
			}
			Eigen::Vector4d vv4 = Eigen::Vector4d(ret.pos.x(), ret.pos.y(), ret.pos.z(), 1);

			ret.error = vv4.transpose() * ret.Q * vv4;

			return ret;
		};

		double                                             t = 0;
		std::priority_queue<pair, std::vector<pair>, pair> pairs;
		for (size_t i = 0; i < heMesh.Vertices().size(); i++)
		{
			V *v0 = heMesh.Vertices()[i];
			for (size_t j = i + 1; j < heMesh.Vertices().size(); j++)
			{
				V *v1 = heMesh.Vertices()[j];

				if (v0->IsConnectedWith(v1) || (v0->pos - v1->pos).norm() < t)
				{
					pairs.push(createPair(v0, v1));
				}
			}
		}

		std::set<V *> delete_vertices;
		while (heMesh.Polygons().size() > targetFaceCount)
		{
			auto p = pairs.top();
			pairs.pop();

			if (delete_vertices.find(p.v0) != delete_vertices.end() || delete_vertices.find(p.v1) != delete_vertices.end())
				continue;

			delete_vertices.insert(p.v0);
			delete_vertices.insert(p.v1);

			if (p.v0->IsConnectedWith(p.v1))
			{
				E *e01 = V::EdgeBetween(p.v0, p.v1);

				V *vn = heMesh.CollapseEdge(e01);

				vn->pos = p.pos;

				Qs[vn] = p.Q;

				auto adjs = vn->AdjVertices();
				for (auto adj : adjs)
				{
					pairs.push(createPair(vn, adj));
				}
			}
			else
			{
			}
		}

		return true;
	}

	double t = 0;
	enum ModelType
	{
		Geometry = 0,
		Texture2D = 1 << 0,
		Normals = 1 << 1
	};
	struct Quadric 
	{
		Eigen::MatrixXd A;
		Eigen::VectorXd b;
		double c;

		Quadric() {}
		Quadric(int d) {
			A = Eigen::MatrixXd::Zero(d, d);
			b = Eigen::VectorXd::Zero(d);
			c = 0;
		}
		Quadric(Eigen::MatrixXd A_, Eigen::VectorXd b_, double c_)
		{
			A = A_;
			b = b_;
			c = c_;
		}
		Quadric operator+(const Quadric &Q) const
		{
			Quadric ret;
			ret.A = A + Q.A;
			ret.b = b + Q.b;
			ret.c = c + Q.c;
			return ret;
		}
	};
	struct Pair
	{
		V              *v0;
		V              *v1;
		Eigen::VectorXd pos;
		double          error;
		Quadric         Q;

		bool operator()(Pair &a, Pair &b)
		{
			return a.error > b.error;
		}
	};
	inline bool QEMwithMaterial(char modelType, int targetFaceCount) {
		int d = 3;
		if (modelType & ModelType::Texture2D)
			d += 2;
		if (modelType & ModelType::Normals)
			d += 3;

		auto getVertex = [&](V *v) -> Eigen::VectorXd {
			Eigen::VectorXd ret(d);

			int i = 0;

			ret.segment(i, 3) = v->pos;
			i += 3;

			if (modelType & ModelType::Texture2D)
			{
				ret.segment(i, 2) = v->uv;
				i += 2;
			}

			if (modelType & ModelType::Normals)
			{
				ret.segment(i, 3) = v->normal;
				i += 3;
			}

			return ret;
		};
		auto setV = [&](V *v, Eigen::VectorXd vertex) {
			int i = 0;

			v->pos = vertex.segment(i, 3);
			i += 3;

			if (modelType & ModelType::Texture2D)
			{
				v->uv = vertex.segment(i, 2);
				i += 2;
			}

			if (modelType & ModelType::Normals)
			{
				v->normal = vertex.segment(i, 3);
				i += 3;
			}
		};
		auto calQ = [&](Eigen::VectorXd p, Eigen::VectorXd q, Eigen::VectorXd r) -> Quadric {
			Eigen::VectorXd e1 = q - p;
			e1 = e1 / e1.norm();
			Eigen::VectorXd e2 = r - p - e1.dot(r - p) * e1;
			e2 = e2 / e2.norm();

			Eigen::MatrixXd A = Eigen::MatrixXd::Identity(d, d) - e1 * e1.transpose() - e2 * e2.transpose();

			Eigen::VectorXd b = p.dot(e1) * e1 + p.dot(e2) * e2 - p;

			double c = p.dot(p) - pow(p.dot(e1), 2) - pow(p.dot(e2), 2);

			return Quadric(A, b, c);
		};

		std::map<V *, Quadric> Qs;
		for (size_t i = 0; i < heMesh.Vertices().size(); i++)
		{
			Qs[heMesh.Vertices()[i]] = Quadric(d);
		}
		for (size_t i = 0; i < heMesh.Polygons().size(); i++)
		{
			P *triangle = heMesh.Polygons()[i];

			V *v0 = triangle->HalfEdge()->Origin();
			V *v1 = triangle->HalfEdge()->End();
			V *v2 = triangle->HalfEdge()->Next()->End();

			Quadric Q = calQ(getVertex(v0),getVertex( v1), getVertex(v2));
			Qs[v0] = Qs[v0] + Q;
			Qs[v1] = Qs[v1] + Q;
			Qs[v2] = Qs[v2] + Q;

			Eigen::Vector3d normal = ((v2->pos - v0->pos).cross(v1->pos - v0->pos)).normalized();
			std::vector<H *> h_pairs = {
			    V::HalfEdgeAlong(v1, v0),
			    V::HalfEdgeAlong(v2, v1),
			    V::HalfEdgeAlong(v0, v2)};
			for (auto h : h_pairs)
			{
				if (h->IsOnBoundary())
				{
					V *origin = h->Origin();
					V *end = h->End();

					Eigen::VectorXd p = getVertex(origin);
					Eigen::VectorXd q = getVertex(end);
					Eigen::VectorXd r = q;
					r.segment(0, 3) += normal;

					Quadric pQ = calQ(p, q, r);

					Qs[origin] = Qs[origin] + pQ;
					Qs[end] = Qs[end] + pQ;
				}
			}
		}

		auto createPair = [&](V *v0, V *v1) -> Pair {
			Pair ret;
			ret.v0 = v0;
			ret.v1 = v1;

			Quadric Q0 = Qs[v0];
			Quadric Q1 = Qs[v1];

			ret.Q = Q0 + Q1;
			
			if (ret.Q.A.determinant() != 0)
			{
				ret.pos = -1 * ret.Q.A.inverse() * ret.Q.b;
			}
			else
			{
				ret.pos = (getVertex(v0) + getVertex(v1)) / 2;
			}
			ret.error = ret.Q.b.transpose() * ret.pos + ret.Q.c;

			return ret;
		};
		std::priority_queue<Pair, std::vector<Pair>, Pair> pairs;
		for (size_t i = 0; i < heMesh.Vertices().size(); i++)
		{
			V *v0 = heMesh.Vertices()[i];
			for (size_t j = i + 1; j < heMesh.Vertices().size(); j++)
			{
				V *v1 = heMesh.Vertices()[j];

				if (v0->IsConnectedWith(v1) || (v0->pos - v1->pos).norm() < t)
				{
					pairs.push(createPair(v0, v1));
				}
			}
		}

		std::set<V *> deleteVs;
		while (heMesh.Polygons().size() > targetFaceCount)
		{
			auto p = pairs.top();
			pairs.pop();

			if (deleteVs.find(p.v0) != deleteVs.end() || deleteVs.find(p.v1) != deleteVs.end())
				continue;

			if (p.v0->IsConnectedWith(p.v1))
			{
				E *e01 = V::EdgeBetween(p.v0, p.v1);

				if (!heMesh.IsCollapsable(e01)) 
					continue;

				V *vn = heMesh.CollapseEdge(e01);

				setV(vn, p.pos);

				Qs[vn] = p.Q;

				auto adjs = vn->AdjVertices();
				for (auto adj : adjs)
				{
					pairs.push(createPair(vn, adj));
				}
			}
			else
			{
			}

			deleteVs.insert(p.v0);
			deleteVs.insert(p.v1);
		}

		return true;
	}

    inline bool DoSimplification(SimplificationMode mode, int targetFaceCount) {
        QEM_DEBUG("DoSimplification(mode=%d, targetFaceCount=%d)", mode, targetFaceCount);
        
		char modelType = 0;
		if (mode == SimplificationMode::PositionNormal || mode == SimplificationMode::PositionNormalUV)
			modelType |= ModelType::Normals;
		if (mode == SimplificationMode::PositionNormalUV)
			modelType |= ModelType::Texture2D;

		return QEMwithMaterial(modelType, targetFaceCount);
    }
};
