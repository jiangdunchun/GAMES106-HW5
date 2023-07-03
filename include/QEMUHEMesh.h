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

    inline bool DoSimplification(SimplificationMode mode, int targetFaceCount) {
        QEM_DEBUG("DoSimplification(mode=%d, targetFaceCount=%d)", mode, targetFaceCount);
        // TODO: implement this

		std::map<V*, Eigen::Matrix4d> Qs;
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
		}

		struct pair
		{
			V *v0;
			V *v1;
			Eigen::Vector3d pos;
			double error;
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

			ret.Q = Q0 + Q1;
			Eigen::Matrix4d d_Qv = ret.Q;
			d_Qv.row(3) << 0.0, 0.0, 0.0, 1.0;

			Eigen::Matrix4d d_Qv_i;
			bool invertible;
			double determinant;
			d_Qv.computeInverseAndDetWithCheck(d_Qv_i, determinant, invertible);
			//if (invertible)
			//{
			//	Eigen::Vector4d vv4 = d_Qv_i * Eigen::Vector4d(0, 0, 0, 1);

			//	ret.pos = Eigen::Vector3d(vv4.x(), vv4.y(), vv4.z()) / vv4.w();
			//}
			//else
			//{
				ret.pos = (v0->pos + v1->pos) / 2;
			//}
			Eigen::Vector4d vv4 = Eigen::Vector4d(ret.pos.x(), ret.pos.y(), ret.pos.z(), 1);

			ret.error = vv4.transpose() * ret.Q * vv4;

			return ret;
		};

		double t = 0;
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

};
