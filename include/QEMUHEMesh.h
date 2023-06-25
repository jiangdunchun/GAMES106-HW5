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

		std::vector<Eigen::Matrix4d> Qs(heMesh.Vertices().size(), Eigen::Matrix4d::Zero());
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

			Qs[heMesh.Index(v0)] += kp;
			Qs[heMesh.Index(v1)] += kp;
			Qs[heMesh.Index(v2)] += kp;
		}

		struct pair
		{
			V *v0;
			V *v1;
		};

		double t = 0.0;

		std::vector<pair> pairs;
		for (size_t i = 0; i < heMesh.Vertices().size(); i++)
		{
			V *v0 = heMesh.Vertices()[i];
			for (size_t j = i + 1; j < heMesh.Vertices().size(); j++)
			{
				V *v1 = heMesh.Vertices()[j];
				
				if (v0->IsConnectedWith(v1) || (v0->pos - v1->pos).norm() < t)
				{
					pairs.push_back({v0, v1});
				}
			}
		}

		int min_i;
		double min_error = DBL_MAX;
		Eigen::Vector3d min_vv;
		for (size_t i = 0; i < pairs.size(); ++i)
		{
			V *v0 = pairs[i].v0;
			V *v1 = pairs[i].v1;

			Eigen::Matrix4d Q0 = Qs[heMesh.Index(v0)];
			Eigen::Matrix4d Q1 = Qs[heMesh.Index(v1)];

			Eigen::Matrix4d Qv = Q0 + Q1;
			Qv.row(3) << 0.0, 0.0, 0.0, 1.0;

			Eigen::Vector3d vv;

			Eigen::Matrix4d Qv_i;
			bool invertible;
			double determinant;
			Qv.computeInverseAndDetWithCheck(Qv_i, determinant, invertible);
			if (invertible)
			{
				Eigen::Vector4d vv4 = Qv_i * Eigen::Vector4d(0, 0, 0, 1);

				vv = Eigen::Vector3d(vv4.x(), vv4.y(), vv4.z()) / vv4.w();
			}
			else
			{
				vv = (v0->pos + v1->pos) / 2;
			}
			Eigen::Vector4d vv4 = Eigen::Vector4d(vv.x(), vv.y(), vv.z(), 1);
			
			double error = vv4.transpose() * Qv * vv4;

			if (error < min_error)
			{
				min_error = error;
				min_i = i;
				min_vv = vv;
			}
		}

		V *v0 = pairs[min_i].v0;
		V *v1 = pairs[min_i].v1;

		Eigen::Vector3d vv = min_vv;

		v0->pos = vv;
		if (v0->IsConnectedWith(v1))
		{
			std::set<P *> removed_p;
			std::set<E *> removed_e;
			std::set<H *> removed_h;
			for (auto adj_v = v1->AdjVertices().begin();
				adj_v != v1->AdjVertices().end();
				adj_v++)
			{
				if (v0 == *adj_v)
					continue;

				if (v0->IsConnectedWith(*adj_v))
				{
					E *edge0 = V::EdgeBetween(v0, *adj_v);

					P *face0 = nullptr;
					if (!edge0->HalfEdge()->IsOnBoundary() && edge0->HalfEdge()->Next()->End() == v1)
					{
						face0 = edge0->HalfEdge()->Polygon();
					}
					else if (!edge0->HalfEdge()->Pair()->IsOnBoundary() && edge0->HalfEdge()->Pair()->Next()->End() == v1)
					{
						face0 = edge0->HalfEdge()->Pair()->Polygon();
					}

					removed_p.insert(face0);

					H *h_0 = nullptr;
					H *h_1 = nullptr;
					for (auto h : face0->AdjHalfEdges())
					{
						removed_h.insert(h);
						if ((h->Origin() == v0 && h->End() == *adj_v) || (h->End() == v0 && h->Origin() == *adj_v))
						{
							h_0 = h;
						}
						if ((h->Origin() == v1 && h->End() == *adj_v) || (h->End() == v1 && h->Origin() == *adj_v))
						{
							h_1 = h;
						}
					}

					removed_e.insert(h_0->Edge());
					removed_e.insert(V::EdgeBetween(v0, v1));

					for (auto adj_h : v1->AdjEdges())
					{
						if (adj_h->HalfEdge()->Origin() == v1)
						{
							adj_h->HalfEdge()->SetOrigin(v0);
						}

						if (adj_h->HalfEdge()->Pair()->Origin() == v1)
						{
							adj_h->HalfEdge()->Pair()->SetOrigin(v0);
						}

						
					}
					h_0->Pair()->SetEdge(h_1->Pair()->Edge());
					h_0->Pair()->SetPair(h_1->Pair());
					h_1->Pair()->SetPair(h_0->Pair());
				}
				else
				{
					H *h_edge1 = V::HalfEdgeAlong(v1, *adj_v);
					if (h_edge1)
					{
						if (h_edge1->Origin() == v1)
						{
							h_edge1->SetOrigin(v0);
						}
						else
						{
							h_edge1->Next()->SetOrigin(v0);
						}
					}
					

					H *hr_edge1 = V::HalfEdgeAlong(*adj_v, v1);
					if (hr_edge1)
					{
						if (hr_edge1->Origin() == v1)
						{
							hr_edge1->SetOrigin(v0);
						}
						else
						{
							hr_edge1->Next()->SetOrigin(v0);
						}
					}
				}
			}

			heMesh.RemoveEdge(V::EdgeBetween(v0, v1));
			heMesh.RemoveVertex(v1);
		}
		else
		{

		}

        return true;
    }

};
