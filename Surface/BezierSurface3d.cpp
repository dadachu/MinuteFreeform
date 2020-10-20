/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#include "BezierSurface3d.h"
#include "Eigen/Dense"
namespace MN {
	void BezierSurface3d::subdivideCpts(const std::vector<Vec3>& cpts, Real t, std::vector<Vec3>& lower, std::vector<Vec3>& upper) {
		// De Casteljou's algorithm
		int size = (int)cpts.size();
		lower.resize(size);
		upper.resize(size);

		Real t1 = 1 - t;

		std::vector<Vec3> copy = cpts;
		lower[0] = copy.front();
		upper[size - 1] = copy.back();

		int cnt = 1;
		for (int i = size - 1; i > 0; i--) {
			std::vector<Vec3> tmp;
			tmp.resize(i);
			for (int j = 0; j < i; j++)
				tmp[j] = copy[j] * t1 + copy[j + 1] * t;
			copy = tmp;

			lower[cnt] = copy.front();
			upper[size - cnt - 1] = copy.back();
			cnt++;
		}
	}
	void BezierSurface3d::uSubdivide(Real u, BezierSurface3d& lower, BezierSurface3d& upper) const {
		ControlPoints lowerCpts = cpts;
		ControlPoints upperCpts = cpts;
		int rowNum = (int)cpts.size();
		int colNum = (int)cpts[0].size();
		std::vector<Vec3> lowerCol;
		std::vector<Vec3> upperCol;
		for (int i = 0; i < colNum; i++) {
			std::vector<Vec3> col;
			col.resize(rowNum);
			for (int j = 0; j < rowNum; j++)
				col[j] = cpts[j][i];
			subdivideCpts(col, u, lowerCol, upperCol);

			for (int j = 0; j < rowNum; j++) {
				lowerCpts[j][i] = lowerCol[j];
				upperCpts[j][i] = upperCol[j];
			}
		}
		lower = BezierSurface3d::create(uDegree, vDegree, lowerCpts, false);
		upper = BezierSurface3d::create(uDegree, vDegree, upperCpts, false);
	}
	void BezierSurface3d::vSubdivide(Real v, BezierSurface3d& lower, BezierSurface3d& upper) const {
		ControlPoints lowerCpts = cpts;
		ControlPoints upperCpts = cpts;
		int rowNum = (int)cpts.size();
		int colNum = (int)cpts[0].size();
		std::vector<Vec3> lowerRow;
		std::vector<Vec3> upperRow;
		for (int i = 0; i < rowNum; i++) {
			std::vector<Vec3> row;
			row.resize(colNum);
			for (int j = 0; j < colNum; j++)
				row[j] = cpts[i][j];
			subdivideCpts(row, v, lowerRow, upperRow);

			for (int j = 0; j < colNum; j++) {
				lowerCpts[i][j] = lowerRow[j];
				upperCpts[i][j] = upperRow[j];
			}
		}
		lower = BezierSurface3d::create(uDegree, vDegree, lowerCpts, false);
		upper = BezierSurface3d::create(uDegree, vDegree, upperCpts, false);
	}
	BezierSurface3d::Ptr BezierSurface3d::subdivide(const Domain& uSubdomain, const Domain& vSubdomain) const {
		// Divide in u direction
		BezierSurface3d lower, upper;
		uSubdivide(uSubdomain.beg(), lower, upper);	// upper
		Real t = (uSubdomain.end() - uSubdomain.beg()) / (1.0 - uSubdomain.beg());
		upper.uSubdivide(t, lower, upper);			// lower

		// Divide in v direction
		lower.vSubdivide(vSubdomain.beg(), lower, upper);	// upper
		t = (vSubdomain.end() - vSubdomain.beg()) / (1.0 - vSubdomain.beg());
		upper.vSubdivide(t, lower, upper);			// lower

		lower.updateDerivMat();
		return std::make_shared<BezierSurface3d>(lower);
	}

	Freeform3ds::ControlPoints BezierSurface3d::degreeElevation(int dir, const ControlPoints& cpts)
	{
		ControlPoints newCpts;
		newCpts.clear();
		if (dir == 1) // v
		{
			for (int i = 0; i < cpts.size(); i++) {
				std::vector<MN::Vec3> line;
				for (int j = 0; j < cpts[0].size(); j++)
					line.push_back(cpts[i][j]);

				std::vector<Vec3> newLine;
				newLine.push_back(line[0]);
				for (int j = 1; j < line.size(); j++)
				{
					double u = double(j) / double(line.size() + 1);
					newLine.push_back(line[j - 1] * u + line[j] * (1.0 - u));
				}
				newLine.push_back(line[line.size() - 1]);
				newCpts.push_back(newLine);
			}
		}
		else // u
		{

			newCpts.resize(cpts.size() + 1);
			for (int i = 0; i < cpts.size() + 1; i++)
				newCpts[i].resize(cpts[0].size());
			for (int i = 0; i < cpts[0].size(); i++)
			{
				std::vector<Vec3> line;
				for (int j = 0; j < cpts.size(); j++)
					line.push_back(cpts[j][i]);

				std::vector<Vec3> newLine;
				newLine.push_back(line[0]);
				for (int j = 1; j < line.size(); j++)
				{
					double u = double(j) / double(line.size() + 1);
					newLine.push_back(line[j - 1] * u + line[j] * (1.0 - u));
				}
				newLine.push_back(line[line.size() - 1]);

				for (int j = 0; j < newLine.size(); j++)
					newCpts[j][i] = newLine[j];
			}


		}
		return newCpts;
		ControlPoints vCpts;
		vCpts.resize(cpts.size());
		for (int i = 0; i < vCpts.size(); i++)
			vCpts[i].resize(cpts[i].size() + 1);
		for (int i = 0; i < cpts.size(); i++)
		{
			vCpts[i][0] = cpts[i][0];
			for (int j = 1; j < cpts[i].size(); j++)
			{
				double u = double(j) / double(cpts.size() + 1);
				vCpts[i][j] = cpts[i][j - 1] * u + cpts[i][j] * (1.0 - u);
			}
			vCpts[i][vCpts[i].size() - 1] = cpts[i][cpts[i].size() - 1];
		}

		ControlPoints nCpts;
		nCpts.resize(cpts.size() + 1);
		for (int i = 0; i < cpts.size() + 1; i++)
			nCpts[i].resize(cpts.size() + 1);
		for (int i = 0; i < nCpts.size(); i++)
		{
			nCpts[0][i] = vCpts[0][i];
			for (int j = 1; j < vCpts.size(); j++)
			{
				double u = double(j) / double(cpts.size() + 1);
				nCpts[j][i] = vCpts[j - 1][i] * u + vCpts[j][i] * (1.0 - u);
			}
			nCpts[nCpts.size() - 1][i] = vCpts[vCpts.size() - 1][i];
		}

		return nCpts;

	}

	Vec2 BezierSurface3d::projection(Vec3 xyz)
	{
		int RES = 100;
		double minDist = 1E+10;
		Vec2 minUV;
		for (int i = 0; i <= RES; i++)
		{
			double u = double(i) / double(RES);
			for (int j = 0; j <= RES; j++)
			{
				double v = double(j) / double(RES);
				auto pt = evaluate(u, v);
				if (xyz.dist(pt) < minDist)
				{
					minUV[0] = u;
					minUV[1] = v;
					minDist = xyz.dist(pt);
				}
			}
		}
		return projection(xyz, minUV);
	}

	Vec2 BezierSurface3d::projection(Vec3 xyz, Vec2 initialGuess)
	{
		Vec2 uv = initialGuess;
		for (int i = 0; i < 10; i++) { //MAGIC NUMBER 100
			Vec3 p0 = evaluate(uv[0], uv[1]);
			Vec3 p0_norm = normal(uv[0], uv[1]);
			Vec3 q = xyz - (p0_norm) * ((xyz - p0).dot(p0_norm));
			Vec3 Su = differentiate(uv[0], uv[1], 1, 0);
			Vec3 Sv = differentiate(uv[0], uv[1], 0, 1);
			Eigen::Matrix2d A;
			Eigen::Vector2d B, X;
			A << Su.dot(Su), Su.dot(Sv), Su.dot(Sv), Sv.dot(Sv);
			B << (q - p0).dot(Su),(q - p0).dot(Sv);
			X = A.colPivHouseholderQr().solve(B);
			double delta_u = X[0];
			double delta_v = X[1];
			uv[0] += delta_u;
			uv[1] += delta_v;
		}
		return uv;
	}

	BezierSurface3d BezierSurface3d::create(int uDegree, int vDegree, const ControlPoints& cpts, bool buildMat) {
		BezierSurface3d surface;
		surface.setDomain(0, Domain::create(0, 1));
		surface.setDomain(1, Domain::create(0, 1));
		surface.setDegree(0, uDegree);
		surface.setDegree(1, vDegree);
		surface.setCpts(cpts);
		if (buildMat)
			surface.updateDerivMat();
		return surface;
	}
	BezierSurface3d::Ptr BezierSurface3d::createPtr(int uDegree, int vDegree, const ControlPoints& cpts, bool buildMat) {
		BezierSurface3d surface = create(uDegree, vDegree, cpts, buildMat);
		return std::make_shared<BezierSurface3d>(surface);
	}
	void BezierSurface3d::updateDerivMat() {
		int rowNum = (int)cpts.size();
		int colNum = (int)cpts[0].size();
		if (uDegree > 0) {
			// derivMatU
			derivMatU.resize(rowNum - 1);
			for (int i = 0; i < rowNum - 1; i++) {
				derivMatU[i].resize(colNum);
				for (int j = 0; j < colNum; j++) {
					derivMatU[i][j] = (cpts[i + 1][j] - cpts[i][j]) * uDegree;
				}
			}
		}
		if (vDegree > 0) {
			// derivMatV
			derivMatV.resize(rowNum);
			for (int i = 0; i < rowNum; i++) {
				derivMatV[i].resize(colNum - 1);
				for (int j = 0; j < colNum - 1; j++) {
					derivMatV[i][j] = (cpts[i][j + 1] - cpts[i][j]) * vDegree;
				}
			}
		}
		if (uDegree > 1) {
			// derivMatUU
			derivMatUU.resize(rowNum - 2);
			for (int i = 0; i < rowNum - 2; i++) {
				derivMatUU[i].resize(colNum);
				for (int j = 0; j < colNum; j++) {
					derivMatUU[i][j] = (derivMatU[i + 1][j] - derivMatU[i][j]) * (uDegree - 1);
				}
			}
		}
		if (uDegree > 0 && vDegree > 0) {
			// derivMatUV
			derivMatUV.resize(rowNum - 1);
			for (int i = 0; i < rowNum - 1; i++) {
				derivMatUV[i].resize(colNum - 1);
				for (int j = 0; j < colNum - 1; j++) {
					derivMatUV[i][j] = (derivMatV[i + 1][j] - derivMatV[i][j]) * uDegree;
				}
			}
		}
		if (vDegree > 1) {
			// derivMatVV
			derivMatVV.resize(rowNum);
			for (int i = 0; i < rowNum; i++) {
				derivMatVV[i].resize(colNum - 2);
				for (int j = 0; j < colNum - 2; j++) {
					derivMatVV[i][j] = (derivMatV[i][j + 1] - derivMatV[i][j]) * (vDegree - 1);
				}
			}
		}
		if (uDegree > 2) {
			// derivMatUUU
			derivMatUUU.resize(rowNum - 3);
			for (int i = 0; i < rowNum - 3; i++) {
				derivMatUUU[i].resize(colNum);
				for (int j = 0; j < colNum; j++) {
					derivMatUUU[i][j] = (derivMatUU[i + 1][j] - derivMatUU[i][j]) * (uDegree - 2);
				}
			}
		}
		if (uDegree > 1 && vDegree > 0) {
			// derivMatUUV
			derivMatUUV.resize(rowNum - 2);
			for (int i = 0; i < rowNum - 2; i++) {
				derivMatUUV[i].resize(colNum - 1);
				for (int j = 0; j < colNum - 1; j++) {
					derivMatUUV[i][j] = (derivMatUU[i][j + 1] - derivMatUU[i][j]) * vDegree;
				}
			}
		}
		if (uDegree > 0 && vDegree > 1) {
			// derivMatUVV
			derivMatUVV.resize(rowNum - 1);
			for (int i = 0; i < rowNum - 1; i++) {
				derivMatUVV[i].resize(colNum - 2);
				for (int j = 0; j < colNum - 2; j++) {
					derivMatUVV[i][j] = (derivMatVV[i + 1][j] - derivMatVV[i][j]) * uDegree;
				}
			}
		}
		if (vDegree > 2) {
			// derivMatVVV
			derivMatVVV.resize(rowNum);
			for (int i = 0; i < rowNum; i++) {
				derivMatVVV[i].resize(colNum - 3);
				for (int j = 0; j < colNum - 3; j++) {
					derivMatVVV[i][j] = (derivMatVV[i][j + 1] - derivMatVV[i][j]) * (vDegree - 2);
				}
			}
		}
	}
	Vec3 BezierSurface3d::evaluate(Real u, Real v) const {
		BasisVector uBasis, vBasis;
		Bezier::calBasisVector(u, uDegree, uBasis);
		Bezier::calBasisVector(v, vDegree, vBasis);
		return tensorProduct(uBasis, cpts, vBasis);
	}
	Vec3 BezierSurface3d::differentiate(Real u, Real v, int uOrder, int vOrder) const {
		BasisVector uBasis, vBasis;
		if (uOrder == 0 && vOrder == 0)
			return evaluate(u, v);
		else if (uOrder == 1 && vOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			return tensorProduct(uBasis, derivMatU, vBasis);
		}
		else if (uOrder == 0 && vOrder == 1) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			return tensorProduct(uBasis, derivMatV, vBasis);
		}
		else if (uOrder == 2 && vOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 2, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			return tensorProduct(uBasis, derivMatUU, vBasis);
		}
		else if (uOrder == 1 && vOrder == 1) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			return tensorProduct(uBasis, derivMatUV, vBasis);
		}
		else if (uOrder == 0 && vOrder == 2) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 2, vBasis);
			return tensorProduct(uBasis, derivMatVV, vBasis);
		}
		else if (uOrder == 3 && vOrder == 0) {
			Bezier::calBasisVector(u, uDegree - 3, uBasis);
			Bezier::calBasisVector(v, vDegree, vBasis);
			return tensorProduct(uBasis, derivMatUUU, vBasis);
		}
		else if (uOrder == 2 && vOrder == 1) {
			Bezier::calBasisVector(u, uDegree - 2, uBasis);
			Bezier::calBasisVector(v, vDegree - 1, vBasis);
			return tensorProduct(uBasis, derivMatUUV, vBasis);
		}
		else if (uOrder == 1 && vOrder == 2) {
			Bezier::calBasisVector(u, uDegree - 1, uBasis);
			Bezier::calBasisVector(v, vDegree - 2, vBasis);
			return tensorProduct(uBasis, derivMatUVV, vBasis);
		}
		else if (uOrder == 0 && vOrder == 3) {
			Bezier::calBasisVector(u, uDegree, uBasis);
			Bezier::calBasisVector(v, vDegree - 3, vBasis);
			return tensorProduct(uBasis, derivMatVVV, vBasis);
		}
		else
			throw(std::runtime_error("Bezier surface differentiation is only allowed up to 2nd derivatives"));
	}
}