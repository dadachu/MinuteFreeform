#include "HermiteCurve2d.h"

namespace MN
{
	HermiteCurve2d::HermiteCurve2d(Vec2 x0, Vec2 v0, Vec2 x1, Vec2 v1)
	{
		this->x0 = x0;
		this->v0  =v0;
		this->x1 = x1;
		this->v1 = v1;
		BezierCurve2d::ControlPoints cp;
		cp.push_back(x0);
		cp.push_back(x0+v0*(1.0/3.0));
		cp.push_back(x1-v1*(1.0/3.0));
		cp.push_back(x1);
		bezCurve = BezierCurve2d::createPtr(3,cp);
	}

	Vec2 HermiteCurve2d::evaluate(double u)
	{
		return bezCurve->evaluate(u);
	}
	Vec2 HermiteCurve2d::differentiate(double u, int order)
	{
		return bezCurve->differentiate(u,order);
	}

}
