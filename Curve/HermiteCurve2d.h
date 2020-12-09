#pragma once
#include "Beziercurve2d.h"
namespace MN
{
	class HermiteCurve2d
	{
		Vec2 x0,x1,v0,v1;
		
	public:
		HermiteCurve2d() = default;
		std::shared_ptr<BezierCurve2d> bezCurve;
	public:
		HermiteCurve2d(Vec2 x0,Vec2 v0,Vec2 x1,Vec2 v1);
		Vec2 evaluate(double u);
		Vec2 differentiate(double u, int order);
	};
}