// Patch_Chapter15.h

// fix InstancingAndCullingDemo.cpp to fit DirectXMath and DirectXCollision.h

// the demo is from the book: 3D Game Programming with DirectX 11, Chapter 15

// the functions is from xnacollision.cpp

// endrollex

// Changed to support Unicode. 2016-10-21 by.progsecu
// 출처 파악제한 패치 소스로부터 유니코드 지원하게 변경
#ifndef PATCH_CHAPTERS_H
#define PATCH_CHAPTERS_H

#include <DirectXCollision.h>
#include <Windows.h>
#include <cfloat>

//

#if defined(DEBUG) | defined(_DEBUG)

	#define XMASSERT(Expression) ((VOID)((Expression) || (XMAssert(#Expression, __FILEW__, __LINE__), 0)))

#else

	#define XMASSERT(Expression) ((VOID)0)

#endif

static VOID XMAssert(_In_z_ CONST CHAR* pExpression, _In_z_ CONST WCHAR* pFileName, UINT LineNumber)

{

	// I can not find this function in xnamath.h series files

	std::string stExp = pExpression;

	std::wstring wstExp(stExp.begin(), stExp.end());

	DXTrace(pFileName, LineNumber, false, wstExp.c_str(), true);

	return;

}

//
namespace XNA

{

	typedef DirectX::BoundingBox AxisAlignedBox;

	typedef DirectX::BoundingFrustum Frustum;

	typedef DirectX::BoundingOrientedBox OrientedBox;

// copy from xnacollision.cpp

static const XMVECTOR g_UnitQuaternionEpsilon =

{

	1.0e-4f, 1.0e-4f, 1.0e-4f, 1.0e-4f

};

//-----------------------------------------------------------------------------

// Return TRUE if the quaterion is a unit quaternion.

//-----------------------------------------------------------------------------

static inline BOOL XMQuaternionIsUnit( FXMVECTOR Q )

{

	XMVECTOR Difference = XMVector4Length( Q ) - XMVectorSplatOne();

	return XMVector4Less( XMVectorAbs( Difference ), g_UnitQuaternionEpsilon );

}

} // namespace

//-----------------------------------------------------------------------------

// Build a frustum from a persepective projection matrix.  The matrix may only

// contain a projection; any rotation, translation or scale will cause the

// constructed frustum to be incorrect.

//-----------------------------------------------------------------------------

static VOID ComputeFrustumFromProjection( BoundingFrustum* pOut, XMMATRIX* pProjection )

{

	XMASSERT( pOut );

	XMASSERT( pProjection );

	// Corners of the projection frustum in homogenous space.

	static XMVECTOR HomogenousPoints[6] =

	{

		{  1.0f,  0.0f, 1.0f, 1.0f },	// right (at far plane)

		{ -1.0f,  0.0f, 1.0f, 1.0f },	// left

		{  0.0f,  1.0f, 1.0f, 1.0f },	// top

		{  0.0f, -1.0f, 1.0f, 1.0f },	// bottom

		{ 0.0f, 0.0f, 0.0f, 1.0f }, 	// near

		{ 0.0f, 0.0f, 1.0f, 1.0f }		// far

	};

	XMVECTOR Determinant;

	XMMATRIX matInverse = XMMatrixInverse( &Determinant, *pProjection );

	// Compute the frustum corners in world space.

	XMVECTOR Points[6];

	for( INT i = 0; i < 6; i++ )

	{

		// Transform point.

		Points[i] = XMVector4Transform( HomogenousPoints[i], matInverse );

	}

	pOut->Origin = XMFLOAT3( 0.0f, 0.0f, 0.0f );

	pOut->Orientation = XMFLOAT4( 0.0f, 0.0f, 0.0f, 1.0f );

	// Compute the slopes.

	Points[0] = Points[0] * XMVectorReciprocal( XMVectorSplatZ( Points[0] ) );

	Points[1] = Points[1] * XMVectorReciprocal( XMVectorSplatZ( Points[1] ) );

	Points[2] = Points[2] * XMVectorReciprocal( XMVectorSplatZ( Points[2] ) );

	Points[3] = Points[3] * XMVectorReciprocal( XMVectorSplatZ( Points[3] ) );

	pOut->RightSlope = XMVectorGetX( Points[0] );

	pOut->LeftSlope = XMVectorGetX( Points[1] );

	pOut->TopSlope = XMVectorGetY( Points[2] );

	pOut->BottomSlope = XMVectorGetY( Points[3] );

	// Compute near and far.

	Points[4] = Points[4] * XMVectorReciprocal( XMVectorSplatW( Points[4] ) );

	Points[5] = Points[5] * XMVectorReciprocal( XMVectorSplatW( Points[5] ) );

	pOut->Near = XMVectorGetZ( Points[4] );

	pOut->Far = XMVectorGetZ( Points[5] );

	return;

}

namespace XNA

{

//-----------------------------------------------------------------------------

// Transform a frustum by an angle preserving transform.

//-----------------------------------------------------------------------------

static VOID TransformFrustum( Frustum* pOut, const Frustum* pIn, FLOAT Scale, FXMVECTOR Rotation, FXMVECTOR Translation )

{

	XMASSERT( pOut );

	XMASSERT( pIn );

	XMASSERT( XMQuaternionIsUnit( Rotation ) );

	// Load the frustum.

	XMVECTOR Origin = XMLoadFloat3( &pIn->Origin );

	XMVECTOR Orientation = XMLoadFloat4( &pIn->Orientation );

	XMASSERT( XMQuaternionIsUnit( Orientation ) );

	// Composite the frustum rotation and the transform rotation.

	Orientation = XMQuaternionMultiply( Orientation, Rotation );

	// Transform the origin.

	Origin = XMVector3Rotate( Origin * XMVectorReplicate( Scale ), Rotation ) + Translation;

	// Store the frustum.

	XMStoreFloat3( &pOut->Origin, Origin );

	XMStoreFloat4( &pOut->Orientation, Orientation );

	// Scale the near and far distances (the slopes remain the same).

	pOut->Near = pIn->Near * Scale;

	pOut->Far = pIn->Far * Scale;

	// Copy the slopes.

	pOut->RightSlope = pIn->RightSlope;

	pOut->LeftSlope = pIn->LeftSlope;

	pOut->TopSlope = pIn->TopSlope;

	pOut->BottomSlope = pIn->BottomSlope;

	return;

}

//-----------------------------------------------------------------------------

// Return TRUE if any of the elements of a 3 vector are equal to 0xffffffff.

// Slightly more efficient than using XMVector3EqualInt.

//-----------------------------------------------------------------------------

static inline BOOL XMVector3AnyTrue( FXMVECTOR V )

{

	XMVECTOR C;

	// Duplicate the fourth element from the first element.

	C = XMVectorSwizzle( V, 0, 1, 2, 0 );

	return XMComparisonAnyTrue( XMVector4EqualIntR( C, XMVectorTrueInt() ) );

}

//-----------------------------------------------------------------------------

// Exact oriented box vs frustum test.

// Return values: 0 = no intersection,

//				  1 = intersection,

//				  2 = box is completely inside frustum

//-----------------------------------------------------------------------------

static INT IntersectOrientedBoxFrustum( const OrientedBox* pVolumeA, const Frustum* pVolumeB )

{

	XMASSERT( pVolumeA );

	XMASSERT( pVolumeB );

	static const XMVECTORI32 SelectY =

	{

		XM_SELECT_0, XM_SELECT_1, XM_SELECT_0, XM_SELECT_0

	};

	static const XMVECTORI32 SelectZ =

	{

		XM_SELECT_0, XM_SELECT_0, XM_SELECT_1, XM_SELECT_0

	};

	XMVECTOR Zero = XMVectorZero();

	// Build the frustum planes.

	XMVECTOR Planes[6];

	Planes[0] = XMVectorSet( 0.0f, 0.0f, -1.0f, pVolumeB->Near );

	Planes[1] = XMVectorSet( 0.0f, 0.0f, 1.0f, -pVolumeB->Far );

	Planes[2] = XMVectorSet( 1.0f, 0.0f, -pVolumeB->RightSlope, 0.0f );

	Planes[3] = XMVectorSet( -1.0f, 0.0f, pVolumeB->LeftSlope, 0.0f );

	Planes[4] = XMVectorSet( 0.0f, 1.0f, -pVolumeB->TopSlope, 0.0f );

	Planes[5] = XMVectorSet( 0.0f, -1.0f, pVolumeB->BottomSlope, 0.0f );

	// Load origin and orientation of the frustum.

	XMVECTOR Origin = XMLoadFloat3( &pVolumeB->Origin );

	XMVECTOR FrustumOrientation = XMLoadFloat4( &pVolumeB->Orientation );

	XMASSERT( XMQuaternionIsUnit( FrustumOrientation ) );

	// Load the box.

	XMVECTOR Center = XMLoadFloat3( &pVolumeA->Center );

	XMVECTOR Extents = XMLoadFloat3( &pVolumeA->Extents );

	XMVECTOR BoxOrientation = XMLoadFloat4( &pVolumeA->Orientation );

	XMASSERT( XMQuaternionIsUnit( BoxOrientation ) );

	// Transform the oriented box into the space of the frustum in order to

	// minimize the number of transforms we have to do.

	Center = XMVector3InverseRotate( Center - Origin, FrustumOrientation );

	BoxOrientation = XMQuaternionMultiply( BoxOrientation, XMQuaternionConjugate( FrustumOrientation ) );

	// Set w of the center to one so we can dot4 with the plane.

	Center = XMVectorInsert( Center, XMVectorSplatOne(), 0, 0, 0, 0, 1);

	// Build the 3x3 rotation matrix that defines the box axes.

	XMMATRIX R = XMMatrixRotationQuaternion( BoxOrientation );

	// Check against each plane of the frustum.

	XMVECTOR Outside = XMVectorFalseInt();

	XMVECTOR InsideAll = XMVectorTrueInt();

	XMVECTOR CenterInsideAll = XMVectorTrueInt();

	for( INT i = 0; i < 6; i++ )

	{

		// Compute the distance to the center of the box.

		XMVECTOR Dist = XMVector4Dot( Center, Planes[i] );

		// Project the axes of the box onto the normal of the plane.  Half the

		// length of the projection (sometime called the "radius") is equal to

		// h(u) * abs(n dot b(u))) + h(v) * abs(n dot b(v)) + h(w) * abs(n dot b(w))

		// where h(i) are extents of the box, n is the plane normal, and b(i) are the

		// axes of the box.

		XMVECTOR Radius = XMVector3Dot( Planes[i], R.r[0] );

		Radius = XMVectorSelect( Radius, XMVector3Dot( Planes[i], R.r[1] ), SelectY );

		Radius = XMVectorSelect( Radius, XMVector3Dot( Planes[i], R.r[2] ), SelectZ );

		Radius = XMVector3Dot( Extents, XMVectorAbs( Radius ) );

		// Outside the plane?

		Outside = XMVectorOrInt( Outside, XMVectorGreater( Dist, Radius ) );

		// Fully inside the plane?

		InsideAll = XMVectorAndInt( InsideAll, XMVectorLessOrEqual( Dist, -Radius ) );

		// Check if the center is inside the plane.

		CenterInsideAll = XMVectorAndInt( CenterInsideAll, XMVectorLessOrEqual( Dist, Zero ) );

	}

	// If the box is outside any of the planes it is outside.

	if ( XMVector4EqualInt( Outside, XMVectorTrueInt() ) )

		return 0;

	// If the box is inside all planes it is fully inside.

	if ( XMVector4EqualInt( InsideAll, XMVectorTrueInt() ) )

		return 2;

	// If the center of the box is inside all planes and the box intersects

	// one or more planes then it must intersect.

	if ( XMVector4EqualInt( CenterInsideAll, XMVectorTrueInt() ) )

		return 1;

	// Build the corners of the frustum.

	XMVECTOR RightTop = XMVectorSet( pVolumeB->RightSlope, pVolumeB->TopSlope, 1.0f, 0.0f );

	XMVECTOR RightBottom = XMVectorSet( pVolumeB->RightSlope, pVolumeB->BottomSlope, 1.0f, 0.0f );

	XMVECTOR LeftTop = XMVectorSet( pVolumeB->LeftSlope, pVolumeB->TopSlope, 1.0f, 0.0f );

	XMVECTOR LeftBottom = XMVectorSet( pVolumeB->LeftSlope, pVolumeB->BottomSlope, 1.0f, 0.0f );

	XMVECTOR Near = XMVectorReplicatePtr( &pVolumeB->Near );

	XMVECTOR Far = XMVectorReplicatePtr( &pVolumeB->Far );

	XMVECTOR Corners[8];

	Corners[0] = RightTop * Near;

	Corners[1] = RightBottom * Near;

	Corners[2] = LeftTop * Near;

	Corners[3] = LeftBottom * Near;

	Corners[4] = RightTop * Far;

	Corners[5] = RightBottom * Far;

	Corners[6] = LeftTop * Far;

	Corners[7] = LeftBottom * Far;

	// Test against box axes (3)

	{

		// Find the min/max values of the projection of the frustum onto each axis.

		XMVECTOR FrustumMin, FrustumMax;

		FrustumMin = XMVector3Dot( Corners[0], R.r[0] );

		FrustumMin = XMVectorSelect( FrustumMin, XMVector3Dot( Corners[0], R.r[1] ), SelectY );

		FrustumMin = XMVectorSelect( FrustumMin, XMVector3Dot( Corners[0], R.r[2] ), SelectZ );

		FrustumMax = FrustumMin;

		for( INT i = 1; i < 8; i++ )

		{

			XMVECTOR Temp = XMVector3Dot( Corners[i], R.r[0] );

			Temp = XMVectorSelect( Temp, XMVector3Dot( Corners[i], R.r[1] ), SelectY );

			Temp = XMVectorSelect( Temp, XMVector3Dot( Corners[i], R.r[2] ), SelectZ );

			FrustumMin = XMVectorMin( FrustumMin, Temp );

			FrustumMax = XMVectorMax( FrustumMax, Temp );

		}

		// Project the center of the box onto the axes.

		XMVECTOR BoxDist = XMVector3Dot( Center, R.r[0] );

		BoxDist = XMVectorSelect( BoxDist, XMVector3Dot( Center, R.r[1] ), SelectY );

		BoxDist = XMVectorSelect( BoxDist, XMVector3Dot( Center, R.r[2] ), SelectZ );

		// The projection of the box onto the axis is just its Center and Extents.

		// if (min > box_max || max < box_min) reject;

		XMVECTOR Result = XMVectorOrInt( XMVectorGreater( FrustumMin, BoxDist + Extents ),

										  XMVectorLess( FrustumMax, BoxDist - Extents ) );

		if( XMVector3AnyTrue( Result ) )

			return 0;

	}

	// Test against edge/edge axes (3*6).

	XMVECTOR FrustumEdgeAxis[6];

	FrustumEdgeAxis[0] = RightTop;

	FrustumEdgeAxis[1] = RightBottom;

	FrustumEdgeAxis[2] = LeftTop;

	FrustumEdgeAxis[3] = LeftBottom;

	FrustumEdgeAxis[4] = RightTop - LeftTop;

	FrustumEdgeAxis[5] = LeftBottom - LeftTop;

	for( INT i = 0; i < 3; i++ )

	{

		for( INT j = 0; j < 6; j++ )

		{

			// Compute the axis we are going to test.

			XMVECTOR Axis = XMVector3Cross( R.r[i], FrustumEdgeAxis[j] );

			// Find the min/max values of the projection of the frustum onto the axis.

			XMVECTOR FrustumMin, FrustumMax;

			FrustumMin = FrustumMax = XMVector3Dot( Axis, Corners[0] );

			for( INT k = 1; k < 8; k++ )

			{

				XMVECTOR Temp = XMVector3Dot( Axis, Corners[k] );

				FrustumMin = XMVectorMin( FrustumMin, Temp );

				FrustumMax = XMVectorMax( FrustumMax, Temp );

			}

			// Project the center of the box onto the axis.

			XMVECTOR Dist = XMVector3Dot( Center, Axis );

			// Project the axes of the box onto the axis to find the "radius" of the box.

			XMVECTOR Radius = XMVector3Dot( Axis, R.r[0] );

			Radius = XMVectorSelect( Radius, XMVector3Dot( Axis, R.r[1] ), SelectY );

			Radius = XMVectorSelect( Radius, XMVector3Dot( Axis, R.r[2] ), SelectZ );

			Radius = XMVector3Dot( Extents, XMVectorAbs( Radius ) );

			// if (center > max + radius || center < min - radius) reject;

			Outside = XMVectorOrInt( Outside, XMVectorGreater( Dist, FrustumMax + Radius ) );

			Outside = XMVectorOrInt( Outside, XMVectorLess( Dist, FrustumMin - Radius ) );

		}

	}

	if ( XMVector4EqualInt( Outside, XMVectorTrueInt() ) )

		return 0;

	// If we did not find a separating plane then the box must intersect the frustum.

	return 1;

}

//-----------------------------------------------------------------------------

// Exact axis alinged box vs frustum test.	Constructs an oriented box and uses

// the oriented box vs frustum test.

//

// Return values: 0 = no intersection,

//				  1 = intersection,

//				  2 = box is completely inside frustum

//-----------------------------------------------------------------------------

static INT IntersectAxisAlignedBoxFrustum( const AxisAlignedBox* pVolumeA, const Frustum* pVolumeB )

{

	XMASSERT( pVolumeA );

	XMASSERT( pVolumeB );

	// Make the axis aligned box oriented and do an OBB vs frustum test.

	OrientedBox BoxA;

	BoxA.Center = pVolumeA->Center;

	BoxA.Extents = pVolumeA->Extents;

	BoxA.Orientation.x = 0.0f;

	BoxA.Orientation.y = 0.0f;

	BoxA.Orientation.z = 0.0f;

	BoxA.Orientation.w = 1.0f;

	return IntersectOrientedBoxFrustum( &BoxA, pVolumeB );

}


//=============================================================================
// Patch_Chapter16 2016-10-22 by.progsecu
// 16장에 필요한 내용 추가 원본은 Common 폴더의 xnacollision.cpp 참조
static const XMVECTOR g_UnitVectorEpsilon =
{
	1.0e-4f, 1.0e-4f, 1.0e-4f, 1.0e-4f
};

//-----------------------------------------------------------------------------
// Return TRUE if the vector is a unit vector (length == 1).
//-----------------------------------------------------------------------------
static inline BOOL XMVector3IsUnit(FXMVECTOR V)
{
	XMVECTOR Difference = XMVector3Length(V) - XMVectorSplatOne();

	return XMVector4Less(XMVectorAbs(Difference), g_UnitVectorEpsilon);
}

//-----------------------------------------------------------------------------
// Compute the intersection of a ray (Origin, Direction) with a triangle 
// (V0, V1, V2).  Return TRUE if there is an intersection and also set *pDist 
// to the distance along the ray to the intersection.
// 
// The algorithm is based on Moller, Tomas and Trumbore, "Fast, Minimum Storage 
// Ray-Triangle Intersection", Journal of Graphics Tools, vol. 2, no. 1, 
// pp 21-28, 1997.
//-----------------------------------------------------------------------------
static BOOL IntersectRayTriangle(FXMVECTOR Origin, FXMVECTOR Direction, FXMVECTOR V0, CXMVECTOR V1, CXMVECTOR V2,
	FLOAT* pDist)
{
	XMASSERT(pDist);
	XMASSERT(XMVector3IsUnit(Direction));

	static const XMVECTOR Epsilon =
	{
		1e-20f, 1e-20f, 1e-20f, 1e-20f
	};

	XMVECTOR Zero = XMVectorZero();

	XMVECTOR e1 = V1 - V0;
	XMVECTOR e2 = V2 - V0;

	// p = Direction ^ e2;
	XMVECTOR p = XMVector3Cross(Direction, e2);

	// det = e1 * p;
	XMVECTOR det = XMVector3Dot(e1, p);

	XMVECTOR u, v, t;

	if (XMVector3GreaterOrEqual(det, Epsilon))
	{
		// Determinate is positive (front side of the triangle).
		XMVECTOR s = Origin - V0;

		// u = s * p;
		u = XMVector3Dot(s, p);

		XMVECTOR NoIntersection = XMVectorLess(u, Zero);
		NoIntersection = XMVectorOrInt(NoIntersection, XMVectorGreater(u, det));

		// q = s ^ e1;
		XMVECTOR q = XMVector3Cross(s, e1);

		// v = Direction * q;
		v = XMVector3Dot(Direction, q);

		NoIntersection = XMVectorOrInt(NoIntersection, XMVectorLess(v, Zero));
		NoIntersection = XMVectorOrInt(NoIntersection, XMVectorGreater(u + v, det));

		// t = e2 * q;
		t = XMVector3Dot(e2, q);

		NoIntersection = XMVectorOrInt(NoIntersection, XMVectorLess(t, Zero));

		if (XMVector4EqualInt(NoIntersection, XMVectorTrueInt()))
			return FALSE;
	}
	else if (XMVector3LessOrEqual(det, -Epsilon))
	{
		// Determinate is negative (back side of the triangle).
		XMVECTOR s = Origin - V0;

		// u = s * p;
		u = XMVector3Dot(s, p);

		XMVECTOR NoIntersection = XMVectorGreater(u, Zero);
		NoIntersection = XMVectorOrInt(NoIntersection, XMVectorLess(u, det));

		// q = s ^ e1;
		XMVECTOR q = XMVector3Cross(s, e1);

		// v = Direction * q;
		v = XMVector3Dot(Direction, q);

		NoIntersection = XMVectorOrInt(NoIntersection, XMVectorGreater(v, Zero));
		NoIntersection = XMVectorOrInt(NoIntersection, XMVectorLess(u + v, det));

		// t = e2 * q;
		t = XMVector3Dot(e2, q);

		NoIntersection = XMVectorOrInt(NoIntersection, XMVectorGreater(t, Zero));

		if (XMVector4EqualInt(NoIntersection, XMVectorTrueInt()))
			return FALSE;
	}
	else
	{
		// Parallel ray.
		return FALSE;
	}

	XMVECTOR inv_det = XMVectorReciprocal(det);

	t *= inv_det;

	// u * inv_det and v * inv_det are the barycentric cooridinates of the intersection.

	// Store the x-component to *pDist
	XMStoreFloat(pDist, t);

	return TRUE;
}

//-----------------------------------------------------------------------------
// Compute the intersection of a ray (Origin, Direction) with an axis aligned 
// box using the slabs method.
//-----------------------------------------------------------------------------
static BOOL IntersectRayAxisAlignedBox(FXMVECTOR Origin, FXMVECTOR Direction, const AxisAlignedBox* pVolume, FLOAT* pDist)
{
	XMASSERT(pVolume);
	XMASSERT(pDist);
	XMASSERT(XMVector3IsUnit(Direction));

	static const XMVECTOR Epsilon =
	{
		1e-20f, 1e-20f, 1e-20f, 1e-20f
	};
	static const XMVECTOR FltMin =
	{
		-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX
	};
	static const XMVECTOR FltMax =
	{
		FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX
	};

	// Load the box.
	XMVECTOR Center = XMLoadFloat3(&pVolume->Center);
	XMVECTOR Extents = XMLoadFloat3(&pVolume->Extents);

	// Adjust ray origin to be relative to center of the box.
	XMVECTOR TOrigin = Center - Origin;

	// Compute the dot product againt each axis of the box.
	// Since the axii are (1,0,0), (0,1,0), (0,0,1) no computation is necessary.
	XMVECTOR AxisDotOrigin = TOrigin;
	XMVECTOR AxisDotDirection = Direction;

	// if (fabs(AxisDotDirection) <= Epsilon) the ray is nearly parallel to the slab.
	XMVECTOR IsParallel = XMVectorLessOrEqual(XMVectorAbs(AxisDotDirection), Epsilon);

	// Test against all three axii simultaneously.
	XMVECTOR InverseAxisDotDirection = XMVectorReciprocal(AxisDotDirection);
	XMVECTOR t1 = (AxisDotOrigin - Extents) * InverseAxisDotDirection;
	XMVECTOR t2 = (AxisDotOrigin + Extents) * InverseAxisDotDirection;

	// Compute the max of min(t1,t2) and the min of max(t1,t2) ensuring we don't
	// use the results from any directions parallel to the slab.
	XMVECTOR t_min = XMVectorSelect(XMVectorMin(t1, t2), FltMin, IsParallel);
	XMVECTOR t_max = XMVectorSelect(XMVectorMax(t1, t2), FltMax, IsParallel);

	// t_min.x = maximum( t_min.x, t_min.y, t_min.z );
	// t_max.x = minimum( t_max.x, t_max.y, t_max.z );
	t_min = XMVectorMax(t_min, XMVectorSplatY(t_min));  // x = max(x,y)
	t_min = XMVectorMax(t_min, XMVectorSplatZ(t_min));  // x = max(max(x,y),z)
	t_max = XMVectorMin(t_max, XMVectorSplatY(t_max));  // x = min(x,y)
	t_max = XMVectorMin(t_max, XMVectorSplatZ(t_max));  // x = min(min(x,y),z)

														// if ( t_min > t_max ) return FALSE;
	XMVECTOR NoIntersection = XMVectorGreater(XMVectorSplatX(t_min), XMVectorSplatX(t_max));

	// if ( t_max < 0.0f ) return FALSE;
	NoIntersection = XMVectorOrInt(NoIntersection, XMVectorLess(XMVectorSplatX(t_max), XMVectorZero()));

	// if (IsParallel && (-Extents > AxisDotOrigin || Extents < AxisDotOrigin)) return FALSE;
	XMVECTOR ParallelOverlap = XMVectorInBounds(AxisDotOrigin, Extents);
	NoIntersection = XMVectorOrInt(NoIntersection, XMVectorAndCInt(IsParallel, ParallelOverlap));

	if (!XMVector3AnyTrue(NoIntersection))
	{
		// Store the x-component to *pDist
		XMStoreFloat(pDist, t_min);
		return TRUE;
	}

	return FALSE;
}
//=============================================================================

//-----------------------------------------------------------------------------
#define XMVectorPermute( a, b, c) XMVectorPermute(a,b, c.i[0], c.i[1], c.i[2], c.i[3] )
static BOOL IntersectTriangleAxisAlignedBox( FXMVECTOR V0, FXMVECTOR V1, FXMVECTOR V2, const AxisAlignedBox* pVolume )
{
    XMASSERT( pVolume );

    static CONST XMVECTORI32 Permute0W1Z0Y0X =
                 {
                    XM_PERMUTE_0W, XM_PERMUTE_1Z, XM_PERMUTE_0Y, XM_PERMUTE_0X
                 };
    static CONST XMVECTORI32 Permute0Z0W1X0Y =
                 {
                    XM_PERMUTE_0Z, XM_PERMUTE_0W, XM_PERMUTE_1X, XM_PERMUTE_0Y
                 };
    static CONST XMVECTORI32 Permute1Y0X0W0Z =
                 {
                    XM_PERMUTE_1Y, XM_PERMUTE_0X, XM_PERMUTE_0W, XM_PERMUTE_0Z
                 };

    XMVECTOR Zero = XMVectorZero();

    // Load the box.
    XMVECTOR Center = XMLoadFloat3( &pVolume->Center );
    XMVECTOR Extents = XMLoadFloat3( &pVolume->Extents );

    XMVECTOR BoxMin = Center - Extents;
    XMVECTOR BoxMax = Center + Extents;

    // Test the axes of the box (in effect test the AAB against the minimal AAB 
    // around the triangle).
    XMVECTOR TriMin = XMVectorMin( XMVectorMin( V0, V1 ), V2 );
    XMVECTOR TriMax = XMVectorMax( XMVectorMax( V0, V1 ), V2 );

    // for each i in (x, y, z) if a_min(i) > b_max(i) or b_min(i) > a_max(i) then disjoint
    XMVECTOR Disjoint = XMVectorOrInt( XMVectorGreater( TriMin, BoxMax ), XMVectorGreater( BoxMin, TriMax ) );
    if( XMVector3AnyTrue( Disjoint ) )
        return FALSE;

    // Test the plane of the triangle.
    XMVECTOR Normal = XMVector3Cross( V1 - V0, V2 - V0 );
    XMVECTOR Dist = XMVector3Dot( Normal, V0 );

    // Assert that the triangle is not degenerate.
    XMASSERT( !XMVector3Equal( Normal, Zero ) );

    // for each i in (x, y, z) if n(i) >= 0 then v_min(i)=b_min(i), v_max(i)=b_max(i)
    // else v_min(i)=b_max(i), v_max(i)=b_min(i)
    XMVECTOR NormalSelect = XMVectorGreater( Normal, Zero );
    XMVECTOR V_Min = XMVectorSelect( BoxMax, BoxMin, NormalSelect );
    XMVECTOR V_Max = XMVectorSelect( BoxMin, BoxMax, NormalSelect );

    // if n dot v_min + d > 0 || n dot v_max + d < 0 then disjoint
    XMVECTOR MinDist = XMVector3Dot( V_Min, Normal );
    XMVECTOR MaxDist = XMVector3Dot( V_Max, Normal );

    XMVECTOR NoIntersection = XMVectorGreater( MinDist, Dist );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( MaxDist, Dist ) );

    // Move the box center to zero to simplify the following tests.
    XMVECTOR TV0 = V0 - Center;
    XMVECTOR TV1 = V1 - Center;
    XMVECTOR TV2 = V2 - Center;

    // Test the edge/edge axes (3*3).
    XMVECTOR e0 = TV1 - TV0;
    XMVECTOR e1 = TV2 - TV1;
    XMVECTOR e2 = TV0 - TV2;

    // Make w zero.
    e0 = XMVectorInsert( e0, Zero, 0, 0, 0, 0, 1 );
    e1 = XMVectorInsert( e1, Zero, 0, 0, 0, 0, 1 );
    e2 = XMVectorInsert( e2, Zero, 0, 0, 0, 0, 1 );

    XMVECTOR Axis;
    XMVECTOR p0, p1, p2;
    XMVECTOR Min, Max;
    XMVECTOR Radius;

    // Axis == (1,0,0) x e0 = (0, -e0.z, e0.y)
    Axis = XMVectorPermute( e0, -e0, Permute0W1Z0Y0X );
    p0 = XMVector3Dot( TV0, Axis );
    // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
    p2 = XMVector3Dot( TV2, Axis );
    Min = XMVectorMin( p0, p2 );
    Max = XMVectorMax( p0, p2 );
    Radius = XMVector3Dot( Extents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (1,0,0) x e1 = (0, -e1.z, e1.y)
    Axis = XMVectorPermute( e1, -e1, Permute0W1Z0Y0X );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( Extents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (1,0,0) x e2 = (0, -e2.z, e2.y)
    Axis = XMVectorPermute( e2, -e2, Permute0W1Z0Y0X );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( Extents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,1,0) x e0 = (e0.z, 0, -e0.x)
    Axis = XMVectorPermute( e0, -e0, Permute0Z0W1X0Y );
    p0 = XMVector3Dot( TV0, Axis );
    // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
    p2 = XMVector3Dot( TV2, Axis );
    Min = XMVectorMin( p0, p2 );
    Max = XMVectorMax( p0, p2 );
    Radius = XMVector3Dot( Extents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,1,0) x e1 = (e1.z, 0, -e1.x)
    Axis = XMVectorPermute( e1, -e1, Permute0Z0W1X0Y );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( Extents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,0,1) x e2 = (e2.z, 0, -e2.x)
    Axis = XMVectorPermute( e2, -e2, Permute0Z0W1X0Y );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( Extents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,0,1) x e0 = (-e0.y, e0.x, 0)
    Axis = XMVectorPermute( e0, -e0, Permute1Y0X0W0Z );
    p0 = XMVector3Dot( TV0, Axis );
    // p1 = XMVector3Dot( V1, Axis ); // p1 = p0;
    p2 = XMVector3Dot( TV2, Axis );
    Min = XMVectorMin( p0, p2 );
    Max = XMVectorMax( p0, p2 );
    Radius = XMVector3Dot( Extents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,0,1) x e1 = (-e1.y, e1.x, 0)
    Axis = XMVectorPermute( e1, -e1, Permute1Y0X0W0Z );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p1;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( Extents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    // Axis == (0,0,1) x e2 = (-e2.y, e2.x, 0)
    Axis = XMVectorPermute( e2, -e2, Permute1Y0X0W0Z );
    p0 = XMVector3Dot( TV0, Axis );
    p1 = XMVector3Dot( TV1, Axis );
    // p2 = XMVector3Dot( V2, Axis ); // p2 = p0;
    Min = XMVectorMin( p0, p1 );
    Max = XMVectorMax( p0, p1 );
    Radius = XMVector3Dot( Extents, XMVectorAbs( Axis ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorGreater( Min, Radius ) );
    NoIntersection = XMVectorOrInt( NoIntersection, XMVectorLess( Max, -Radius ) );

    return XMVector4NotEqualInt( NoIntersection, XMVectorTrueInt() );
}

} // namespace
#endif