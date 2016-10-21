// Patch_Chapter15.h

// fix InstancingAndCullingDemo.cpp to fit DirectXMath and DirectXCollision.h

// the demo is from the book: 3D Game Programming with DirectX 11, Chapter 15

// the functions is from xnacollision.cpp

// endrollex

// Changed to support Unicode. 2016-10-21 by.progsecu

#include "DirectXCollision.h"
#include <Windows.h>
#include <cfloat>
//

#if defined(DEBUG) | defined(_DEBUG)

	#define XMASSERT(Expression) ((VOID)((Expression) || (XMAssert(#Expression, __FILEW__, __LINE__), 0)))

#else

	#define XMASSERT(Expression) ((VOID)0)

#endif

VOID XMAssert(_In_z_ CONST CHAR* pExpression, _In_z_ CONST WCHAR* pFileName, UINT LineNumber)

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

VOID ComputeFrustumFromProjection( BoundingFrustum* pOut, XMMATRIX* pProjection )

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

VOID TransformFrustum( Frustum* pOut, const Frustum* pIn, FLOAT Scale, FXMVECTOR Rotation, FXMVECTOR Translation )

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

INT IntersectOrientedBoxFrustum( const OrientedBox* pVolumeA, const Frustum* pVolumeB )

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

INT IntersectAxisAlignedBoxFrustum( const AxisAlignedBox* pVolumeA, const Frustum* pVolumeB )

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

} // namespace

