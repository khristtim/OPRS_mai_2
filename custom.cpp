//---------------------------------------------------------------------------

#include "custom.h"
#include <math.h>

//---------------------------------------------------------------------------
// Задача Аренсторфа (начальные условия 1)

const long double TArenstorfModel::m  = 0.012277471L;

TArenstorfModel::TArenstorfModel( long double t0, long double t1, long double step )
    : TModel()
{
    t0_   = t0;
    t1_   = t1;
    step_ = step;

    X0_.resize(4);
    X0_[0] = 0.994L;
    X0_[1] = 0.0L;
    X0_[2] = 0.0L;
    X0_[3] = -2.00158510637908252240537862224L;
}

//---------------------------------------------------------------------------

TVector TArenstorfModel::getRight( const TVector& X, long double /*t*/ )
{
    TVector Y(4);

    const long double mu  = m;
    const long double mu1 = 1.0L - m;

    const long double y1 = X[0];
    const long double y2 = X[1];
    const long double v1 = X[2];
    const long double v2 = X[3];

    const long double d1 = powl( (y1 + mu)*(y1 + mu) + y2*y2, 1.5L );
    const long double d2 = powl( (y1 - mu1)*(y1 - mu1) + y2*y2, 1.5L );

    Y[0] = v1;
    Y[1] = v2;
    Y[2] = y1 + 2.0L*v2 - mu1*(y1 + mu)/d1 - mu*(y1 - mu1)/d2;
    Y[3] = y2 - 2.0L*v1 - mu1*y2/d1 - mu*y2/d2;

    return Y;
}

//---------------------------------------------------------------------------
// Задача Аренсторфа (начальные условия 2)

TArenstorfModel2::TArenstorfModel2( long double t0, long double t1, long double step )
    : TArenstorfModel( t0, t1, step )
{
    X0_[0] = 0.994L;
    X0_[1] = 0.0L;
    X0_[2] = 0.0L;
    X0_[3] = -2.0317326295573368357302057924L;
}
//---------------------------------------------------------------------------
