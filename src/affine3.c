#include <vsu/num.h>
#include <stdarg.h>

struct affine3
affine3mul(struct affine3 A, struct affine3 B)
{
    struct affine3 new;
    new.a11 = A.a11*B.a11 + A.a12*B.a21 + A.a13*B.a31;
    new.a12 = A.a11*B.a12 + A.a12*B.a22 + A.a13*B.a32;
    new.a13 = A.a11*B.a13 + A.a12*B.a23 + A.a13*B.a33;
    new.a21 = A.a21*B.a11 + A.a22*B.a21 + A.a23*B.a31;
    new.a22 = A.a21*B.a12 + A.a22*B.a22 + A.a23*B.a32;
    new.a23 = A.a21*B.a13 + A.a22*B.a23 + A.a23*B.a33;
    new.a31 = A.a31*B.a11 + A.a32*B.a21 + A.a33*B.a31;
    new.a32 = A.a31*B.a12 + A.a32*B.a22 + A.a33*B.a32;
    new.a33 = A.a31*B.a13 + A.a32*B.a23 + A.a33*B.a33;
    new.b1 = A.a11*B.b1 + A.a12*B.b2 + A.a13*B.b3 + A.b1*B.alpha;
    new.b2 = A.a21*B.b1 + A.a22*B.b2 + A.a23*B.b3 + A.b2*B.alpha;
    new.b3 = A.a31*B.b1 + A.a32*B.b2 + A.a33*B.b3 + A.b3*B.alpha;
    new.alpha = A.alpha*B.alpha;
    return new;
}

struct affine3
affine3mul_n(size_t n, struct affine3 transforms, ...)
{
    va_list ap;
    va_start(ap, transforms);

    struct affine3 T = AFFINE3_ID;
    for (; n > 0; --n, transforms=va_arg(ap, struct affine3))
        /* TODO: it is probably better to pass it via pointers */
        T = affine3mul(T, transforms);

    va_end(ap);
    return T;
}

struct affine3
affine3rotx(RR phi)
{
    RR c = cosl(phi);
    RR s = sinl(phi);
    struct affine3 T = {
        .a11=1, .a12=0, .a13=0,
        .a21=0, .a22=c, .a23=-s,
        .a31=0, .a32=s, .a33=c,
        .b1=0, .b2=0, .b3=0, .alpha=1
    };
    return T;
}

struct affine3
affine3roty(RR phi)
{
    RR c = cosl(phi);
    RR s = sinl(phi);
    /* signs are inverted, since rotation around oY
     * corresponds to a rotation around the origin
     * in the plane ZX
     */
    struct affine3 T = {
        .a11=c,  .a12=0, .a13=s,
        .a21=0,  .a22=1, .a23=0,
        .a31=-s, .a32=0, .a33=c,
        .b1=0, .b2=0, .b3=0, .alpha=1
    };
    return T;
}

struct affine3
affine3rotz(RR phi)
{
    RR c = cosl(phi);
    RR s = sinl(phi);
    /* signs are inverted, since rotation around oY
     * corresponds to a rotation around the origin
     * in the plane ZX
     */
    struct affine3 T = {
        .a11=c, .a12=-s, .a13=0, .b1=0,
        .a21=s, .a22=c,  .a23=0, .b2=0,
        .a31=0, .a32=0,  .a33=1, .b3=0,
                                 .alpha=1
    };
    return T;
}

struct affine3
affine3scale(RR x, RR y, RR z)
{
    struct affine3 T = {
        .a11=x, .a12=0, .a13=0, .b1=0,
        .a21=0, .a22=y, .a23=0, .b2=0,
        .a31=0, .a32=0, .a33=z, .b3=0,
                                .alpha=1
    };
    return T;
}

struct affine3
affine3tr(RR x, RR y, RR z)
{
    struct affine3 T = {
        .a11=1, .a12=0, .a13=0, .b1=x,
        .a21=0, .a22=1, .a23=0, .b2=y,
        .a31=0, .a32=0, .a33=1, .b3=z,
                                .alpha=1
    };
    return T;
}

void
affine3apply_rr(struct affine3 T, RR *X, RR *Y, RR *Z)
{
    RR x = *X, y = *Y, z = *Z;
    *X = (T.a11*x + T.a12*y + T.a13*z + T.b1)/T.alpha;
    *Y = (T.a21*x + T.a22*y + T.a23*z + T.b2)/T.alpha;
    *Z = (T.a31*x + T.a32*y + T.a33*z + T.b2)/T.alpha;
}
