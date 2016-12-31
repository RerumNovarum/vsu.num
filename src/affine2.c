#include <vsu/num.h>
#include <stdarg.h>

void
affine2hom(struct affine2 T, RR *X, RR *Y, RR *ALPHA)
{
    RR x = *X, y = *Y, alpha = *ALPHA;
    *X = T.a11*x + T.a12*y + T.b1;
    *Y = T.a21*x + T.a22*y + T.b2;
    *ALPHA = alpha * T.alpha;
}

struct vec2rr
affine2rr_immut(struct affine2 T, struct vec2rr v)
{
    RR x = v.x, y = v.y;
    v.x = (T.a11*x + T.a12*y + T.b1)/T.alpha;
    v.y = (T.a21*x + T.a22*y + T.b2)/T.alpha;
    return v;
}

void
affine2rr(struct affine2 T, RR *X, RR *Y)
{
    RR x = *X, y = *Y;
    *X = (T.a11*x + T.a12*y + T.b1)/T.alpha;
    *Y = (T.a21*x + T.a22*y + T.b2)/T.alpha;
}

struct affine2
affine2mul(struct affine2 A, struct affine2 B)
{
    struct affine2 new;
    new.a11 = A.a11*B.a11 + A.a12*B.a21;
    new.a12 = A.a11*B.a12 + A.a12*B.a22;
    new.a21 = A.a21*B.a11 + A.a22*B.a21;;
    new.a22 = A.a21*B.a12 + A.a22*B.a22;
    new.b1 = A.a11*B.b1 + A.a12*B.b2 + A.b1*B.alpha;
    new.b2 = A.a21*B.b1 + A.a22*B.b2 + A.b2*B.alpha;
    new.alpha = A.alpha*B.alpha;
    return new;
}

struct affine2
affine2mul_n(size_t n, struct affine2 transforms, ...)
{
    va_list ap;
    va_start(ap, transforms);

    struct affine2 T = AFFINE2_ID;
    for (; n > 0; --n, transforms=va_arg(ap, struct affine2))
        /* TODO: it is probably better to pass it via pointers */
        T = affine2mul(T, transforms);

    va_end(ap);
    return T;
}

struct affine2
affine2rot(RR phi)
{
    RR c = cosl(phi);
    RR s = sinl(phi);
    return (struct affine2) { c, -s, s, c, 0, 0, 1 };
}
struct affine2
affine2tr(RR x, RR y)
{
    return (struct affine2) { 1, 0, 0, 1, x, y, 1 };
}
struct affine2
affine2scale(RR x, RR y)
{
    return (struct affine2) { x, 0, 0, y, 0, 0, 1 };
}
