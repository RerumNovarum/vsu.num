#include <string.h>
#include <err.h>
#include <vsu/num.h>

int main()
{
    char *sn = "12345";
    NN n = 12345;
    NN cn = num_sgetn(sn, strlen(sn));
    if (cn == n)
        fputs("num_sgetn seems ok\n", stderr);
    else
        err(1, "num_sgetn failed on %s; computed=%zu", sn, cn);

    char *sz = "-321";
    ZZ z = -321;
    ZZ cz = num_sgetz(sz, strlen(sz));
    if (cz == z)
        fputs("num_sgetz seems ok\n", stderr);
    else
        err(1, "num_sgetz failed on %s; computed=%Ld", sz, cz);

    char *sr = "-1.3125";
    RR r = -1.3125;
    RR cr = num_sgetr(sr, strlen(sr));
    if (cr == r)
        fputs("num_sgetr seems ok\n", stderr);
    else
        err(1, "num_sgetr failed on %s; computed=%Lf", sr, cr);
    char *sc = "1+2i";
    CC c = 1+2*I;
    CC cc = num_sgetc(sc, strlen(sc));
    if (cc == c)
        fputs("num_sgetc seems ok\n", stderr);
    else
        err(1, "num_sgetc failed on %s; computed=%Lf+%Lfi", sc, REAL(cc), IMAG(cc));
    return 0;
}
