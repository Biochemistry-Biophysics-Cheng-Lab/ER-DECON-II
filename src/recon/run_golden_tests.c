/*
 * run_golden_tests.c -- driver for the delivered golden unit tests.
 * The Test*() functions live in JacobiNL{2,3}DduTests.c (delivered); no
 * main() was delivered with them.  JacobiNL4Ddu is out of scope (its
 * Decon4DParams struct was never delivered) and is not built.
 *
 * Exit status 0 = all pass.  Until Res3D, LineSearch3D2n and Derivative are
 * reconstructed every test reports "function call reported nonzero status".
 *
 * ERDecon reference build, 2026-09-18.
 */
#include "JacobiNL3DduTests.h"
#include "JacobiNL2DduTests.h"
#include <stdio.h>

int
main(void)
{
    int status = 0;
    int s;

    (void) printf("== JacobiNL3Ddu golden tests ==\n");
    s = TestJacobiNL3Ddu();
    (void) printf("TestJacobiNL3Ddu -> %d\n", s);
    status |= (s != 0);

    (void) printf("== JacobiNL2Ddu golden tests ==\n");
    s = TestJacobiNL2Ddu();
    (void) printf("TestJacobiNL2Ddu -> %d\n", s);
    status |= (s != 0);

    (void) printf(status ? "SOME TESTS FAILED\n" : "ALL TESTS PASSED\n");
    return status;
}
