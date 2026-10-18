#include "fft_double.h"

#ifdef FFTW2
#ifdef FFT_THREADED

static int g_thread_count;

void
set_fft_double_thread_count(int nthreads)
{
    g_thread_count = nthreads;
}


int
get_fft_double_thread_count(void)
{
    return g_thread_count;
}

#endif
#endif
