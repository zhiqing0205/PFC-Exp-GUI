// pfc_common.h - Shared infrastructure for all PFC simulation models.
// Provides MPI stubs, deterministic RNG, argument parsing helpers,
// grid size macros, and common utility functions.
#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <math.h>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <system_error>

#if defined(PFC_USE_MPI)
#include <fftw3-mpi.h>
#else
#include <fftw3.h>

// Minimal MPI stubs for single-process builds.
using MPI_Comm = int;
using MPI_Datatype = int;
using MPI_Op = int;

struct MPI_Status {
    int dummy = 0;
};

static constexpr MPI_Comm MPI_COMM_WORLD = 0;
static constexpr int MPI_PROC_NULL = -1;

static constexpr MPI_Datatype MPI_INT = 0;
static constexpr MPI_Datatype MPI_FLOAT = 1;
static constexpr MPI_Datatype MPI_DOUBLE = 2;

static constexpr MPI_Op MPI_MAX = 0;

static inline int MPI_Init(int*, char***) { return 0; }
static inline int MPI_Finalize() { return 0; }
static inline int MPI_Comm_rank(MPI_Comm, int* rank) {
    if (rank) *rank = 0;
    return 0;
}
static inline int MPI_Comm_size(MPI_Comm, int* size) {
    if (size) *size = 1;
    return 0;
}
static inline int MPI_Bcast(void*, int, MPI_Datatype, int, MPI_Comm) { return 0; }
static inline int MPI_Barrier(MPI_Comm) { return 0; }

static inline size_t MPI_Type_size_bytes(MPI_Datatype datatype) {
    switch (datatype) {
        case MPI_INT:
            return sizeof(int);
        case MPI_FLOAT:
            return sizeof(float);
        case MPI_DOUBLE:
            return sizeof(double);
        default:
            return 1;
    }
}

static inline int MPI_Allreduce(const void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op,
                                MPI_Comm) {
    if (!sendbuf || !recvbuf || count <= 0) return 0;
    const size_t bytes = static_cast<size_t>(count) * MPI_Type_size_bytes(datatype);
    std::memcpy(recvbuf, sendbuf, bytes);
    return 0;
}

static inline int MPI_Gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int,
                             MPI_Datatype, int root, MPI_Comm) {
    if (root != 0) return 0;
    if (!sendbuf || !recvbuf || sendcount <= 0) return 0;
    const size_t bytes = static_cast<size_t>(sendcount) * MPI_Type_size_bytes(sendtype);
    std::memcpy(recvbuf, sendbuf, bytes);
    return 0;
}

static inline int MPI_Gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf,
                              const int* recvcounts, const int* displs, MPI_Datatype recvtype, int root,
                              MPI_Comm) {
    if (root != 0) return 0;
    if (!sendbuf || !recvbuf || sendcount <= 0) return 0;
    (void)recvcounts;
    (void)recvtype;
    const size_t type_size = MPI_Type_size_bytes(sendtype);
    const int displacement = displs ? displs[0] : 0;
    std::memcpy(static_cast<char*>(recvbuf) + static_cast<size_t>(displacement) * type_size, sendbuf,
                static_cast<size_t>(sendcount) * type_size);
    return 0;
}

static inline int MPI_Sendrecv(const void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int, void* recvbuf,
                               int recvcount, MPI_Datatype recvtype, int source, int, MPI_Comm, MPI_Status*) {
    if (dest == MPI_PROC_NULL || source == MPI_PROC_NULL) return 0;
    if (dest != 0 || source != 0) return 0;
    if (!sendbuf || !recvbuf || sendcount <= 0 || recvcount <= 0) return 0;
    const size_t send_bytes = static_cast<size_t>(sendcount) * MPI_Type_size_bytes(sendtype);
    const size_t recv_bytes = static_cast<size_t>(recvcount) * MPI_Type_size_bytes(recvtype);
    std::memcpy(recvbuf, sendbuf, std::min(send_bytes, recv_bytes));
    return 0;
}

// FFTW MPI API shims (serial execution).
static inline int fftw_mpi_init() { return 0; }

static inline ptrdiff_t fftw_mpi_local_size_3d(ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2, MPI_Comm, ptrdiff_t* local_n0,
                                               ptrdiff_t* local_0_start) {
    if (local_n0) *local_n0 = n0;
    if (local_0_start) *local_0_start = 0;
    return n0 * n1 * n2;
}

static inline fftw_plan fftw_mpi_plan_dft_3d(ptrdiff_t n0, ptrdiff_t n1, ptrdiff_t n2, fftw_complex* in, fftw_complex* out,
                                             MPI_Comm, int sign, unsigned flags) {
    return fftw_plan_dft_3d(static_cast<int>(n0), static_cast<int>(n1), static_cast<int>(n2), in, out, sign, flags);
}
#endif

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <ctime>
#include <chrono>

// ---------------------------------------------------------------------------
// Grid size macros (compile-time configurable via CMake)
// ---------------------------------------------------------------------------
#ifndef PFC_GRID_L
#define PFC_GRID_L 256
#endif
#ifndef PFC_GRID_M
#define PFC_GRID_M 256
#endif
#ifndef PFC_GRID_N
#define PFC_GRID_N 1
#endif

// ---------------------------------------------------------------------------
// Common constants
// ---------------------------------------------------------------------------
#define TTT 100
#define pi 3.1415926
#define RD 1
#define FIXED_SEED 20200604u
#define Radius 140
#define Num_euc 4
#define rad pi/180.

// ---------------------------------------------------------------------------
// Deterministic RNG helpers
// ---------------------------------------------------------------------------
// NOTE: std::rand() is implementation-defined and yields different sequences on
// different platforms. We use a small deterministic generator so that the same
// --seed produces comparable results on Windows/macOS/Linux.
static inline uint64_t splitmix64_next(uint64_t& state) {
    uint64_t z = (state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static inline double rng_u01(uint64_t& state) {
    // [0, 1) with 53 bits of precision.
    constexpr double inv = 1.0 / 9007199254740992.0;  // 2^53
    return static_cast<double>(splitmix64_next(state) >> 11) * inv;
}

static inline double rng_u11(uint64_t& state) {
    // (-1, 1)
    return 2.0 * rng_u01(state) - 1.0;
}

// ---------------------------------------------------------------------------
// Argument parsing helpers
// ---------------------------------------------------------------------------
static inline bool PFC_ParseIntArg(const char* value, int& out) {
    if (!value) return false;
    char* end = nullptr;
    errno = 0;
    long v = strtol(value, &end, 10);
    if (errno != 0 || end == value || *end != '\0') return false;
    if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max()) return false;
    out = static_cast<int>(v);
    return true;
}

static inline bool PFC_ParseFloatArg(const char* value, float& out) {
    if (!value) return false;
    char* end = nullptr;
    errno = 0;
    float v = strtof(value, &end);
    if (errno != 0 || end == value || *end != '\0') return false;
    out = v;
    return true;
}

static inline bool PFC_ParseUintArg(const char* value, unsigned int& out) {
    if (!value) return false;
    char* end = nullptr;
    errno = 0;
    unsigned long long v = strtoull(value, &end, 10);
    if (errno != 0 || end == value || *end != '\0') return false;
    if (v > std::numeric_limits<unsigned int>::max()) return false;
    out = static_cast<unsigned int>(v);
    return true;
}

static inline bool PFC_HasPrefix(std::string_view s, std::string_view prefix) {
    return s.size() >= prefix.size() && s.substr(0, prefix.size()) == prefix;
}

// Helper to require a value for a CLI option.
// Usage: const char* v = PFC_RequireValue(argc, argv, i, "--option", myid);
static inline const char* PFC_RequireValue(int argc, char** argv, int& i, const char* opt, int myid) {
    if (i + 1 >= argc) {
        if (myid == 0) {
            std::cerr << "Missing value after " << opt << "\n";
        }
        return nullptr;
    }
    return argv[++i];
}

// ---------------------------------------------------------------------------
// Common utility functions (present identically in all three models)
// ---------------------------------------------------------------------------
static inline float pfc_ks(int i, int j, int k, float* dkrefl, float* dkrefm, float* dkrefn) {
    return dkrefl[i]*dkrefl[i] + dkrefm[j]*dkrefm[j] + dkrefn[k]*dkrefn[k];
}

static inline float pfc_min3(float a, float b, float c) {
    float aa = fabs(a), bb = fabs(b), cc = fabs(c);
    return ((aa < bb ? aa : bb) < cc) ? (aa < bb ? aa : bb) : cc;
}

// ---------------------------------------------------------------------------
// Output directory setup
// ---------------------------------------------------------------------------
static inline bool PFC_SetupOutputDir(const std::string& output_dir, int myid) {
    if (output_dir.empty()) return true;
    int ok = 1;
    if (myid == 0) {
        std::error_code ec;
        std::filesystem::create_directories(output_dir, ec);
        if (ec) {
            std::cerr << "Failed to create --outdir '" << output_dir << "': " << ec.message() << std::endl;
            ok = 0;
        }
    }
    MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (!ok) return false;
    std::error_code ec;
    std::filesystem::current_path(output_dir, ec);
    if (ec) {
        std::cerr << "Failed to chdir to --outdir '" << output_dir << "': " << ec.message() << std::endl;
        return false;
    }
    return true;
}

// ---------------------------------------------------------------------------
// Spherical harmonics and associated Legendre (used by misfit and elastic)
// ---------------------------------------------------------------------------
static inline unsigned long pfc_facto(int h) {
    int g = 1;
    for (int i = 1; i <= h; i++) { g = g * i; }
    return g;
}

static inline float pfc_lg_poly(int l, float x) {
    float plm2 = 1, plm1 = x, pl = 0;
    if (l == 0) return 1;
    if (l == 1) return x;
    else
        for (int i = 2; i <= l; i++) {
            pl = ((2 * i - 1) * x * float(plm1) - (i - 1) * plm2) / float(i);
            plm2 = plm1;
            plm1 = pl;
        }
    return pl;
}

static inline float pfc_plgndr(int l, int m, float x) {
    float fact, pll, pmm, pmmp1, somx2;
    if (m < 0) m = -m;
    if (m > l || fabs(x) > 1.0) std::cout << "Bad arguments in routine plgndr" << std::endl;
    if (m == 0) return pfc_lg_poly(l, x);
    pmm = 1.0;
    if (m > 0) {
        somx2 = sqrt((1.0 - x) * (1.0 + x));
        fact = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm *= -fact * somx2;
            fact += 2.0;
        }
    }
    if (l == m) return pmm;
    else {
        pmmp1 = x * (2 * m + 1) * pmm;
        if (l == (m + 1)) return pmmp1;
        else {
            for (int ll = m + 2; ll <= l; ll++) {
                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                pmm = pmmp1;
                pmmp1 = pll;
            }
            return pll;
        }
    }
}

static inline float pfc_Ylmtf_r(int l, int m, float theta, float fai) {
    if (m < 0) m = -m;
    float fnorm = pow(-1.0, m) * sqrt((2 * l + 1) / 4. / pi) * sqrt(float(pfc_facto(l - m)) / float(pfc_facto(l + m)));
    float plm = pfc_plgndr(l, m, cos(theta));
    float comp_m = cos(m * fai);
    return fnorm * plm * comp_m;
}

static inline float pfc_Ylmtf_c(int l, int m, float theta, float fai) {
    float fnorm, plm, comp_m;
    if (m < 0) {
        m = -m;
        fnorm = pow(-1.0, m + m) * sqrt((2 * l + 1) / 4. / pi) * sqrt(float(pfc_facto(l - m)) / float(pfc_facto(l + m)));
        plm = pfc_plgndr(l, m, cos(theta));
        comp_m = -sin(m * fai);
    } else {
        fnorm = pow(-1.0, m) * sqrt((2 * l + 1) / 4. / pi) * sqrt(float(pfc_facto(l - m)) / float(pfc_facto(l + m)));
        plm = pfc_plgndr(l, m, cos(theta));
        comp_m = sin(m * fai);
    }
    return fnorm * plm * comp_m;
}

// ---------------------------------------------------------------------------
// Model entry point declarations
// ---------------------------------------------------------------------------
int run_misfit(int argc, char** argv);
int run_cvd(int argc, char** argv);
int run_elastic(int argc, char** argv);
