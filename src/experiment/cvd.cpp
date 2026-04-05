// CVD (Controlled Vapor Deposition) PFC simulation model
// Adapted from ref/cvd.cpp for unified CLI framework
// Date: 2014.09.25

#include "pfc_common.h"

static constexpr int L = PFC_GRID_L;
static constexpr int M = PFC_GRID_M;
static constexpr int N = PFC_GRID_N;

using namespace std;

// Use shared utilities from pfc_common.h
#define ks pfc_ks
#define min3 pfc_min3
#define Ylmtf_r pfc_Ylmtf_r
#define Ylmtf_c pfc_Ylmtf_c
#define plgndr pfc_plgndr
#define lg_poly pfc_lg_poly
#define facto pfc_facto

static int total_steps = 2000;
static std::string output_dir;
static unsigned int base_seed = FIXED_SEED;

static float a_100=1.0, a_110=1.414*a_100, a_111=1.732*a_100; 
static float ai=10;





static float a_lattice=1.0,b_lattice=sqrt(3.0)/2.0; //lattice A and B for square

static float u0 = 0.01, sig = 0.03, con0 = 0.2, dt = 0.05, dx = 0.125;
//static float u0=0.01,sig=0.03,con0=0.2.,dt=0.01,dx=0.125;
static float gaa1=0.8, gaa2=gaa1*sqrt(2.); //for A  //alpha parts in Gauss terms
static float gab1=0.8, gab2=gaa2*sqrt(2.); //for B
static float niu=1.4,ksi=1.,epsilona0=0.0,epsilonb0=0.0,epsilona=0.0,epsilonb=0.0;
static float  ka0=2.0*pi/a_lattice,ka1=2.0*pi,ka2=sqrt(2.0)*ka1,sigMa1=0.55,sigMa2=0.55; //for square A

static float kb0=2.5*pi/b_lattice,kb1=5.0*pi/sqrt(3.0),kb2=kb1*sqrt(2.0),sigMb1=0.55,sigMb2=0.55; //for square B
//static float kb0=2.0*pi/b_lattice,kb1=1.0,roub1=1.0,sigMb1=0.55,sigMb2=0.55; //for triangle B

static float para_c0=0.5,dynamic_mn=1.,dynamic_mc=1.,para_alpha=0.5,para_w=0.02;





static float axx=0.0*a_lattice,b0=sqrt(axx/dx/dx/dx/dt),b1=1.0*a_lattice,lam=a_lattice/sqrt(2.);  //noise

static int kd=0,mod=200;



static float f_l,f_s,phi_l,phi_s,con_l,con_s;
static float phimm=-1.0,phi_rat=0.618,phimin=0.0;
static double enERgy[TTT]; 
static double energySum=0.0;
static float krefl[L],krefm[M],krefn[N];
// ks provided by pfc_common.h
// min3 provided by pfc_common.h

static float frame_x=L*dx,frame_y=M*dx,frame_z=N*dx; 
static float nx=frame_x/a_lattice,ny=frame_y/a_lattice,nz=frame_z/a_lattice;
static int topmax=N-30,topmin=30;
// facepositon removed (was fixed-size with Tend)

// Forward declarations
static void INITIAL_C2(fftw_complex *dphi,fftw_complex *dcon,int dlocal_0_start, int dlocal_n0);
static void INITIAL_con(fftw_complex *dcon,int dlocal_n0); 
static void UPDATECOEF(fftw_complex *dphi,fftw_complex *dcon,int dlocal_n0);  
static void C2AK(float *dC2AHat, int dlocal_n0, int dlocal_0_start);
static void C2BK(float *dC2BHat, int dlocal_n0, int dlocal_0_start);
static void UPDATEINN3(fftw_complex *din1, fftw_complex *dphi,int dlocal_n0);
static void UPDATEINFMIX(fftw_complex *din2, fftw_complex *dcon,int dlocal_n0);
static void UPDATEINFCMIX(fftw_plan pea,fftw_plan peb,fftw_complex *din3,float *dC2AHat ,float *dC2BHat, fftw_complex *dphiHat, fftw_complex *dc2adenp,fftw_complex *dc2bdenp,fftw_complex *dc2aden,fftw_complex *dc2bden,fftw_complex *dcon,int dlocal_n0);
static void UPDATEINC2MIX(fftw_complex *dinn, fftw_complex *din1,fftw_complex *din2,fftw_complex *din3,int dlocal_n0);
static void UPDATEINC2MIXC(fftw_complex *dinc, fftw_complex *dphi, fftw_complex *dc2aden,fftw_complex *dc2bden,fftw_complex *dcon,int dlocal_n0);
static void MAINEQUATIONC2MIX(fftw_complex *dphiHat,fftw_complex *dconHat, fftw_complex *doutn,fftw_complex *doutc, int dlocal_n0, int dlocal_0_start); 
static void MAINEQUATIONC2_NOI(fftw_complex *dphiHat,fftw_complex *dconHat, fftw_complex *doutn,fftw_complex *doutc,fftw_complex *dnoisf, int dlocal_n0, int dlocal_0_start); 
static void UPDATEINC2GRA(fftw_plan pe1,fftw_plan pe2,fftw_plan pe3, fftw_complex *dgra1,fftw_complex *dgra2,fftw_complex *dgra3, fftw_complex *dgrao1,fftw_complex *dgrao2,fftw_complex *dgrao3,fftw_complex *dconHat, fftw_complex *dgrac, float *dkrefl,float *dkrefm,float *dkrefn,int dlocal_n0);
static void DAT(fftw_complex *dphi,fftw_complex *dcon,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start);
static void Gauss(fftw_complex *dnois,int dlocal_n0,int dlocal_0_start,int dmyid);
static void ENERGYC2(fftw_plan pe, fftw_complex *dphi, fftw_complex *dcon,fftw_complex *dconHat,fftw_complex *din3,fftw_complex *din,fftw_complex *dgrac,fftw_complex *dout,double *denei,double &dEEsum,double &dDSsum,int dlocal_n0,int dlocal_0_start,int dmyid, int dnumprocs);
static void INITIALREADPHI(fftw_complex *dphi, fftw_complex *dcon,int dlocal_n0,int dmyid);
static void atompositionCFGC2(fftw_complex *dphi,int dlocal_n0,int dlocal_0_start, int dmyid,int dnumprocs);
static void GR(float *drbuf_x, float *drbuf_y, float *drbuf_z, int dsize_buf);
static void Q6_GLO(float *drbuf_x, float *drbuf_y, float *drbuf_z, int dsize_buf);
static void OUTPUTPHI(fftw_complex *dphi,fftw_complex *dcon,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start);
static void SLICE(fftw_complex *dphi,fftw_complex *dcon,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start);
static void SLICEHat(float *dphi,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start,char *slicename);
static void LINEF(fftw_complex *dphi,fftw_complex *dcon,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start,char *LINEname);
static void OUTPUT(fftw_complex *dphi,fftw_complex *dcon,double *denei,float *dphiave,float *dconave,float *denei_ave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs);
static void interface_energy(fftw_complex *dphi,double *denei,float *dphiave,float *denei_ave,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start);
static void sVTKwriteBianry(fftw_complex *dphi, int dlocal_n0, int dlocal_0_start, int dmyid, int dnumprocs);
static void PVTKWRITE(int dlocal_n0, int dnumprocs, int dstartz, int dendz, char *dfilename);
static void sVTKwriteBianryc(fftw_complex *dcon, int dlocal_n0, int dlocal_0_start, int dmyid, int dnumprocs);
static void PVTKWRITEC(int dlocal_n0, int dnumprocs, int dstartz, int dendz, char *dfilename);
static void SMOOTHP(float *dphi,float *dphiave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename);
static void SMOOTHC(float *dcon,float *dconave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename);
static void SMOOTHE(float *dphi,float *denei_ave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename);
static void SMOOTHCON(fftw_complex *dcon,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs);
static void INTERFACEM(fftw_complex *dphi,int dicd,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,int &dintertop,int &dinterlow);
static void inirotationzb(float vtr2[9], float rotmatrix[9]);
static void ROT_XYZ(float dang_base[9], float dang[9], float *dtheta, char *daxis);
static void READPHIVTIN(fftw_complex *dphi, int dlocal_n0,int dmyid);
static void READPHIVTIC(fftw_complex *dcon, int dlocal_n0,int dmyid);

static void PrintUsage(const char* argv0) {
    cerr
        << "Usage: " << (argv0 ? argv0 : "pfc-exp-cli") << " --model cvd [options]\n"
        << "Options:\n"
        << "  --u0 <float>        Density/phi mean (default: 0.01)\n"
        << "  --con0 <float>      Concentration mean (default: 0.2)\n"
        << "  --sig <float>       Kernel decay sigma (default: 0.03)\n"
        << "  --dt <float>        Time step (default: 0.05)\n"
        << "  --dx <float>        Space step (default: 0.125)\n"
        << "  --steps <int>       Total iterations (default: 2000)\n"
        << "  --mod <int>         Output/analysis interval (default: 200)\n"
        << "  --seed <uint>       Random seed base (default: 20200604)\n"
        << "  --outdir <path>     Output directory (default: current dir)\n"
        << "  -h, --help          Show this help\n";
}

static int ParseArgs(int argc, char** argv, int myid) {
    for (int i = 1; i < argc; ++i) {
        const char* raw = argv[i];
        if (!raw) continue;
        std::string_view arg(raw);

        auto require_value = [&](const char* opt) -> const char* {
            if (i + 1 >= argc) {
                if (myid == 0) cerr << "Missing value after " << opt << "\n";
                return nullptr;
            }
            return argv[++i];
        };

        if (arg == "-h" || arg == "--help") {
            if (myid == 0) PrintUsage(argv[0]);
            MPI_Finalize();
return 0;
        }
        if (arg == "--u0" || PFC_HasPrefix(arg, "--u0=")) {
            const char* v = (arg == "--u0") ? require_value("--u0") : (raw + 5);
            if (!PFC_ParseFloatArg(v, u0)) goto bad;
            continue;
        }
        if (arg == "--con0" || PFC_HasPrefix(arg, "--con0=")) {
            const char* v = (arg == "--con0") ? require_value("--con0") : (raw + 7);
            if (!PFC_ParseFloatArg(v, con0)) goto bad;
            continue;
        }
        if (arg == "--sig" || PFC_HasPrefix(arg, "--sig=")) {
            const char* v = (arg == "--sig") ? require_value("--sig") : (raw + 6);
            if (!PFC_ParseFloatArg(v, sig)) goto bad;
            continue;
        }
        if (arg == "--dt" || PFC_HasPrefix(arg, "--dt=")) {
            const char* v = (arg == "--dt") ? require_value("--dt") : (raw + 5);
            if (!PFC_ParseFloatArg(v, dt)) goto bad;
            continue;
        }
        if (arg == "--dx" || PFC_HasPrefix(arg, "--dx=")) {
            const char* v = (arg == "--dx") ? require_value("--dx") : (raw + 5);
            if (!PFC_ParseFloatArg(v, dx)) goto bad;
            continue;
        }
        if (arg == "--steps" || PFC_HasPrefix(arg, "--steps=")) {
            const char* v = (arg == "--steps") ? require_value("--steps") : (raw + 8);
            if (!PFC_ParseIntArg(v, total_steps)) goto bad;
            continue;
        }
        if (arg == "--mod" || PFC_HasPrefix(arg, "--mod=")) {
            const char* v = (arg == "--mod") ? require_value("--mod") : (raw + 6);
            if (!PFC_ParseIntArg(v, mod)) goto bad;
            continue;
        }
        if (arg == "--seed" || PFC_HasPrefix(arg, "--seed=")) {
            const char* v = (arg == "--seed") ? require_value("--seed") : (raw + 7);
            if (!PFC_ParseUintArg(v, base_seed)) goto bad;
            continue;
        }
        if (arg == "--outdir" || PFC_HasPrefix(arg, "--outdir=")) {
            const char* v = (arg == "--outdir") ? require_value("--outdir") : (raw + 9);
            if (!v || *v == '\0') goto bad;
            output_dir = v;
            continue;
        }
        // Skip unknown args silently (may be from dispatcher)
        continue;
    }
    if (total_steps <= 0 || dt <= 0.0f || dx <= 0.0f || mod <= 0) {
        if (myid == 0) cerr << "Invalid parameter values\n";
        goto bad;
    }
    return 0;
bad:
    if (myid == 0) PrintUsage(argv[0]);
    return 2;
}

static void RecomputeDerivedParams() {
    frame_x = L * dx;
    frame_y = M * dx;
    frame_z = N * dx;
    nx = frame_x / a_lattice;
    ny = frame_y / a_lattice;
    nz = frame_z / a_lattice;
    b0 = (dx > 0.0f && dt > 0.0f && axx >= 0.0f) ? sqrt(axx / dx / dx / dx / dt) : 0.0f;
}

int run_cvd(int argc, char **argv) {


double EEsum,bufEEsum,EEsumTotal,DSsum,EE_sum=1.0E8, *enei=0 ;
float *phiHatPow=0, *phiHati=0;
float *enei_ave,*phiave,*conave;
float *C2Hat=0,*C2AHat=0,*C2BHat=0;
static char *Pname="phi_",*Pnamec="con_", *model="pr", *Lname="Line";
char filename[20];
fftw_complex *phi=0, *con=0, *phiHat=0, *conHat=0,*in=0,*inn=0,*inc=0, *gra1=0, *grao1=0,*gra2=0, *grao2=0,*gra3=0, *grao3=0,*grac=0;
fftw_complex *in1=0,*in2=0,*in3=0, *out=0,*outn=0,*outc=0, *nois=0, *noisf=0, *V=0,*c2bden=0,*c2adenp=0,*c2aden=0,*c2bdenp=0;
fftw_plan p1, p2,p3,p4, p5,p6,p7,p8,p9,p10,p11,pe12,pe13,pe14;
ptrdiff_t alloc_local,local_n0,local_0_start;
int myid,numprocs;

MPI_Status stat;
float  Vel,bufVel,VelTotal;
int icd=25;
int intertop=L-1-icd, interlow=1+icd;
//int dmyid;
MPI_Init(&argc, &argv);
fftw_mpi_init();
MPI_Comm_rank(MPI_COMM_WORLD, & myid);
MPI_Comm_size(MPI_COMM_WORLD, & numprocs);

    // Parse CLI arguments
    for (int i = 1; i < argc; ++i) {
        if (!argv[i]) continue;
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            PrintUsage(argv[0]);
            return 0;
        }
    }
    const int parse_rc = ParseArgs(argc, argv, myid);
    if (parse_rc != 0) {
        MPI_Finalize();
        return (parse_rc == 1) ? 0 : parse_rc;
    }
    RecomputeDerivedParams();

    if (!PFC_SetupOutputDir(output_dir, myid)) {
        MPI_Finalize();
        return 2;
    }

    ofstream checkpoint_log;
    if (myid == 0) {
        cerr << "CVD Run config: "
             << "u0=" << u0 << " con0=" << con0 << " sig=" << sig
             << " dt=" << dt << " dx=" << dx
             << " steps=" << total_steps << " mod=" << mod
             << " seed=" << base_seed
             << " grid=(" << L << "," << M << "," << N << ")"
             << " outdir=" << (output_dir.empty() ? "." : output_dir)
             << endl;

        ofstream run_cfg("run_config.txt", ios::out | ios::trunc);
        if (run_cfg.is_open()) {
            run_cfg << "model cvd\n";
            run_cfg << "u0 " << u0 << "\n";
            run_cfg << "con0 " << con0 << "\n";
            run_cfg << "sig " << sig << "\n";
            run_cfg << "dt " << dt << "\n";
            run_cfg << "dx " << dx << "\n";
            run_cfg << "steps " << total_steps << "\n";
            run_cfg << "mod " << mod << "\n";
            run_cfg << "seed " << base_seed << "\n";
            run_cfg << "L " << L << "\n";
            run_cfg << "M " << M << "\n";
            run_cfg << "N " << N << "\n";
            run_cfg.close();
        }

        checkpoint_log.open("checkpoint_timestamps.txt", ios::out | ios::trunc);
        if (checkpoint_log.is_open()) {
            checkpoint_log << "# Checkpoint timestamps - Step number and milliseconds since epoch" << endl;
        }
    }


alloc_local = fftw_mpi_local_size_3d(L, M, N,MPI_COMM_WORLD,&local_n0, &local_0_start);
phi = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
con = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
phiHat = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
conHat = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);

gra1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
gra2 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
gra3 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
grao1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
grao2 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
grao3 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
grac = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);

in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
inn = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
inc = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
in1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
in2 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
in3 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
outn = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
outc = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
phiHatPow = (float *) fftw_malloc(sizeof(float) * alloc_local);
phiHati = (float *) fftw_malloc(sizeof(float) * alloc_local);
enei = (double *) fftw_malloc(sizeof(double) * alloc_local);
C2Hat = (float *) fftw_malloc(sizeof(float) * alloc_local);
C2AHat = (float *) fftw_malloc(sizeof(float) * alloc_local);
C2BHat = (float *) fftw_malloc(sizeof(float) * alloc_local);
nois = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
noisf = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
enei_ave = (float *) fftw_malloc(sizeof(float) * alloc_local);
phiave = (float *) fftw_malloc(sizeof(float) * alloc_local);
conave = (float *) fftw_malloc(sizeof(float) * alloc_local);
c2bden= (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
c2adenp= (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
c2aden= (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);
c2bdenp= (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * alloc_local);

p1 = fftw_mpi_plan_dft_3d(L, M, N, phi, phiHat,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE);
p2 = fftw_mpi_plan_dft_3d(L, M, N, phiHat, phi,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE);
p3 = fftw_mpi_plan_dft_3d(L, M, N, in, out,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE);
p4 = fftw_mpi_plan_dft_3d(L, M, N, out, in,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE);
p5 = fftw_mpi_plan_dft_3d(L, M, N, nois, noisf,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE);
p6 = fftw_mpi_plan_dft_3d(L, M, N, con, conHat,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE);
p7 = fftw_mpi_plan_dft_3d(L, M, N, conHat, con,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE);
p8 = fftw_mpi_plan_dft_3d(L, M, N, inn, outn,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE);
p9 = fftw_mpi_plan_dft_3d(L, M, N, inc, outc,MPI_COMM_WORLD,FFTW_FORWARD,FFTW_ESTIMATE);
p10 = fftw_mpi_plan_dft_3d(L, M, N, c2adenp, c2aden,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE);
p11 = fftw_mpi_plan_dft_3d(L, M, N, c2bdenp, c2bden,MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE);

pe12 = fftw_mpi_plan_dft_3d(L, M, N, gra1, grao1, MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE);
pe13 = fftw_mpi_plan_dft_3d(L, M, N, gra2, grao2, MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE);
pe14 = fftw_mpi_plan_dft_3d(L, M, N, gra3, grao3, MPI_COMM_WORLD,FFTW_BACKWARD,FFTW_ESTIMATE);


for (int l=0; l <L/2+1; l ++) 	krefl[l]=l*2.0*pi/(L*dx);	
for (int l=L/2+1; l <L; l ++)	krefl[l]=(l-L)*2.0*pi/(L*dx);

for (int m=0; m <M/2+1; m ++) 	krefm[m]=m*2.0*pi/(M*dx);	
for (int m=M/2+1; m <M; m ++)	krefm[m]=(m-M)*2.0*pi/(M*dx);

for (int n=0; n <N/2+1; n ++) 	krefn[n]=n*2.0*pi/(N*dx);	
for (int n=N/2+1; n <N; n ++)	krefn[n]=(n-N)*2.0*pi/(N*dx);



kd=-1;

//for(int tt=0; tt<1; tt++){      //*
    
//sig=0.01*tt;

C2AK(C2AHat,local_n0,local_0_start);
C2BK(C2BHat,local_n0,local_0_start);



//for(int ttt=50;ttt<(TTT+1);ttt++){ // ttt

	//con0=1.0E-10+0.01*ttt;

INITIAL_C2(phi,con,local_0_start,local_n0);  

//READPHIVTIN(phi,local_n0,myid);
//READPHIVTIC(con,local_n0,myid);



kd=-1;




for(int t=0;t<total_steps; t++){  //*
	kd=kd+1;

    // Report progress to GUI
    if (myid == 0) {
        fprintf(stderr, "PFC_PROGRESS step=%d total=%d\n", t+1, total_steps);
        fflush(stderr);
    }


	//for(int tttt=0;tttt<50;tttt++){ //***
	UPDATEINN3(in1,phi,local_n0);
    UPDATEINFMIX(in2,con,local_n0);
	fftw_execute(p1);//phi-->phiHat     
    UPDATEINFCMIX(p10,p11,in3,C2AHat,C2BHat,phiHat,c2adenp,c2bdenp,c2aden,c2bden,con,local_n0);
    UPDATEINC2MIX(inn,in1,in2,in3,local_n0);
	UPDATEINC2MIXC(inc,phi,c2aden,c2bden,con,local_n0);
	fftw_execute(p8);//out(non_phi^3Hat)
	fftw_execute(p9);//out(non_con^3Hat)
	fftw_execute(p6);//con-->conHat
    MAINEQUATIONC2MIX(phiHat,conHat,outn,outc,local_n0,local_0_start);

   // Gauss(nois,local_n0,local_0_start,myid);   //����
   // fftw_execute(p5);//nois--->noisf           //����
   // MAINEQUATIONC2_NOI(phiHat,conHat,outn,outc,noisf,local_n0,local_0_start);

    UPDATEINC2GRA(pe12, pe13, pe14, gra1, gra2, gra3, grao1, grao2, grao3, conHat, grac, krefl, krefm, krefn, local_n0);
 
	fftw_execute(p2);//phiHat-->phi
	fftw_execute(p7);//conHat-->con
    UPDATECOEF(phi,con,local_n0);//phi<--XXX/N^3

//-------------------------------------------------------------------------------------------------------
	UPDATEINFCMIX(p10,p11,in3,C2AHat,C2BHat,phiHat,c2adenp,c2bdenp,c2aden,c2bden,con,local_n0);
//	         } //***tttt

	if ((t % mod) == 0) {
		kd = t;
		sVTKwriteBianry(phi, local_n0, local_0_start, myid, numprocs);
		PVTKWRITE(local_n0, numprocs, 0, N, Pname);
		sVTKwriteBianryc(con, local_n0, local_0_start, myid, numprocs);
		PVTKWRITEC(local_n0, numprocs, 0, N, Pnamec);

		if (myid == 0 && checkpoint_log.is_open()) {
		    const auto now = std::chrono::system_clock::now();
		    const auto timestamp_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
		        now.time_since_epoch()).count();
		    checkpoint_log << "step " << t << " timestamp_ms " << timestamp_ms << endl;
		    checkpoint_log.flush();
		}
	}

      ENERGYC2(p4, phi, con,conHat,in3,in,grac, out,enei,EEsum,DSsum,local_n0,local_0_start,myid, numprocs);



/*if(fabs(energySum-EE_sum)>1.0E-5)
	{  
	    EE_sum=energySum;		
	}
	else 
	{
		break;
	}*/
} //* t
//��������Сֵ*************************************


 //enERgy[ttt]=energySum;

             // char *file2;
             // file2="phi_";
           // if (myid==0) {PVTKWRITE(local_n0,numprocs,file2);} 
             //  sVTKwriteBianry(phi,local_n0,local_0_start,myid,numprocs,file2);

            //  file2="con_";
           // if (myid==0) { PVTKWRITE(local_n0,numprocs,file2);} 
             //  sVTKwriteBianry(con,local_n0,local_0_start,myid,numprocs,file2); 




//} //ttt con




 //  sprintf(filename,"%s%d%s", "B_",int(100*sig),".txt");

	//ofstream outf(filename,ios_base::out);
 //   outf.setf(ios_base::fixed,ios_base::floatfield);
 //   outf.precision(5);
	////outf<<"sig="<<sig<<endl;
 //   for(int lll=50;lll<TTT+1;lll++) { outf<<0.01*lll<<" "<<enERgy[lll]<<endl;}
 //   outf.close(); 

//energy out<<











//}  //* sig







          

fftw_destroy_plan(p1);
fftw_destroy_plan(p2);
fftw_destroy_plan(p3);
fftw_destroy_plan(p4);
fftw_destroy_plan(p5);	
fftw_destroy_plan(p6);
fftw_destroy_plan(p7);
fftw_destroy_plan(p8);
fftw_destroy_plan(p9);
fftw_destroy_plan(p10);	
fftw_destroy_plan(p11);	
fftw_destroy_plan(pe12);
fftw_destroy_plan(pe13);
fftw_destroy_plan(pe14);

free(phi);
free(phiHat);
free(in);
free(out);
free(C2Hat);
free(C2AHat);
free(C2BHat);
free(nois);
free(noisf);
free(enei);
free(phiHatPow);
free(phiHati);
free(enei_ave);
free(phiave);
free(con);
free(conHat);
free(inn);
free(inc);
free(in1);
free(in2);
free(in3);
free(outn);
free(outc);
free(conave);
free(c2bden);
free(c2adenp);
free(c2aden);
free(c2bdenp);
free(gra1);
free(gra2);
free(gra3);
free(grao1);
free(grao2);
free(grao3);
free(grac);

MPI_Finalize();
return 0;
}

static void INITIAL_con(fftw_complex *dcon,int dlocal_n0){



   for(int l=0; l <dlocal_n0; l ++){
	   for (int m=0; m <M; m ++){
		for (int n=0; n <N; n ++) {

                dcon[n+N*(m+M*l)][0]=con0;
		dcon[n+N*(m+M*l)][1]=1.0E-10;

                                   }}}

}










static void INITIAL_C2(fftw_complex *dphi,fftw_complex *dcon,int dlocal_0_start, int dlocal_n0){
    float x, y, z;
	float A1=0.15,A2=0.106; 
	 float r1,w1,w2,w3;
     w1=0.;w2=0.;w3=0.*pi/180;


   for(int l=0; l <dlocal_n0; l ++){
	   for (int m=0; m <M; m ++){
		for (int n=0; n <N; n ++) {
                  x=(dlocal_0_start-L/2+l)*dx;
                  y=(m-M/2)*dx;
                  z=(n-N/2)*dx;
                  dphi[n+N*(m+M*l)][1]=0;
				  dcon[n+N*(m+M*l)][1]=1.0E-10;
		//-----------------------------------------------------------------------------------------------------
	
            //����
			//dphi[n+N*(m+M*l)][0]=u0+A1*(cos(ka0*l*dx)+2*cos(ka0*l*dx/2)*cos(ka0*sqrt(3.)*m*dx/2));
            
		    //�ķ�
               // dphi[n+N*(m+M*l)][0]=2*A1*(cos(ka0*x*cos(w3)+ka0*y*sin(w3))+cos(-x*ka0*sin(w3)+ka0*y*cos(w3)))+4*A2*cos(ka0*x*cos(w3)+ka0*y*sin(w3))*cos(-x*ka0*sin(w3)+y*ka0*cos(w3))+u0; 
		

		


//-----------------------------------------------------------------------------------------------------
		//for B
	
            //����
			//dphi[n+N*(m+M*l)][0]=u0+A1*(cos(kb0*l*dx)+2*cos(kb0*l*dx/2)*cos(kb0*sqrt(3.)*m*dx/2));
           
			//�ķ�
				  if (m < M/2) {
					  dphi[n + N*(m + M*l)][0] = 2 * A1*(cos(ka0*x*cos(w3) + ka0*y*sin(w3)) + cos(-x*ka0*sin(w3) + ka0*y*cos(w3))) + 4 * A2*cos(ka0*x*cos(w3) + ka0*y*sin(w3))*cos(-x*ka0*sin(w3) + y*ka0*cos(w3)) + u0;
					  dcon[n + N*(m + M*l)][0] = con0;
					  				  }
				  else {
					  dphi[n + N*(m + M*l)][0] = 2 * A1*(cos(kb0*x*cos(w3) + kb0*y*sin(w3)) + cos(-x*kb0*sin(w3) + kb0*y*cos(w3))) + 4 * A2*cos(kb0*x*cos(w3) + kb0*y*sin(w3))*cos(-x*kb0*sin(w3) + y*kb0*cos(w3)) + u0;
					  dcon[n + N*(m + M*l)][0] =1-con0;
				  }
       // dphi[n+N*(m+M*l)][0]=2*A1*(cos(kb0*x*cos(w3)+kb0*y*sin(w3))+cos(-x*kb0*sin(w3)+kb0*y*cos(w3)))+4*A2*cos(kb0*x*cos(w3)+kb0*y*sin(w3))*cos(-x*kb0*sin(w3)+y*kb0*cos(w3))+u0; 	

       		//dcon[n+N*(m+M*l)][0]=con0;
        //-----------------------------------------------------------------------------------------------------


			}
	   }
	}  




}

static void UPDATECOEF(fftw_complex *dphi,fftw_complex *dcon,int dlocal_n0) {
        for (int l=0; l <dlocal_n0; l ++){
		for (int m=0; m <M; m ++){
			for (int n=0; n <N; n ++){
		           dphi[n+N*(m+M*l)][0]=dphi[n+N*(m+M*l)][0]/(L*M*N);
                   dphi[n+N*(m+M*l)][1]=0;
				   dcon[n+N*(m+M*l)][0]=dcon[n+N*(m+M*l)][0]/(L*M*N);
                   dcon[n+N*(m+M*l)][1]=0.0;
			}
		}
	}
}


static void C2AK(float *dC2AHat , int dlocal_n0, int dlocal_0_start){
	int lnew;
         float kc1,aa1,aa2,bb1,bb2;
         aa1=pow(sig/sigMa1,2);  aa2=pow(sig/sigMa2,2);
                bb1=1./(2.*gaa1*gaa1); bb2=1./(2.*gaa2*gaa2);
               kc1=((bb1*ka1-bb2*ka2)+sqrt(pow(bb1*ka1-bb2*ka2,2)-(bb1-bb2)*(aa1-aa2+bb1*ka1*ka1-bb2*ka2*ka2)))/(bb1-bb2);
               if (gaa1==gaa2) kc1=(aa2-aa1+bb1*(ka2*ka2-ka1*ka1))/(2*bb1*(ka2-ka1)); 


	for (int l=0; l <dlocal_n0; l ++){
		for (int m=0; m <M; m ++){
		   for (int n=0; n <N; n ++){
             lnew=l+dlocal_0_start;

             if(sqrt(ks(lnew,m,n,krefl,krefm,krefn))<=kc1) dC2AHat[n+N*(m+M*l)]=exp(-pow(sig/sigMa1,2))*exp(-pow( sqrt(ks(lnew,m,n,krefl,krefm,krefn))-ka1,2)/(2*gaa1*gaa1));
             if(sqrt(ks(lnew,m,n,krefl,krefm,krefn))>kc1) dC2AHat[n+N*(m+M*l)]=exp(-pow(sig/sigMa2,2))*exp(-pow(sqrt(ks(lnew,m,n,krefl,krefm,krefn))-ka2,2)/(2*gaa2*gaa2));

	
			// dC2AHat[n+N*(m+M*l)]=exp(-pow(sig/sigMa1,2))*exp(-pow( sqrt(ks(lnew,m,n,krefl,krefm,krefn))-ka1,2)/(2*gaa1*gaa1)); //+exp(-pow(sig/sigMa2,2))*exp(-pow( sqrt(ks(lnew,m,n,krefl,krefm,krefn))-ka2,2)/(2*gaa2*gaa2) )+epsilona;
	
		   }
		}
	}
}


static void C2BK(float *dC2BHat, int dlocal_n0, int dlocal_0_start){
	int lnew;

         float kc2,aa1,aa2,bb1,bb2;
         aa1=pow(sig/sigMb1,2);  aa2=pow(sig/sigMb2,2);
                bb1=1./(2.*gab1*gab1); bb2=1./(2.*gab2*gab2); 
                kc2=((bb1*kb1-bb2*kb2)+sqrt(pow(bb1*kb1-bb2*kb2,2)-(bb1-bb2)*(aa1-aa2+bb1*kb1*kb1-bb2*kb2*kb2)))/(bb1-bb2);
              if (gab1==gab2) kc2=(aa2-aa1+bb1*(kb2*kb2-kb1*kb1))/(2*bb1*(kb2-kb1)); 

	for (int l=0; l <dlocal_n0; l ++){
		for (int m=0; m <M; m ++){
		   for (int n=0; n <N; n ++){
             lnew=l+dlocal_0_start;

           if(sqrt(ks(lnew,m,n,krefl,krefm,krefn))<=kc2) dC2BHat[n+N*(m+M*l)]=exp(-pow(sig/sigMb1,2))*exp(-pow(sqrt(ks(lnew,m,n,krefl,krefm,krefn))-kb1,2)/(2*gab1*gab1));
           if(sqrt(ks(lnew,m,n,krefl,krefm,krefn))>kc2)  dC2BHat[n+N*(m+M*l)]=exp(-pow(sig/sigMb2,2))*exp(-pow(sqrt(ks(lnew,m,n,krefl,krefm,krefn))-kb2,2)/(2*gab2*gab2));

	
	// dC2BHat[n+N*(m+M*l)]=exp(-pow(sig/sigMb1,2))*exp(-pow( sqrt(ks(lnew,m,n,krefl,krefm,krefn))-kb1,2)/(2*gab1*gab1) ); //+exp(-pow(sig/sigMb2,2))*exp(-pow( sqrt(ks(lnew,m,n,krefl,krefm,krefn))-kb2,2)/(2*gab2*gab2) )+epsilonb;

	
		   }
		}
	}
}

static void UPDATEINN3(fftw_complex *din1, fftw_complex *dphi,int dlocal_n0){
        for (int l=0; l <dlocal_n0; l ++){
	        for (int m=0; m <M; m ++){
			for (int n=0; n <N; n ++){
			   din1[n+N*(m+M*l)][0]=-niu*pow(dphi[n+N*(m+M*l)][0],2)/2.0+ksi*pow(dphi[n+N*(m+M*l)][0],3)/3.0;
               din1[n+N*(m+M*l)][1]=0.0; //-niu*pow(dphi[n+N*(m+M*l)][1],2)/2.0+ksi*pow(dphi[n+N*(m+M*l)][1],3)/3.0;
			 //  din[n+N*(m+M*l)][0]=para_w*(dcon[n+N*(m+M*l)][0]*log(dcon[n+N*(m+M*l)][0]/para_c0)+(1-dcon[n+N*(m+M*l)][0])*log((1-dcon[n+N*(m+M*l)][0])/(1-para_c0)));
			}
		}
	}

}

static void UPDATEINFMIX(fftw_complex *din2, fftw_complex *dcon,int dlocal_n0){
        for (int l=0; l <dlocal_n0; l ++){
	        for (int m=0; m <M; m ++){
			for (int n=0; n <N; n ++){
			  // din[n+N*(m+M*l)][0]=-niu*pow(dphi[n+N*(m+M*l)][0],2)/2+ksi*pow(dphi[n+N*(m+M*l)][0],3)/3;
              // din[n+N*(m+M*l)][1]=-niu*pow(dphi[n+N*(m+M*l)][1],2)/2+ksi*pow(dphi[n+N*(m+M*l)][1],3)/3;
			   din2[n+N*(m+M*l)][0]=para_w*(dcon[n+N*(m+M*l)][0]*log(dcon[n+N*(m+M*l)][0]/para_c0)+(1.0-dcon[n+N*(m+M*l)][0])*log((1.0-dcon[n+N*(m+M*l)][0])/(1.0-para_c0)));
			   din2[n+N*(m+M*l)][1]=0.0; //para_w*(dcon[n+N*(m+M*l)][1]*log(1.0E-10+fabs(dcon[n+N*(m+M*l)][1]/para_c0))+(1.0-dcon[n+N*(m+M*l)][1])*log(1.0E-10+fabs((1.0-dcon[n+N*(m+M*l)][1])/(1.0-para_c0))));
			}
		}
	}

}


static void UPDATEINFCMIX(fftw_plan pea,fftw_plan peb,fftw_complex *din3,float *dC2AHat ,float *dC2BHat , fftw_complex *dphiHat, fftw_complex *dc2adenp,fftw_complex *dc2bdenp,fftw_complex *dc2aden,fftw_complex *dc2bden,fftw_complex *dcon,int dlocal_n0){
    
	for (int l=0; l <dlocal_n0; l ++){
	   for (int m=0; m <M; m ++){
	       for (int n=0; n <N; n ++){
                    
		   dc2adenp[n+N*(m+M*l)][0]=dC2AHat[n+N*(m+M*l)]*dphiHat[n+N*(m+M*l)][0];
           dc2adenp[n+N*(m+M*l)][1]=dC2AHat[n+N*(m+M*l)]*dphiHat[n+N*(m+M*l)][1];
		   dc2bdenp[n+N*(m+M*l)][0]=dC2BHat[n+N*(m+M*l)]*dphiHat[n+N*(m+M*l)][0];
           dc2bdenp[n+N*(m+M*l)][1]=dC2BHat[n+N*(m+M*l)]*dphiHat[n+N*(m+M*l)][1];
				}
			}
		}  
	
    fftw_execute(pea); //enerHat-->ener
	fftw_execute(peb); //enerHat-->ener

	for (int l=0; l <dlocal_n0; l ++){
	     for (int m=0; m <M; m ++){
			for (int n=0; n <N; n ++){
               dc2aden[n+N*(m+M*l)][0]=dc2aden[n+N*(m+M*l)][0]/(L*M*N);
               dc2aden[n+N*(m+M*l)][1]=dc2aden[n+N*(m+M*l)][1]/(L*M*N);
			   dc2bden[n+N*(m+M*l)][0]=dc2bden[n+N*(m+M*l)][0]/(L*M*N);
			   dc2bden[n+N*(m+M*l)][1]=dc2bden[n+N*(m+M*l)][1]/(L*M*N);

			   din3[n+N*(m+M*l)][0]=-(1.0-3.0*pow(dcon[n+N*(m+M*l)][0],2)+2.0*pow(dcon[n+N*(m+M*l)][0],3))*(dc2aden[n+N*(m+M*l)][0])-(1.0-3.0*pow((1.0-dcon[n+N*(m+M*l)][0]),2)+2.0*pow((1.0-dcon[n+N*(m+M*l)][0]),3))*(dc2bden[n+N*(m+M*l)][0]);
			   din3[n+N*(m+M*l)][1]=0.0;//-(1.0-3.0*pow(dcon[n+N*(m+M*l)][1],2)+2.0*pow(dcon[n+N*(m+M*l)][1],3))*(dc2aden[n+N*(m+M*l)][1])-(1.0-3.0*pow((1.0-dcon[n+N*(m+M*l)][1]),2)+2.0*pow((1.0-dcon[n+N*(m+M*l)][1]),3))*(dc2bden[n+N*(m+M*l)][1]);
			
	        }
	    }
	}

}

static void UPDATEINC2MIX(fftw_complex *dinn, fftw_complex *din1,fftw_complex *din2,fftw_complex *din3,int dlocal_n0){
        for (int l=0; l <dlocal_n0; l ++){
	        for (int m=0; m <M; m ++){
			for (int n=0; n <N; n ++){
				
			   dinn[n+N*(m+M*l)][0]=din1[n+N*(m+M*l)][0]+din2[n+N*(m+M*l)][0]+din3[n+N*(m+M*l)][0];
               dinn[n+N*(m+M*l)][1]=0.0; //din1[n+N*(m+M*l)][1]+din2[n+N*(m+M*l)][1]+din3[n+N*(m+M*l)][1];
		  //     dinn[n+N*(m+M*l)][0]=din1[n+N*(m+M*l)][0]+din3[n+N*(m+M*l)][0];
          //     dinn[n+N*(m+M*l)][1]=din1[n+N*(m+M*l)][1]+din3[n+N*(m+M*l)][1];
			}
		}
	}

}


static void UPDATEINC2MIXC(fftw_complex *dinc, fftw_complex *dphi, fftw_complex *dc2aden,fftw_complex *dc2bden,fftw_complex *dcon,int dlocal_n0){
    float cmix1,cmix11,cmix2,cmix22,cmix3,cmix33,cmix4,cmix44;    

	for (int l=0; l <dlocal_n0; l ++){
	        for (int m=0; m <M; m ++){
			for (int n=0; n <N; n ++){
				
			   cmix1=para_w*(dphi[n+N*(m+M*l)][0]+1.0)*log(dcon[n+N*(m+M*l)][0]/para_c0);
               cmix11=0.0; //para_w*(dphi[n+N*(m+M*l)][1]+1.0)*log(1.0E-10+fabs(dcon[n+N*(m+M*l)][1]/para_c0));
               cmix2=para_w*(dphi[n+N*(m+M*l)][0]+1.0)*log((1.0-dcon[n+N*(m+M*l)][0])/(1.0-para_c0));
               cmix22=0.0; //para_w*(dphi[n+N*(m+M*l)][1]+1.0)*log(1.0E-10+fabs((1.0-dcon[n+N*(m+M*l)][1])/(1.0-para_c0)));
			   cmix3=3.0*dphi[n+N*(m+M*l)][0]*(dc2aden[n+N*(m+M*l)][0]-dc2bden[n+N*(m+M*l)][0])*pow(dcon[n+N*(m+M*l)][0],2);
               cmix33=0.0; //3.0*dphi[n+N*(m+M*l)][1]*(dc2aden[n+N*(m+M*l)][1]-dc2bden[n+N*(m+M*l)][1])*pow(dcon[n+N*(m+M*l)][1],2);
			   cmix4=3.0*dphi[n+N*(m+M*l)][0]*(-1.0*dc2aden[n+N*(m+M*l)][0]+dc2bden[n+N*(m+M*l)][0])*dcon[n+N*(m+M*l)][0];
               cmix44=0.0; //3.0*dphi[n+N*(m+M*l)][1]*(-1.0*dc2aden[n+N*(m+M*l)][1]+dc2bden[n+N*(m+M*l)][1])*dcon[n+N*(m+M*l)][1];

			   dinc[n+N*(m+M*l)][0]=cmix1-cmix2-cmix3-cmix4;
               dinc[n+N*(m+M*l)][1]=0.0;//cmix11-cmix22-cmix33-cmix44;
			 
			}
		}
	}

}


static void UPDATEINC2GRA(fftw_plan pe1,fftw_plan pe2,fftw_plan pe3, fftw_complex *dgra1,fftw_complex *dgra2,fftw_complex *dgra3, fftw_complex *dgrao1,fftw_complex *dgrao2,fftw_complex *dgrao3,fftw_complex *dconHat, fftw_complex *dgrac, float *dkrefl,float *dkrefm,float *dkrefn,int dlocal_n0){
    float cmix1,cmix11,cmix2,cmix22,cmix3,cmix33,cmix4,cmix44;    

	for (int l=0; l <dlocal_n0; l ++){
	        for (int m=0; m <M; m ++){
			for (int n=0; n <N; n ++){

				dgra1[n+N*(m+M*l)][0]=1.0*dkrefl[l]*dconHat[n+N*(m+M*l)][1];
                dgra1[n+N*(m+M*l)][1]=-1.0*dkrefl[l]*dconHat[n+N*(m+M*l)][0];

				dgra2[n+N*(m+M*l)][0]=1.0*dkrefm[m]*dconHat[n+N*(m+M*l)][1];
                dgra2[n+N*(m+M*l)][1]=-1.0*dkrefm[m]*dconHat[n+N*(m+M*l)][0];

				dgra3[n+N*(m+M*l)][0]=1.0*dkrefn[n]*dconHat[n+N*(m+M*l)][1];
                dgra3[n+N*(m+M*l)][1]=-1.0*dkrefn[n]*dconHat[n+N*(m+M*l)][0];
				
			  
			}
		}
	}

	fftw_execute(pe1);
	fftw_execute(pe2);
	fftw_execute(pe3);

	for (int l=0; l <dlocal_n0; l ++){
	        for (int m=0; m <M; m ++){
			for (int n=0; n <N; n ++){

				dgrao1[n+N*(m+M*l)][0]=dgrao1[n+N*(m+M*l)][0]/L/M/N;
                dgrao1[n+N*(m+M*l)][1]=dgrao1[n+N*(m+M*l)][1]/L/M/N;

				dgrao2[n+N*(m+M*l)][0]=dgrao2[n+N*(m+M*l)][0]/L/M/N;
                dgrao2[n+N*(m+M*l)][1]=dgrao2[n+N*(m+M*l)][1]/L/M/N;

				dgrao3[n+N*(m+M*l)][0]=dgrao3[n+N*(m+M*l)][0]/L/M/N;
                dgrao3[n+N*(m+M*l)][1]=dgrao3[n+N*(m+M*l)][1]/L/M/N;

                dgrac[n+N*(m+M*l)][0]=pow(dgrao1[n+N*(m+M*l)][0],2)+pow(dgrao2[n+N*(m+M*l)][0],2)+pow(dgrao2[n+N*(m+M*l)][0],2);
				dgrac[n+N*(m+M*l)][1]=0.0;
				
			  
			}
		}
	}



}


static void MAINEQUATIONC2MIX(fftw_complex *dphiHat,fftw_complex *dconHat, fftw_complex *doutn,fftw_complex *doutc, int dlocal_n0, int dlocal_0_start){
	float Rhat1, Rhat2,Rhat3, Rhat4;
        int lnew;
	for (int l=0; l <dlocal_n0; l ++){
		for (int m=0; m <M; m ++){
		   for (int n=0; n <N; n ++){
                        lnew=l+dlocal_0_start;
	                Rhat1=dphiHat[n+N*(m+M*l)][0]-dynamic_mn*dt*ks(lnew,m,n,krefl,krefm,krefn)*doutn[n+N*(m+M*l)][0];
	                Rhat2=dphiHat[n+N*(m+M*l)][1]-dynamic_mn*dt*ks(lnew,m,n,krefl,krefm,krefn)*doutn[n+N*(m+M*l)][1];
                        dphiHat[n+N*(m+M*l)][0]=Rhat1/(1.0+dynamic_mn*dt*ks(lnew,m,n,krefl,krefm,krefn));
	                    dphiHat[n+N*(m+M*l)][1]=Rhat2/(1.0+dynamic_mn*dt*ks(lnew,m,n,krefl,krefm,krefn));

                    Rhat3=dconHat[n+N*(m+M*l)][0]-dynamic_mc*dt*ks(lnew,m,n,krefl,krefm,krefn)*doutc[n+N*(m+M*l)][0];
	                Rhat4=dconHat[n+N*(m+M*l)][1]-dynamic_mc*dt*ks(lnew,m,n,krefl,krefm,krefn)*doutc[n+N*(m+M*l)][1];
                        dconHat[n+N*(m+M*l)][0]=Rhat3/(1.0+dynamic_mc*dt*para_alpha*ks(lnew,m,n,krefl,krefm,krefn)*ks(lnew,m,n,krefl,krefm,krefn));
	                    dconHat[n+N*(m+M*l)][1]=Rhat4/(1.0+dynamic_mc*dt*para_alpha*ks(lnew,m,n,krefl,krefm,krefn)*ks(lnew,m,n,krefl,krefm,krefn));


			}
		}
	}
}

static void Gauss(fftw_complex *dnois,int dlocal_n0,int dlocal_0_start, int dmyid){
   float iseed_radia=20.;
   float r,v1,v2,f,ra,x,y,z;
		 srand((unsigned int)time(0)*(dmyid+1));
for(int l=0; l <dlocal_n0; l ++){
 for (int m=0; m <M; m ++){
     for (int n=0; n <N; n ++) {
	//	 if(((dlocal_0_start+l+1-L/2)*(dlocal_0_start+l+1-L/2)+(m+1-M/2)*(m+1-M/2))<32*32){
                  x=(dlocal_0_start+l+1-L/2)*dx;
                  y=(m+1-M/2)*dx;
                  z=0.0;//(n+1-N/2)*dx;
                  ra=sqrt(x*x+y*y+z*z);
        
				  
        r=2.; 
        while (r>1) 
		{ 
         v1=2*float(rand())/RAND_MAX-1;
         v2=2*float(rand())/RAND_MAX-1;
         r=v1*v1+v2*v2;
		}
          f=sqrt(-2.*log(r)/r);
          dnois[n+N*(m+M*l)][0]=f*v1;
          dnois[n+N*(m+M*l)][1]=0.0;
	//	 }

}
}
}
}

static void MAINEQUATIONC2_NOI(fftw_complex *dphiHat,fftw_complex *dconHat, fftw_complex *doutn,fftw_complex *doutc,fftw_complex *dnoisf, int dlocal_n0, int dlocal_0_start){
	float Rhat1, Rhat2,Rhat3, Rhat4;
	float detrx,detrxn,detry,detryn;
        int lnew;
	for (int l=0; l <dlocal_n0; l ++){
		for (int m=0; m <M; m ++){
		   for (int n=0; n <N; n ++){
                        lnew=l+dlocal_0_start;
					detrx=0.0;
	                detrxn=0.0;
	                detry=0.0;
	                detryn=0.0;
					if(lnew==0){detrx=1.0;}
					if(lnew==L/2){detrxn=1.0;}
					if(m==0){detry=1.0;}
					if(m==M/2){detryn=1.0;}
					/*b0=sqrt(axx/dx/dx/dt);
					if(dsmooth_phi[n+N*(m+M*l)][0]<phi_solid){b0=0.0;}*/

	                Rhat1=dphiHat[n+N*(m+M*l)][0]-dynamic_mn*dt*ks(lnew,m,n,krefl,krefm,krefn)*doutn[n+N*(m+M*l)][0]+dt*b0*dnoisf[n+N*(m+M*l)][0]*sqrt(ks(lnew,m,n,krefl,krefm,krefn))*sqrt(abs(0.5+0.5*(detrx+detrxn)*(detry+detryn)))*(0.5-0.5*tanh((sqrt(ks(lnew,m,n,krefl,krefm,krefn))-2.0*pi)/0.1));
	                Rhat2=dphiHat[n+N*(m+M*l)][1]-dynamic_mn*dt*ks(lnew,m,n,krefl,krefm,krefn)*doutn[n+N*(m+M*l)][1]+dt*b0*dnoisf[n+N*(m+M*l)][1]*sqrt(ks(lnew,m,n,krefl,krefm,krefn))*sqrt(abs(0.5-0.5*(detrx+detrxn)*(detry+detryn)))*(0.5-0.5*tanh((sqrt(ks(lnew,m,n,krefl,krefm,krefn))-2.0*pi)/0.1));
                        dphiHat[n+N*(m+M*l)][0]=Rhat1/(1.0+dynamic_mn*dt*ks(lnew,m,n,krefl,krefm,krefn));
	                    dphiHat[n+N*(m+M*l)][1]=Rhat2/(1.0+dynamic_mn*dt*ks(lnew,m,n,krefl,krefm,krefn));

                    Rhat3=dconHat[n+N*(m+M*l)][0]-dynamic_mc*dt*ks(lnew,m,n,krefl,krefm,krefn)*doutc[n+N*(m+M*l)][0];//+dt*b0*dnoisf[n+N*(m+M*l)][0]*sqrt(ks(lnew,m,n,krefl,krefm,krefn))*sqrt(abs(0.5+0.5*(detrx+detrxn)*(detry+detryn)))*(0.5-0.5*tanh((sqrt(ks(lnew,m,n,krefl,krefm,krefn))-1.0*pi)/0.1));
	                Rhat4=dconHat[n+N*(m+M*l)][1]-dynamic_mc*dt*ks(lnew,m,n,krefl,krefm,krefn)*doutc[n+N*(m+M*l)][1];//+dt*b0*dnoisf[n+N*(m+M*l)][1]*sqrt(ks(lnew,m,n,krefl,krefm,krefn))*sqrt(abs(0.5-0.5*(detrx+detrxn)*(detry+detryn)))*(0.5-0.5*tanh((sqrt(ks(lnew,m,n,krefl,krefm,krefn))-1.0*pi)/0.1));
                        dconHat[n+N*(m+M*l)][0]=Rhat3/(1.0+dynamic_mc*dt*para_alpha*ks(lnew,m,n,krefl,krefm,krefn)*ks(lnew,m,n,krefl,krefm,krefn));
	                    dconHat[n+N*(m+M*l)][1]=Rhat4/(1.0+dynamic_mc*dt*para_alpha*ks(lnew,m,n,krefl,krefm,krefn)*ks(lnew,m,n,krefl,krefm,krefn));
			}
		}
	}
}


static void DAT(fftw_complex *dphi,fftw_complex *dcon,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start){
  char filename1[20],filename2[20],filename3[20],filename4[20];
  sprintf(filename4,"%s%d%s", "phi_", kd,".dat");


//OUTPUT phi.dat
  float *temp2, *rbuf2;
  if ( dmyid == 0) {rbuf2= (float *) fftw_malloc(sizeof(float)*L*M*N);} 
  temp2= (float *) fftw_malloc(sizeof(float) * M*N*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  for (int m=0; m<M; m++) for (int n=0; n <N; n ++) {temp2[n+N*(m+M*l)]=dphi[n+N*(m+M*l)][0];}

  MPI_Gather(temp2,N*M*dlocal_n0,MPI_FLOAT,rbuf2,N*M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
    ofstream outf4(filename4,ios_base::out);
    outf4.setf(ios_base::fixed,ios_base::floatfield);
    outf4.precision(4);
	outf4<<"variable=\"x\",\"y\",\"z\",\"function\""<<endl; //dat
	outf4<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",k="<<N<<",f=point"<<endl;//dat
    for(int n=0; n <N; n ++)for(int m=0; m <M; m ++)for(int l=0; l <L; l ++) {outf4<<rbuf2[n+N*(m+M*l)]<<endl;}
    outf4.close(); 
free(temp2);
free(rbuf2);   
  } 
}


static void ENERGYC2(fftw_plan pe, fftw_complex *dphi, fftw_complex *dcon,fftw_complex *dconHat,fftw_complex *din3,fftw_complex *din,fftw_complex *dgrac,fftw_complex *dout,double *denei,double &dEEsum,double &dDSsum,int dlocal_n0,int dlocal_0_start,int dmyid, int dnumprocs){
        dEEsum=0.; 
         dDSsum=0.;
        int lnew;
		float energy_mix,energy_ex,energy_id;
	
	for (int l=0; l <dlocal_n0; l ++){
	    for (int m=0; m <M; m ++){
		for (int n=0; n <N; n ++){
           

		   energy_id=pow(dphi[n+N*(m+M*l)][0],2)/2.0-niu*pow(dphi[n+N*(m+M*l)][0],3)/6.0+ksi*pow(dphi[n+N*(m+M*l)][0],4)/12.0;
           energy_mix=(dphi[n+N*(m+M*l)][0]+1.0)*para_w*(dcon[n+N*(m+M*l)][0]*log(dcon[n+N*(m+M*l)][0]/para_c0)+(1.0-dcon[n+N*(m+M*l)][0])*log((1.0-dcon[n+N*(m+M*l)][0])/(1.0-para_c0)));
		   energy_ex=0.5*dphi[n+N*(m+M*l)][0]*(din3[n+N*(m+M*l)][0])+para_alpha*dgrac[n+N*(m+M*l)][0];
	
		              denei[n+N*(m+M*l)]=energy_id+energy_mix+energy_ex;//����CPUÿһ���������        
		              dEEsum += denei[n+N*(m+M*l)];//*dx*dx*dx;                   
                      dDSsum += dphi[n+N*(m+M*l)][0]; //����CPU�ܶ����

	        }
	    }
	}


   double * rcounts,* rcounts1;
       rcounts=new double[dnumprocs];
	   rcounts1=new double[dnumprocs];

   MPI_Gather(&dEEsum,1,MPI_DOUBLE,rcounts ,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
   MPI_Gather(&dDSsum,1,MPI_DOUBLE,rcounts1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
 if(dmyid==0) {
  
  double sum_enei=0;
  double sum_dens=0;
  for(int i=0;i<dnumprocs;i++){
     sum_enei+= rcounts[i];
     sum_dens+=rcounts1[i];}
     energySum=sum_enei/L/M/N;

  
    ofstream wafile("energy.txt",ios_base::out|ios_base::app);
	wafile << kd << " " << sum_enei/L/M/N <<" "<<sum_dens/L/M/N  <<endl;
}


MPI_Bcast(&energySum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);









delete [] rcounts;
delete [] rcounts1;
rcounts=0;
rcounts1=0;
}

static void INITIALREADPHI(fftw_complex *dphi, fftw_complex *dcon,int dlocal_n0,int dmyid){
 char filenamen[20],filenamec[20];
 sprintf(filenamen, "%s%d%s%d%s", "phi_", 0,"_",dmyid,".txt");
 ifstream infile1(filenamen,ios_base::in);
 for(int n=0; n <N; n ++) for (int m=0; m <M; m ++) for (int l=0; l <dlocal_n0; l ++) {
               infile1 >> dphi[n+N*(m+M*l)][0];
 }
 for(int n=0; n <N; n ++) for (int m=0; m <M; m ++) for (int l=0; l <dlocal_n0; l ++) { 
               dphi[n+N*(m+M*l)][1]=0;
 }
 infile1.close();


 sprintf(filenamec, "%s%d%s%d%s", "con_", 0,"_",dmyid,".txt");
 ifstream infile2(filenamec,ios_base::in);
 for(int n=0; n <N; n ++) for (int m=0; m <M; m ++) for (int l=0; l <dlocal_n0; l ++) {
               infile2 >> dcon[n+N*(m+M*l)][0];
 }
 for(int n=0; n <N; n ++) for (int m=0; m <M; m ++) for (int l=0; l <dlocal_n0; l ++) { 
               dcon[n+N*(m+M*l)][1]=0.0;
 }
 infile2.close();


}

static void atompositionCFGC2(fftw_complex *dphi,int dlocal_n0,int dlocal_0_start, int dmyid,int dnumprocs){
   MPI_Status stat;
   int left=dmyid-1,right=dmyid+1,
	   erro0=0,erro1=0,erro2=0,erro3=0,erro4=0,erro5=0,erro6=0,erro7=0,erro8=0,erro9=0,
	   erro10=0,erro11=0,erro12=0,erro13=0,erro14=0,erro15=0,erro16=0,erro17=0,erro18=0;
   int head=0,tail=dnumprocs-1;
   if (dmyid>0) {left=dmyid-1;}else{left=MPI_PROC_NULL;}
   if (dmyid<dnumprocs-1){right=dmyid+1;}else{right=MPI_PROC_NULL;}

////////find the max/min value
   float phic, phimax=phimm, bufphimax;
   for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {
	   if(dphi[n+N*(m+M*l)][0]>phimax) {phimax=dphi[n+N*(m+M*l)][0];}
   }

   float pii;
   MPI_Allreduce(&phimax,&pii,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);//�����н��̽���һ��ȡ���ֵ��?������Ϊpii
   phic=u0+(pii-u0)*phi_rat; 

   ofstream ophimax("phimax.txt",ios_base::out|ios_base::app);
   ophimax<<kd<<" , "<< phic<<" , "<<  phimax<<" , "<< dmyid<<endl;

   ophimax<<"phimax is done"<<" "<<"myid="<<dmyid<<endl;
 
   float *buf;
   buf = (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+2));
   for(int l=0; l <dlocal_n0+2; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=0;
   }
   for(int l=1; l <dlocal_n0+1; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=dphi[n+N*(m+M*(l-1))][0];//ƽ����l����
   }
   ophimax<<" buf is done "<<endl;
    if(dmyid==head) {MPI_Sendrecv(&buf[N*M],N*M,MPI_FLOAT,tail,40,&buf[0],N*M,MPI_FLOAT,tail,40,MPI_COMM_WORLD,&stat); }
    if(dmyid==tail) {MPI_Sendrecv(&buf[(dlocal_n0)*N*M],N*M,MPI_FLOAT,head,40,&buf[(dlocal_n0+1)*N*M],N*M,MPI_FLOAT,head,40,MPI_COMM_WORLD,&stat);}

    MPI_Sendrecv(&buf[(dlocal_n0)*M*N],N*M,MPI_FLOAT,right,20,&buf[0],N*M,MPI_FLOAT,left,20,MPI_COMM_WORLD,&stat); //->������
    MPI_Sendrecv(&buf[N*M],N*M,MPI_FLOAT,left,20,&buf[N*M*(dlocal_n0+1)],N*M,MPI_FLOAT,right,20,MPI_COMM_WORLD,&stat);      //<-������  

    float * temp;
    temp =  (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
    for(int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       temp[n+N*(m+M*l)]=buf[n+N*(m+M*l)];
    }
   int xtemp=0,ytemp=0,ztemp=0,num_atoms=0;
   int lf=0,lb=0,mf=0,mb=0,nf=0,nb=0,
	   lf2=0,lb2=0,mf2=0,mb2=0,nf2=0,nb2=0;
  // float pc[L*M][3];

   float *pc_x,*pc_y,*pc_z;
   pc_x=(float*) malloc(L*M*100*sizeof(float));
   pc_y=(float*) malloc(L*M*100*sizeof(float));
   pc_z=(float*) malloc(L*M*100*sizeof(float));

   ophimax<<" pc_xyz is initializing "<<endl;   


  for (int l=1; l <dlocal_n0+1; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {

          if(buf[n+N*(m+M*l)]>=phic) {
                            lf=l-1;
							lb=l+1;
							lf2=l-2;
							lb2=l+2;
                            if( (buf[n+N*(m+M*l)] >= buf[n+N*(m+M*lf)]) && (buf[n+N*(m+M*l)] >= buf[n+N*(m+M*lb)]) )// &&
							//	(buf[n+N*(m+M*l)] >= buf[n+N*(m+M*lf2)]) && (buf[n+N*(m+M*l)] >= buf[n+N*(m+M*lb2)]) )
							{
                               xtemp=(dlocal_0_start+l);
                               mf=m-1;							   
                               mb=m+1;
							   mf2=m-2;
							   mb2=m+2;
                               if ((buf[n+N*(m+M*l)] >= buf[n+N*(mf+M*l)]) && (buf[n+N*(m+M*l)] >= buf[n+N*(mb+M*l)]) )// &&
							//	   (buf[n+N*(m+M*l)] >= buf[n+N*(mf2+M*l)]) && (buf[n+N*(m+M*l)] >= buf[n+N*(mb2+M*l)]) )
							   {     
                                  ytemp=m+1;
                                  nf=n-1;								  
                                  nb=n+1;
								  nf2=n-2;
								  nb2=n+2;
                                  if ((buf[n+N*(m+M*l)] >= buf[nf+N*(m+M*l)]) && (buf[n+N*(m+M*l)] >= buf[nb+N*(m+M*l)]) )//&&
							//		  (buf[n+N*(m+M*l)] >= buf[nf2+N*(m+M*l)]) && (buf[n+N*(m+M*l)] >= buf[nb2+N*(m+M*l)]) )
								  { 
                                      ztemp=n+1;								
                                      pc_x[num_atoms]=(ytemp*dx)/(frame_y);
                                      pc_y[num_atoms]=(ztemp*dx)/(frame_z);
                                      pc_z[num_atoms]=(xtemp*dx)/(frame_x);
                                      num_atoms++;
							      
								  }
							   
							   }
					
							}
		
		  }
   }

   ophimax<<" pc_xyz is done "<<endl; 

  if (num_atoms==0) {pc_x[0]=-1; pc_y[0]=-1; pc_z[0]=-1;num_atoms=1;}// add false data if num_atoms is zero
    cout<<num_atoms<<" "<<kd<<" "<<dmyid<<endl;
   int * rcounts, *displs;


   rcounts = (int *) malloc(dnumprocs*sizeof(int));
   displs = (int *) malloc(dnumprocs*sizeof(int));


   int size_pc=num_atoms;

   MPI_Gather(&size_pc,1,MPI_INT,rcounts,1,MPI_INT,0,MPI_COMM_WORLD);

   MPI_Bcast(rcounts,dnumprocs,MPI_INT,0,MPI_COMM_WORLD);
    cout<<"second"<<"myid="<<dmyid<<endl;
    MPI_Barrier(MPI_COMM_WORLD); 
   displs[0]=0;
   for(int i=1; i<dnumprocs; ++i){displs[i]=displs[i-1]+rcounts[i-1];}
   int size_buf=displs[dnumprocs-1]+rcounts[dnumprocs-1];
   
  // float rbuf[size_buf][3];
     float *rbuf_x,*rbuf_y,*rbuf_z;
     rbuf_x=(float*) malloc(size_buf*sizeof(float));
     rbuf_y=(float*) malloc(size_buf*sizeof(float));
     rbuf_z=(float*) malloc(size_buf*sizeof(float));

   MPI_Gatherv(pc_x,size_pc,MPI_FLOAT,rbuf_x,rcounts,displs,MPI_FLOAT,0,MPI_COMM_WORLD);   
   MPI_Gatherv(pc_y,size_pc,MPI_FLOAT,rbuf_y,rcounts,displs,MPI_FLOAT,0,MPI_COMM_WORLD); 
   MPI_Gatherv(pc_z,size_pc,MPI_FLOAT,rbuf_z,rcounts,displs,MPI_FLOAT,0,MPI_COMM_WORLD); 

      ophimax<<" pc_xyz has been gathered "<<endl; 

    cout<<"third"<<"myid="<<dmyid<<endl;
   if (dmyid==0) {
    int num_atoms_total=0;
    for (int i=0; i<size_buf; i++) {if (rbuf_x[i]>=0) {num_atoms_total++; }}
    cout<<"num_atoms_total="<<num_atoms_total<<endl;
    char filename[20];
    sprintf(filename, "%s%d%s", "phiball_", kd,".cfg");
    ofstream outfile(filename,ios_base::out);
    outfile.setf(ios_base::fixed,ios_base::floatfield);
    outfile.precision(4);
    outfile << "Number of particles = " << num_atoms_total << endl;
    outfile << "A = "<< ai <<" Angstrom (basic length-scale)" << endl;   
    outfile << "H0(1,1) = "<<frame_y<<" A "<<endl;
    outfile << "H0(1,2) = "<<0<<" A "<<endl;
    outfile << "H0(1,3) = "<<0<<" A "<<endl;
    outfile << "H0(2,1) = "<<0<<" A "<<endl;
    outfile << "H0(2,2) = "<<frame_z<<" A "<<endl;
    outfile << "H0(2,3) = "<<0<<" A "<<endl;  
    outfile << "H0(3,1) = "<<0<<" A "<<endl;
    outfile << "H0(3,2) = "<<0<<" A "<<endl;
    outfile << "H0(3,3) = "<<frame_x<<" A "<<endl;
    outfile << ".NO_VELOCITY."<<endl;
    outfile << "entry_count = 3"<<endl;
    outfile << "63.546"<<endl;
    outfile << "Cu"<<endl;

    for (int i=0; i<size_buf;i++) {
        if (rbuf_x[i]>=0)   outfile <<rbuf_x[i]<<" "<<rbuf_y[i]<<" "<<rbuf_z[i]<<endl;} 
    outfile.close();
	ophimax<<" atomeye is done "<<endl;
      ophimax.close();

	GR( rbuf_x,  rbuf_y,  rbuf_z, size_buf);
	if(kd>=0 && kd%20==0) Q6_GLO( rbuf_x,  rbuf_y,  rbuf_z, size_buf);
   }
free(temp);
free(pc_x);
free(pc_y);
free(pc_z);
free(rbuf_x);
free(rbuf_y);
free(rbuf_z);
free(buf);
temp=0;
pc_x=0;
pc_y=0;
pc_z=0;
rbuf_x=0;
rbuf_y=0;
rbuf_z=0;
buf=0;
}


static void GR(float *drbuf_x, float *drbuf_y, float *drbuf_z, int dsize_buf){
	const float delr=a_lattice/25.;
	const int maxbin=static_cast<int>(L*dx/delr), nstep=1, dnum=dsize_buf, CNmax=24;
	unsigned bin, *hist=0; //*Nb=0
	float  drx,dry,drz,rjx,rjy,rjz,rr, r3[9]={0}, *gr=0;
//	Nb= (unsigned *) fftw_malloc(sizeof(unsigned) * dnum);
	hist=(unsigned *) fftw_malloc(sizeof(unsigned) * maxbin);
	gr=(float *) fftw_malloc(sizeof(float) * maxbin);

	int num_t=0;
	for (int tt=0; tt<dsize_buf; tt++) {if (drbuf_x[tt]>=0) {num_t++; }}
    cout<<"dsize_buf="<<dsize_buf<<" , num_t="<<num_t<<endl;
    
    for (int ii=0; ii<maxbin; ii++){
	         gr[ii]=0; hist[ii]=0;	}

      cout<<"initial hist, gr"<<endl;

	for(int i=0; i<dsize_buf-1; i++){
	  for(int j=i+1; j<dsize_buf; j++){
		r3[0]=(drbuf_x[i]-drbuf_x[j])*frame_x;
		r3[1]=(drbuf_y[i]-drbuf_y[j])*frame_y;
		r3[2]=(drbuf_z[i]-drbuf_z[j])*frame_z;
		r3[3]=r3[0]+frame_x;
		r3[4]=r3[1]+frame_y;
		r3[5]=r3[2]+frame_z;
		r3[6]=r3[0]-frame_x;
		r3[7]=r3[1]-frame_y;
		r3[8]=r3[2]-frame_z;
		drx=min3(r3[0], r3[3], r3[6]); 
		dry=min3(r3[1], r3[4], r3[7]); 
		drz=min3(r3[2], r3[5], r3[8]); 
		rr=sqrt(drx*drx+dry*dry+drz*drz)+1e-6;
	
		bin=static_cast<int>(rr/delr);
		if(rr>a_lattice/sqrt(2.)*0.5 && bin< maxbin) { hist[bin]=hist[bin]+2; }
	  }
	}
	char histname[20];
    sprintf(histname, "%s%d%s", "hist_", kd,".txt");
    ofstream outhist(histname,ios_base::out);
    outhist.setf(ios_base::fixed,ios_base::floatfield);
    outhist.precision(4);
	for (int ii=0; ii<maxbin; ii++){ 
		outhist << ii << " " << hist[ii]<< " " <<maxbin<<endl;
	}
	outhist.close();
	
//	nstep=1;
	float rho=float(num_t)/float(L*M*N*dx*dx*dx);
	float preco=4.*pi*rho/3.;
    cout<<rho<<" "<<preco<<endl;

    char filename1[20];
    sprintf(filename1, "%s%d%s", "Grq6_", kd,".txt");
    ofstream outgr(filename1,ios_base::out);
    outgr.setf(ios_base::fixed,ios_base::floatfield);
    outgr.precision(4);
	for (int jj=0; jj<maxbin; jj++){ 
	float rlower=float(jj)*delr;
	float rupper=rlower+delr;
	float nideal=preco*(rupper*rupper*rupper-rlower*rlower*rlower);
	   gr[jj]=float(hist[jj])/float(nstep)/float(num_t)/nideal;
	   outgr << jj*delr << " " << gr[jj] << " "<< hist[jj] <<" "<< nideal << endl;}
	outgr.close();

free(gr);
free(hist);
hist=0;
gr=0;
}

static void Q6_GLO(float *drbuf_x, float *drbuf_y, float *drbuf_z, int dsize_buf){
	const int nstep=1, dnum=dsize_buf, CNmax=24;
    unsigned *Nb=0; 
	float drx,dry,drz,rjx,rjy,rjz,rr, r3[9]={0}; 
	Nb= (unsigned *) fftw_malloc(sizeof(unsigned) * dnum);

	float *CoX=0, *CoY=0, *CoZ=0, *rij=0, *fai=0,*theta=0;
	rij=(float *) fftw_malloc(sizeof(float) * dnum*CNmax);
	CoX=(float *) fftw_malloc(sizeof(float) * dnum*CNmax);
	CoY=(float *) fftw_malloc(sizeof(float) * dnum*CNmax);
	CoZ=(float *) fftw_malloc(sizeof(float) * dnum*CNmax);
	theta=(float *) fftw_malloc(sizeof(float) * dnum*CNmax);
	fai=(float *) fftw_malloc(sizeof(float) * dnum*CNmax);

	int num_t=0;
	for (int tt=0; tt<dsize_buf; tt++) {if (drbuf_x[tt]>=0) {num_t++; }}
    cout<<"dsize_buf="<<dsize_buf<<" , num_t="<<num_t<<endl;
    
	for (int ii=0; ii<dnum; ii++){
	         Nb[ii]=0;   }
	for (int ii=0; ii<dnum*CNmax; ii++){
		 theta[ii]=0.; fai[ii]=0.;}
		 
      cout<<"initial Nb"<<endl;

	for(int i=0; i<dsize_buf; i++){
	  for(int j=0; j<dsize_buf; j++){
	   if(i!=j)
	   {
		r3[0]=(drbuf_x[i]-drbuf_x[j])*frame_x;
		r3[1]=(drbuf_y[i]-drbuf_y[j])*frame_y;
		r3[2]=(drbuf_z[i]-drbuf_z[j])*frame_z;
		r3[3]=r3[0]+frame_x;
		r3[4]=r3[1]+frame_y;
		r3[5]=r3[2]+frame_z;
		r3[6]=r3[0]-frame_x;
		r3[7]=r3[1]-frame_y;
		r3[8]=r3[2]-frame_z;
		drx=min3(r3[0], r3[3], r3[6]); 
		dry=min3(r3[1], r3[4], r3[7]); 
		drz=min3(r3[2], r3[5], r3[8]); 
		rr=sqrt(drx*drx+dry*dry+drz*drz)+1e-6; 		
		
//		if(rr<=a_lattice*sqrt(2.)/2.*1.2+1e-3 && rr>=a_lattice*sqrt(2.)/2.*0.9-1e-3)
		if(rr<=a_lattice*sqrt(3.)/2.*1.2+1e-3 && rr>=a_lattice*sqrt(3.)/2.*0.9-1e-3)
		{ 
          Nb[i]=Nb[i]+1;  
		  if(Nb[i]<=CNmax) { 
			rij[Nb[i]-1+CNmax*i]=rr;
            for(int ii=0; ii<7; ii=ii+3){
			   if(fabs(r3[ii])-drx<1e-5) rjx=drbuf_x[i]*frame_x-r3[ii];}		  
            for(int ii=1; ii<8; ii=ii+3){
			   if(fabs(r3[ii])-dry<1e-5) rjy=drbuf_y[i]*frame_y-r3[ii];}		  
            for(int ii=2; ii<9; ii=ii+3){
			   if(fabs(r3[ii])-drz<1e-5) rjz=drbuf_z[i]*frame_z-r3[ii];}		  
            CoX[Nb[i]-1+CNmax*i]=rjx;
	    	CoY[Nb[i]-1+CNmax*i]=rjy;
			CoZ[Nb[i]-1+CNmax*i]=rjz;
		  }
		}
	   }
	  }
	}


/*	char filename[20];
    sprintf(filename, "%s%d%s", "Nb_", kd,".txt");
    ofstream outNb(filename,ios_base::out);
    outNb.setf(ios_base::fixed,ios_base::floatfield);
    outNb.precision(4);
	for (int ii=0; ii<dnum; ii++){ 
		outNb << ii << " " << Nb[ii] <<endl;
	}
	outNb.close();
*/   

	cout<<kd<<" initial Cox, rij, ready for theta & fai "<<endl;
    for(int i=0; i<dnum; i++){
	 if(Nb[i]>0  && Nb[i]<=CNmax){
		for (int j=0; j<Nb[i]; j++){		
			theta[j+CNmax*i]=acos( (CoZ[j+CNmax*i]-drbuf_z[i]*frame_z)/rij[j+CNmax*i] );		
			if(CoY[j+CNmax*i]-drbuf_y[i]*frame_y ==0.) {fai[j+CNmax*i]=0.;} //CoX[j+CNmax*i]-drbuf_x[i]*frame_x<1e-5
			else 
		    fai[j+CNmax*i]=atan( (CoY[j+CNmax*i]-drbuf_y[i]*frame_y)/(CoX[j+CNmax*i]-drbuf_x[i]*frame_x) );
//		outrxyz<<i<<" "<<j<<" "<<(CoY[j+CNmax*i]-drbuf_y[i]*frame_y)<<" "<<(CoX[j+CNmax*i]-drbuf_x[i]*frame_x)<<" "<<fai[j+CNmax*i]<<endl;
		}
	 }
	}
//	outrxyz.close();

/*	char angname[20];
    sprintf(angname, "%s%d%s", "ang_", kd,".txt");
    ofstream outang(angname,ios_base::out);
    outang.setf(ios_base::fixed,ios_base::floatfield);
    outang.precision(6);
	for(int i=0; i<dnum; i++){
	  for(int j=0; j<Nb[i]; j++){
	    outang<<i<<" "<<j<<" "<<rij[j+CNmax*i]<<" "<<theta[j+CNmax*i]<<" "<<fai[j+CNmax*i]<<endl;
	  }
	}
    outang.close();
*/
  	char Qlname[20];
    sprintf(Qlname, "%s%d%s", "Ql_", kd,".txt");
    ofstream outQl(Qlname,ios_base::out);
    outQl.setf(ios_base::fixed,ios_base::floatfield);
    outQl.precision(6);

    const int l=6;
    float Qlm_m=0., *qlmb_m=0, *Ql_L=0;
	qlmb_m=(float *) fftw_malloc(sizeof(float) * dsize_buf);
	Ql_L=(float *) fftw_malloc(sizeof(float) * dsize_buf);
    for(int i=0; i<dnum; i++)  { qlmb_m[i]=0.;  Ql_L[i]=0.; }

	fftw_complex *qlm=0, *qlmb=0;
	qlm= (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * dnum);
	qlmb= (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * dnum);    

   for(int m=-l; m<=l; m++)
	{	  
	  for (int ii=0; ii<dnum; ii++)
	  {
	     qlm[ii][0]=0.; qlm[ii][1]=0.;
		qlmb[ii][0]=0.; qlmb[ii][1]=0.;
	   }	  
	 for(int i=0; i<dnum; i++){
	      if(Nb[i]>0 && Nb[i]<=CNmax){
		for (int j=0; j<Nb[i]; j++){
	       qlm[i][0]+=Ylmtf_r(l, m, theta[j+CNmax*i] ,fai[j+CNmax*i] );
		   qlm[i][1]+=Ylmtf_c(l, m, theta[j+CNmax*i] ,fai[j+CNmax*i] );
		   if( j==int(Nb[i])-1 ) {
	         qlmb[i][0]=qlm[i][0]/float(Nb[i]);
	         qlmb[i][1]=qlm[i][1]/float(Nb[i]);  }
		}
		  }
	 }
    float SumNb=0., Qlm_r=0., Qlm_c=0.;
	for(int i=0; i<dnum; i++){
		qlmb_m[i]+=qlmb[i][0]*qlmb[i][0]+qlmb[i][1]*qlmb[i][1];
	       Qlm_r+=Nb[i]*qlmb[i][0];
		   Qlm_c+=Nb[i]*qlmb[i][1];
		   SumNb+=Nb[i];
	}
	Qlm_r/=SumNb;
	Qlm_c/=SumNb;	
    Qlm_m+=(Qlm_r*Qlm_r+Qlm_c*Qlm_c);
//	outQl << m<<" "<<Qlm_r << " " << Qlm_c <<" "<<Qlm_m<<endl;
  }
	const float delq=0.05;
	const int maxbin=30;
    unsigned bin=0;
	unsigned long *histq=0;
	histq=(unsigned long *) fftw_malloc(sizeof(unsigned long) * maxbin);
	for(int i=0; i<maxbin; i++){ histq[i]=0; }

    float xishu=4*pi/float(2*l+1);
    float Ql=sqrt(xishu*Qlm_m); //GLOBAL
//  outQl <<xishu << " " << Ql <<endl;
    for(int i=0; i<dnum; i++)
  {  
	    Ql_L[i]=sqrt(xishu*qlmb_m[i]); 
		float Ql=sqrt(xishu*Qlm_m);
        outQl << i <<" "<<Nb[i]<< " " << Ql_L[i] <<" "<<Ql<<endl; //LOCAL;

		bin=int(Ql_L[i]/delq);
		if(bin<maxbin) histq[bin]=histq[bin]+1;
  }
	outQl.close();
    cout<<kd<<", Ql is done"<<endl;
  	char histqname[20];
    sprintf(histqname, "%s%d%s", "HistQ_", kd,".txt");
    ofstream outHistQ(histqname,ios_base::out);
    outHistQ.setf(ios_base::fixed,ios_base::floatfield);
    outHistQ.precision(4);
	for (int i=0; i<maxbin; i++)
       outHistQ<<(i+1)*delq<<" "<<histq[i]<<endl;
    outHistQ.close();
	
free(rij);
free(CoX);
free(CoY);
free(CoZ);
free(Nb);
free(theta);
free(fai);
free(qlm);
free(qlmb);
free(qlmb_m);
free(Ql_L);
//free(gr);
//free(hist);
Nb=0;
//hist=0;
CoX=0;
CoY=0;
CoZ=0;
//gr=0;
theta=0;
rij=0;
qlm=0;
qlmb=0;
qlmb_m=0;
Ql_L=0;
}

// facto provided by pfc_common.h

// Ylmtf_r provided by pfc_common.h

// Ylmtf_c provided by pfc_common.h

// lg_poly provided by pfc_common.h

// plgndr provided by pfc_common.h

//void sVTKwriteBianry(fftw_complex *dphi,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename){
//   char filename[20];
//   sprintf(filename, "%s%d%s%d%s%d%s", dfilename,dmyid,"_",int(100*sig),"_",int(100*con0),".vti");
//
//    int startn,endn;
//   if (dmyid==0){startn=dlocal_0_start; endn=dlocal_0_start+dlocal_n0-1;}
//   else{startn=dlocal_0_start-1; endn=dlocal_0_start+dlocal_n0-1;}
//
//   
//
//    MPI_Status stat;
//   int left=dmyid-1,right=dmyid+1;
//   float *buf, *temp;
//   buf = (float *) fftw_malloc(sizeof(float) * M*N*(dlocal_n0+1));
//   temp= (float *) fftw_malloc(sizeof(float) * M*N*(dlocal_n0));
//   if (dmyid>0) {left=dmyid-1;}else{left=MPI_PROC_NULL;}
//   if (dmyid<dnumprocs-1){right=dmyid+1;}else{right=MPI_PROC_NULL;}
//   for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {
//               temp[n+N*(m+M*l)]=dphi[n+N*(m+M*l)][0];
//        }
//   MPI_Sendrecv(&temp[N*M*(dlocal_n0-1)],N*M,MPI_FLOAT,right,10,&buf[0],N*M,MPI_FLOAT,left,10,MPI_COMM_WORLD,&stat);
//
//
//
	float treal_val;
//   int    nvar=1; /* Plot only one variable set */
//   int    bytes[100], off[100];
//   off[0]=0;
//   for (int i=0; i<nvar; i++){
//    if (dmyid==0){bytes[i]=(N*M*dlocal_n0)*sizeof(float);} else {bytes[i]=(N*M*(dlocal_n0+1))*sizeof(float);}
//    if (i<nvar-1) off[i+1]=off[i]+sizeof(int)+bytes[i];
//    bytes[i]=bytes[i]+sizeof(int);
//   }
//
//
//  FILE * outf_tfield;
//  outf_tfield=fopen(filename,"w");
//  fprintf(outf_tfield,"<?xml version=\"1.0\"?>\n");
//  fprintf(outf_tfield,"<VTKFile type=\"ImageData\" version=\"0.1\"  byte_order=\"LittleEndian\">\n");
//  fprintf(outf_tfield,"<ImageData WholeExtent=\"%d %d %d %d %d %d\"   Origin=\"%d %d %d\" Spacing=\"%d %d %d\">\n",startn,endn,0,M-1,0,N-1
//   ,0,0,0,1,1,1);
//  fprintf(outf_tfield,"<Piece Extent=\"%d %d %d %d %d %d\">\n",startn,endn,0,M-1,0,N-1); 
//  fprintf(outf_tfield,"<PointData Scalars=\"phi\">\n");
//  fprintf(outf_tfield,"<DataArray type=\"Float32\" Name=\"phi\"   format=\"appended\" offset=\"%d\" />\n",off[0]); 
//  fprintf(outf_tfield,"</PointData>\n");
//  fprintf(outf_tfield,"<CellData>\n");
//  fprintf(outf_tfield,"</CellData>\n");
//  fprintf(outf_tfield,"</Piece>\n");
//  fprintf(outf_tfield,"</ImageData>\n");
//  fprintf(outf_tfield,"<AppendedData encoding=\"raw\">\n");
//  fprintf(outf_tfield,"_");
//
//  /* Write arrival time field */
//  fwrite(&bytes[0],sizeof(int),1,outf_tfield);
//  if (dmyid==0){
//                 for(int n=0; n<N; n++){
//                    for(int m=0; m<M; m++){
//                        for(int l=0; l<dlocal_n0; l++){
//                            treal_val = dphi[n+N*(m+M*l)][0];
//                            fwrite(&treal_val,sizeof(float),1,outf_tfield);
//                         }
//                     }
//                 }
//  } else {
//                 for(int l=1; l <dlocal_n0+1; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
//                   buf[n+N*(m+M*l)]=temp[n+N*(m+M*(l-1))];
//                 }
//                 for(int n=0; n<N; n++){
//                    for(int m=0; m<M; m++){
//                        for(int l=0; l<dlocal_n0+1; l++){
//                            treal_val = buf[n+N*(m+M*l)];
//                            fwrite(&treal_val,sizeof(float),1,outf_tfield);
//                         }
//                     }
//                 }
//  }
//  fprintf(outf_tfield,"\n");
//  fprintf(outf_tfield,"</AppendedData>\n");
//  fprintf(outf_tfield,"</VTKFile>\n");
//  fclose(outf_tfield);
//  
//  free(buf);
//  free(temp);
//
//}
//
//
//
//void PVTKWRITE(int dlocal_n0,int dnumprocs,char *dfilename){
//   char filename[20];
//   char filenames[20];   
//   sprintf(filename, "%s%d%s%d%s", dfilename, int(100*sig),"_",int(100*con0),".pvti");
//   ofstream outfile(filename,ios_base::out);
//   outfile.setf(ios_base::fixed,ios_base::floatfield);
//   outfile.precision(4);
//   outfile << "<?xml version=\"1.0\"?>" << endl;
//   outfile << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
//   outfile << "  <PImageData WholeExtent=\"0 "<<L/RD-1<<" 0 "<<M/RD-1<<" 0 "<<N/RD-1<<"\" GhostLevel=\"#\" Origin=\"0  0  0\"  Spacing=\"1  1  1\">"<<endl;   
//   outfile << "    <PPointData>"<< endl;
//   outfile << "      <PDataArray type=\"Float32\" Name=\"Phi\"/>" << endl;
////   outfile << "      <PDataArray type=\"Float32\" Name=\"enei\"/>" << endl;
//   outfile << "    </PPointData>"<<endl;
//   
//   sprintf(filenames, "%s%d%s%d%s%d%s", dfilename,0,"_",int(100*sig),"_",int(100*con0),".vti");
//   outfile << "    <Piece Extent=\"0 "<<(dlocal_n0/RD)-1<<" 0 "<<M/RD-1<<" 0 "<<N/RD-1<<"\" Source=\""<<filenames<<"\"/>"<<endl;
//   for (int i=1; i<dnumprocs; i ++) {
//   sprintf(filenames, "%s%d%s%d%s%d%s", dfilename,i,"_",int(100*sig),"_",int(100*con0),".vti");
//   outfile << "    <Piece Extent=\""<<i*(dlocal_n0/RD)-1<<"  "<<(i+1)*(dlocal_n0/RD)-1<<" 0 "<<M/RD-1<<" 0 "<<N/RD-1<<"\" Source=\""<<filenames<<"\"/>"<<endl;
//   }
//
//   outfile << "  </PImageData>"<<endl;
//   outfile << "</VTKFile>"<<endl;
//   outfile.close();
//
//}

static void PVTKWRITE(int dlocal_n0, int dnumprocs, int dstartz, int dendz, char *dfilename) {
	char filename[20];
	char filenames[20];
	sprintf(filename, "%s%d%s", dfilename, kd, ".pvti");
	ofstream outfile(filename, ios_base::out);
	outfile.setf(ios_base::fixed, ios_base::floatfield);
	outfile.precision(4);
	outfile << "<?xml version=\"1.0\"?>" << endl;
	outfile << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
	outfile << "  <PImageData WholeExtent=\"0 " << L / RD - 1 << " 0 " << M / RD - 1 << " " << dstartz << " " << dendz - 1 << "\" GhostLevel=\"#\" Origin=\"0  0  0\"  Spacing=\"1  1  1\">" << endl;
	outfile << "    <PPointData>" << endl;
	outfile << "      <PDataArray type=\"Float32\" Name=\"phi\"/>" << endl;
	//   outfile << "      <PDataArray type=\"Float32\" Name=\"enei\"/>" << endl;
	outfile << "    </PPointData>" << endl;

	sprintf(filenames, "%s%d%s%d%s", dfilename, kd, "_", 0, ".vti");
	outfile << "    <Piece Extent=\"0 " << (dlocal_n0 / RD) - 1 << " 0 " << M / RD - 1 << " " << dstartz << " " << dendz - 1 << "\" Source=\"" << filenames << "\"/>" << endl;
	for (int i = 1; i<dnumprocs; i++) {
		sprintf(filenames, "%s%d%s%d%s", dfilename, kd, "_", i, ".vti");
		outfile << "    <Piece Extent=\"" << i*(dlocal_n0 / RD) - 1 << "  " << (i + 1)*(dlocal_n0 / RD) - 1 << " 0 " << M / RD - 1 << " " << dstartz << " " << dendz - 1 << "\" Source=\"" << filenames << "\"/>" << endl;
	}

	outfile << "  </PImageData>" << endl;
	outfile << "</VTKFile>" << endl;
	outfile.close();

}

static void sVTKwriteBianry(fftw_complex *dphi, int dlocal_n0, int dlocal_0_start, int dmyid, int dnumprocs) {
	char filename[20];
	sprintf(filename, "%s%d%s%d%s", "phi_", kd, "_", dmyid, ".vti");


	int startn, endn;
	if (dmyid == 0) { startn = dlocal_0_start; endn = dlocal_0_start + dlocal_n0 - 1; }
	else { startn = dlocal_0_start - 1; endn = dlocal_0_start + dlocal_n0 - 1; }


	MPI_Status stat;
	int left = dmyid - 1, right = dmyid + 1;
	float *buf, *temp;
	buf = (float *)fftw_malloc(sizeof(float) * M*N*(dlocal_n0 + 1));
	temp = (float *)fftw_malloc(sizeof(float) * M*N*(dlocal_n0));
	if (dmyid>0) { left = dmyid - 1; }
	else { left = MPI_PROC_NULL; }
	if (dmyid<dnumprocs - 1) { right = dmyid + 1; }
	else { right = MPI_PROC_NULL; }
	for (int l = 0; l <dlocal_n0; l++) for (int m = 0; m <M; m++) for (int n = 0; n <N; n++) {
		temp[n + N*(m + M*l)] = dphi[n + N*(m + M*l)][0];
	}
	MPI_Sendrecv(&temp[N*M*(dlocal_n0 - 1)], N*M, MPI_FLOAT, right, 10, &buf[0], N*M, MPI_FLOAT, left, 10, MPI_COMM_WORLD, &stat);
	float treal_val;
	int    nvar = 1; /* Plot only one variable set */
	int    bytes[100], off[100];
	off[0] = 0;
	for (int i = 0; i<nvar; i++) {
		if (dmyid == 0) { bytes[i] = (N*M*dlocal_n0) * sizeof(float); }
		else { bytes[i] = (N*M*(dlocal_n0 + 1)) * sizeof(float); }
		if (i<nvar - 1) off[i + 1] = off[i] + sizeof(int) + bytes[i];
		bytes[i] = bytes[i] + sizeof(int);
	}


	FILE * outf_tfield;
	outf_tfield = fopen(filename, "w");
	fprintf(outf_tfield, "<?xml version=\"1.0\"?>\n");
	fprintf(outf_tfield, "<VTKFile type=\"ImageData\" version=\"0.1\"  byte_order=\"LittleEndian\">\n");
	fprintf(outf_tfield, "<ImageData WholeExtent=\"%d %d %d %d %d %d\"   Origin=\"%d %d %d\" Spacing=\"%d %d %d\">\n", startn, endn, 0, M - 1, 0, N - 1
		, 0, 0, 0, 1, 1, 1);
	fprintf(outf_tfield, "<Piece Extent=\"%d %d %d %d %d %d\">\n", startn, endn, 0, M - 1, 0, N - 1);
	fprintf(outf_tfield, "<PointData Scalars=\"phi\">\n");
	fprintf(outf_tfield, "<DataArray type=\"Float32\" Name=\"phi\"   format=\"appended\" offset=\"%d\" />\n", off[0]);
	fprintf(outf_tfield, "</PointData>\n");
	fprintf(outf_tfield, "<CellData>\n");
	fprintf(outf_tfield, "</CellData>\n");
	fprintf(outf_tfield, "</Piece>\n");
	fprintf(outf_tfield, "</ImageData>\n");
	fprintf(outf_tfield, "<AppendedData encoding=\"raw\">\n");
	fprintf(outf_tfield, "_");

	/* Write arrival time field */
	fwrite(&bytes[0], sizeof(int), 1, outf_tfield);
	if (dmyid == 0) {
		for (int n = 0; n<N; n++) {
			for (int m = 0; m<M; m++) {
				for (int l = 0; l<dlocal_n0; l++) {
					treal_val = dphi[n + N*(m + M*l)][0];
					fwrite(&treal_val, sizeof(float), 1, outf_tfield);
				}
			}
		}
	}
	else {
		for (int l = 1; l <dlocal_n0 + 1; l++) for (int m = 0; m <M; m++) for (int n = 0; n <N; n++) {
			buf[n + N*(m + M*l)] = temp[n + N*(m + M*(l - 1))];
		}
		for (int n = 0; n<N; n++) {
			for (int m = 0; m<M; m++) {
				for (int l = 0; l<dlocal_n0 + 1; l++) {
					treal_val = buf[n + N*(m + M*l)];
					fwrite(&treal_val, sizeof(float), 1, outf_tfield);
				}
			}
		}
	}
	fprintf(outf_tfield, "\n");
	fprintf(outf_tfield, "</AppendedData>\n");
	fprintf(outf_tfield, "</VTKFile>\n");
	fclose(outf_tfield);

	free(buf);
	free(temp);

}

static void PVTKWRITEC(int dlocal_n0, int dnumprocs, int dstartz, int dendz, char *dfilename) {
	char filename[20];
	char filenames[20];
	sprintf(filename, "%s%d%s", dfilename, kd, ".pvti");
	ofstream outfile(filename, ios_base::out);
	outfile.setf(ios_base::fixed, ios_base::floatfield);
	outfile.precision(4);
	outfile << "<?xml version=\"1.0\"?>" << endl;
	outfile << "<VTKFile type=\"PImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
	outfile << "  <PImageData WholeExtent=\"0 " << L / RD - 1 << " 0 " << M / RD - 1 << " " << dstartz << " " << dendz - 1 << "\" GhostLevel=\"#\" Origin=\"0  0  0\"  Spacing=\"1  1  1\">" << endl;
	outfile << "    <PPointData>" << endl;
	outfile << "      <PDataArray type=\"Float32\" Name=\"con\"/>" << endl;
	//   outfile << "      <PDataArray type=\"Float32\" Name=\"enei\"/>" << endl;
	outfile << "    </PPointData>" << endl;

	sprintf(filenames, "%s%d%s%d%s", dfilename, kd, "_", 0, ".vti");
	outfile << "    <Piece Extent=\"0 " << (dlocal_n0 / RD) - 1 << " 0 " << M / RD - 1 << " " << dstartz << " " << dendz - 1 << "\" Source=\"" << filenames << "\"/>" << endl;
	for (int i = 1; i<dnumprocs; i++) {
		sprintf(filenames, "%s%d%s%d%s", dfilename, kd, "_", i, ".vti");
		outfile << "    <Piece Extent=\"" << i*(dlocal_n0 / RD) - 1 << "  " << (i + 1)*(dlocal_n0 / RD) - 1 << " 0 " << M / RD - 1 << " " << dstartz << " " << dendz - 1 << "\" Source=\"" << filenames << "\"/>" << endl;
	}

	outfile << "  </PImageData>" << endl;
	outfile << "</VTKFile>" << endl;
	outfile.close();

}

static void sVTKwriteBianryc(fftw_complex *dcon, int dlocal_n0, int dlocal_0_start, int dmyid, int dnumprocs) {
	char filename[20];
	sprintf(filename, "%s%d%s%d%s", "con_", kd, "_", dmyid, ".vti");


	int startn, endn;
	if (dmyid == 0) { startn = dlocal_0_start; endn = dlocal_0_start + dlocal_n0 - 1; }
	else { startn = dlocal_0_start - 1; endn = dlocal_0_start + dlocal_n0 - 1; }


	MPI_Status stat;
	int left = dmyid - 1, right = dmyid + 1;
	float *buf, *temp;
	buf = (float *)fftw_malloc(sizeof(float) * M*N*(dlocal_n0 + 1));
	temp = (float *)fftw_malloc(sizeof(float) * M*N*(dlocal_n0));
	if (dmyid>0) { left = dmyid - 1; }
	else { left = MPI_PROC_NULL; }
	if (dmyid<dnumprocs - 1) { right = dmyid + 1; }
	else { right = MPI_PROC_NULL; }
	for (int l = 0; l <dlocal_n0; l++) for (int m = 0; m <M; m++) for (int n = 0; n <N; n++) {
		temp[n + N*(m + M*l)] = dcon[n + N*(m + M*l)][0];
	}
	MPI_Sendrecv(&temp[N*M*(dlocal_n0 - 1)], N*M, MPI_FLOAT, right, 10, &buf[0], N*M, MPI_FLOAT, left, 10, MPI_COMM_WORLD, &stat);
	float treal_val;
	int    nvar = 1; /* Plot only one variable set */
	int    bytes[100], off[100];
	off[0] = 0;
	for (int i = 0; i<nvar; i++) {
		if (dmyid == 0) { bytes[i] = (N*M*dlocal_n0) * sizeof(float); }
		else { bytes[i] = (N*M*(dlocal_n0 + 1)) * sizeof(float); }
		if (i<nvar - 1) off[i + 1] = off[i] + sizeof(int) + bytes[i];
		bytes[i] = bytes[i] + sizeof(int);
	}


	FILE * outf_tfield;
	outf_tfield = fopen(filename, "w");
	fprintf(outf_tfield, "<?xml version=\"1.0\"?>\n");
	fprintf(outf_tfield, "<VTKFile type=\"ImageData\" version=\"0.1\"  byte_order=\"LittleEndian\">\n");
	fprintf(outf_tfield, "<ImageData WholeExtent=\"%d %d %d %d %d %d\"   Origin=\"%d %d %d\" Spacing=\"%d %d %d\">\n", startn, endn, 0, M - 1, 0, N - 1
		, 0, 0, 0, 1, 1, 1);
	fprintf(outf_tfield, "<Piece Extent=\"%d %d %d %d %d %d\">\n", startn, endn, 0, M - 1, 0, N - 1);
	fprintf(outf_tfield, "<PointData Scalars=\"con\">\n");
	fprintf(outf_tfield, "<DataArray type=\"Float32\" Name=\"con\"   format=\"appended\" offset=\"%d\" />\n", off[0]);
	fprintf(outf_tfield, "</PointData>\n");
	fprintf(outf_tfield, "<CellData>\n");
	fprintf(outf_tfield, "</CellData>\n");
	fprintf(outf_tfield, "</Piece>\n");
	fprintf(outf_tfield, "</ImageData>\n");
	fprintf(outf_tfield, "<AppendedData encoding=\"raw\">\n");
	fprintf(outf_tfield, "_");

	/* Write arrival time field */
	fwrite(&bytes[0], sizeof(int), 1, outf_tfield);
	if (dmyid == 0) {
		for (int n = 0; n<N; n++) {
			for (int m = 0; m<M; m++) {
				for (int l = 0; l<dlocal_n0; l++) {
					treal_val = dcon[n + N*(m + M*l)][0];
					fwrite(&treal_val, sizeof(float), 1, outf_tfield);
				}
			}
		}
	}
	else {
		for (int l = 1; l <dlocal_n0 + 1; l++) for (int m = 0; m <M; m++) for (int n = 0; n <N; n++) {
			buf[n + N*(m + M*l)] = temp[n + N*(m + M*(l - 1))];
		}
		for (int n = 0; n<N; n++) {
			for (int m = 0; m<M; m++) {
				for (int l = 0; l<dlocal_n0 + 1; l++) {
					treal_val = buf[n + N*(m + M*l)];
					fwrite(&treal_val, sizeof(float), 1, outf_tfield);
				}
			}
		}
	}
	fprintf(outf_tfield, "\n");
	fprintf(outf_tfield, "</AppendedData>\n");
	fprintf(outf_tfield, "</VTKFile>\n");
	fclose(outf_tfield);

	free(buf);
	free(temp);

}

static void OUTPUTPHI(fftw_complex *dphi,fftw_complex *dcon,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start){
 char filename[20];
 sprintf(filename, "%s%d%s%d%s", "phi_", kd,"_",dmyid,".txt");
 ofstream outfile(filename,ios_base::out);
 outfile.setf(ios_base::fixed,ios_base::floatfield);
 outfile.precision(4);
 for(int n=0; n <N; n ++) for (int m=0; m <M; m ++) for (int l=0; l <dlocal_n0; l ++) {
               outfile<<dphi[n+N*(m+M*l)][0]<<endl;
 }
 outfile.close();
}


static void SLICE(fftw_complex *dphi,fftw_complex *dcon,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start){
  char filename1[20],filename2[20],filename3[20],filename4[20];
  sprintf(filename1,"%s%d%s", "slice(100)_", kd,".dat");
  sprintf(filename2,"%s%d%s", "slice(010)_", kd,".dat");

  sprintf(filename3,"%s%d%s", "phi_001", kd,".dat");
  sprintf(filename4,"%s%d%s", "con_001", kd,".dat");
 
/*
// out 100 plane
 int proc_mid=(dnumprocs+1)/2;
  if (dmyid==proc_mid) {     
    ofstream outf1(filename1,ios_base::out);
    outf1.setf(ios_base::fixed,ios_base::floatfield);
    outf1.precision(4);
    int lx=0;    
	outf1<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outf1<<"zone t=\"big zone\""<<",i="<<M<<",j="<<N<<",f=point"<<endl;//dat
    for(int n=0; n <N; n ++)for(int m=0; m <M; m ++)  { outf1<<dphi[n+N*(m+M*lx)][0]<<endl; }
    outf1.close(); 
  }
*/
 /*
  // output 010 plane   
  int my;
  float *temp, *rbuf;
  if ( dmyid == 0) {rbuf= (float *) fftw_malloc(sizeof(float)*L*N);} 
  temp= (float *) fftw_malloc(sizeof(float) * N*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
     my=M/2;//extract m=M/2 
     for (int n=0; n <N; n ++) {temp[n+N*l]=dphi[n+N*(my+M*l)][0];}
  }
  MPI_Gather(temp,N*dlocal_n0,MPI_FLOAT,rbuf,N*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
    ofstream outf3(filename2,ios_base::out);
    outf3.setf(ios_base::fixed,ios_base::floatfield);
    outf3.precision(4);
	outf3<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outf3<<"zone t=\"big zone\""<<",i="<<L<<",j="<<N<<",f=point"<<endl;//dat
    for(int n=0; n <N; n ++)for(int l=0; l <L; l ++) {outf3<<rbuf[n+N*l]<<endl;}    
    outf3.close();
free(temp);
free(rbuf);   
  } 
*/

// output 001 plane 

  float *temp1, *rbuf1;
  if ( dmyid == 0) {rbuf1= (float *) fftw_malloc(sizeof(float)*L*M);} 
  temp1= (float *) fftw_malloc(sizeof(float) * M*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
    int nz=N/2;//
     for (int m=0; m <M; m ++) {temp1[m+M*l]=dphi[nz+N*(m+M*l)][0];}
  }
  MPI_Gather(temp1,M*dlocal_n0,MPI_FLOAT,rbuf1,M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
    ofstream outf2(filename3,ios_base::out);
    outf2.setf(ios_base::fixed,ios_base::floatfield);
    outf2.precision(4);
	outf2<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outf2<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",f=point"<<endl;//dat
    for(int m=0; m <M; m ++)for(int l=0; l <L; l ++) {outf2<<rbuf1[m+M*l]<<endl;}
    outf2.close(); 
free(temp1);
free(rbuf1);   
  } 


  float *temp2, *rbuf2;
  if ( dmyid == 0) {rbuf2= (float *) fftw_malloc(sizeof(float)*L*M);} 
  temp2= (float *) fftw_malloc(sizeof(float) * M*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
    int nz=N/2;//
     for (int m=0; m <M; m ++) {temp2[m+M*l]=dcon[nz+N*(m+M*l)][0];}
  }
  MPI_Gather(temp2,M*dlocal_n0,MPI_FLOAT,rbuf2,M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
    ofstream outf2(filename4,ios_base::out);
    outf2.setf(ios_base::fixed,ios_base::floatfield);
    outf2.precision(4);
	outf2<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outf2<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",f=point"<<endl;//dat
    for(int m=0; m <M; m ++)for(int l=0; l <L; l ++) {outf2<<rbuf2[m+M*l]<<endl;}
    outf2.close(); 
free(temp2);
free(rbuf2);   
  } 

//     char* slicenameE="Eave_", *slicenamePa="phiave", *slicenameP="phi",
//	     *linenamePa="Lave", *linenameP="L";
//     LINEF( dphi, dmyid, dnumprocs, dlocal_n0, dlocal_0_start, linenameP);

}



static void SLICEHat(float *dphi,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start,char *slicename){

  char filename1[28],filename2[28],filename3[28];
  sprintf(filename1,"%s%s%d%s", slicename, "(100)_", kd,".dat");
  sprintf(filename2,"%s%s%d%s", slicename, "(010)_", kd,".dat");
  sprintf(filename3,"%s%s%d%s", slicename, "(001)_", kd,".dat");

  //  float x0=(nx/2.+0.25)*ax, y0=ny/2.*ay, z0=(nz/2.+0.05)*az;
// out 1120 plane
  int proc_mid=0;
  if (dmyid==proc_mid) {     
    ofstream outfile1(filename1,ios_base::out);
    outfile1.setf(ios_base::fixed,ios_base::floatfield);
    outfile1.precision(4);
    int lx=130;    
	outfile1<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outfile1<<"zone t=\"big zone\""<<",i="<<M<<",j="<<N<<",f=point"<<endl;//dat
    for(int n=0; n <N; n ++) for(int m=0; m <M; m ++)  { outfile1<<dphi[n+N*(m+M*lx)]<<endl; }
    outfile1.close(); 
  }
// output 1100 plane   

  float *temp, *rbuf;
  if ( dmyid == 0) {rbuf= (float *) fftw_malloc(sizeof(float)*L*N);} 
  temp= (float *) fftw_malloc(sizeof(float) * N*(dlocal_n0)); 
  int my=167;
  for (int l=0; l <dlocal_n0; l ++) for (int n=0; n <N; n ++) {temp[n+N*l]=dphi[n+N*(my+M*l)];}

  MPI_Gather(temp,N*dlocal_n0,MPI_FLOAT,rbuf,N*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
    ofstream outfile3(filename2,ios_base::out);
    outfile3.setf(ios_base::fixed,ios_base::floatfield);
    outfile3.precision(4);
	outfile3<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outfile3<<"zone t=\"big zone\""<<",i="<<L<<",j="<<N<<",f=point"<<endl;//dat
    for(int n=0; n <N; n ++) for(int l=0; l <L; l ++)  { outfile3<<rbuf[n+N*l]<<endl;}
    outfile3.close();

free(temp);
free(rbuf);   
  } 


// output 0001 plane 

  float *temp1, *rbuf1;
  if ( dmyid == 0) {rbuf1= (float *) fftw_malloc(sizeof(float) *L*M);} 
  temp1= (float *) fftw_malloc(sizeof(float) *M*(dlocal_n0)); 
  int mz=87;
  for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) {temp1[m+M*l]=dphi[mz+N*(m+M*l)];}

  MPI_Gather(temp1,M*dlocal_n0,MPI_FLOAT,rbuf1,M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
    ofstream outfile2(filename3,ios_base::out);
    outfile2.setf(ios_base::fixed,ios_base::floatfield);
    outfile2.precision(4);
	outfile2<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outfile2<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",f=point"<<endl;//dat
    for(int m=0; m <M; m ++) for(int l=0; l <L; l ++) { outfile2<<rbuf[m+M*l]<<endl; }
    outfile2.close(); 

free(temp1);
free(rbuf1);   
  } 
}

static void LINEF(fftw_complex *dphi,fftw_complex *dcon,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start,char *LINEname){
  char filename1[28],filename2[28],filename3[28];
  sprintf(filename1,"%s%s%d%s", LINEname, "[010]_", kd,".txt");
  sprintf(filename2,"%s%s%d%s", LINEname, "[100]_", kd,".txt");
  sprintf(filename3,"%s%s%d%s", LINEname, "[001]_", kd,".txt");

  int proc_mid=(dnumprocs)/2;
  //For[010] direction
    if(dmyid==proc_mid) {     
    ofstream outfile1(filename1,ios_base::out);
    outfile1.setf(ios_base::fixed,ios_base::floatfield);
    outfile1.precision(4);
    int lx=0,nz=N/2;    
    for(int m=0; m <M; m ++)  { outfile1<<m<<" "<<dphi[nz+N*(m+M*lx)][0]<<endl; }
    outfile1.close(); 
  }

//For[100] direction
  
int my,nz;
  float *temp, *rbuf;
  if ( dmyid == 0) {rbuf= (float *) fftw_malloc(sizeof(float)*L);} 
  temp= (float *) fftw_malloc(sizeof(float)*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
     my=M/2;
	 nz=N/2;
     temp[l]=dphi[nz+N*(my+M*l)][0];
  }
  MPI_Gather(temp,dlocal_n0,MPI_FLOAT,rbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
    ofstream outf2(filename2,ios_base::out);
    outf2.setf(ios_base::fixed,ios_base::floatfield);
    outf2.precision(4);
    for(int l=0; l<L; l ++) {outf2<<l<<" "<<rbuf[l]<<endl;}    
    outf2.close();
free(temp);
free(rbuf);   
  } 
/*
//For[001] direction

  if (dmyid==proc_mid) {     
    ofstream outfile3(filename3,ios_base::out);
    outfile3.setf(ios_base::fixed,ios_base::floatfield);
    outfile3.precision(4);
    int lx=0,my=M/2;    
    for(int n=0; n <N; n ++)  { outfile3<<n<<" "<<dphi[n+N*(my+M*lx)][0]<<endl; }
    outfile3.close(); 
  }*/
}

static void OUTPUT(fftw_complex *dphi,fftw_complex *dcon,double *denei,float *dphiave,float *dconave,float *denei_ave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs){
	int startz=0, endz=N;
	float *phi,*deneif,*phiavee,*enei_avee,*con,*conavee;

  phi= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
  con= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
  phiavee= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
  enei_avee= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
  conavee= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));              
  deneif= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));

   for(int l=0; l <dlocal_n0; l ++) for(int m=0; m <M; m ++) for(int n=0; n <N; n ++) {phi[n+N*(m+M*l)]=dphi[n+N*(m+M*l)][0];
   con[n+N*(m+M*l)]=dcon[n+N*(m+M*l)][0];
  deneif[n+N*(m+M*l)]=static_cast<float>(denei[n+N*(m+M*l)]);}              
  //output phi
   char *filename="phi";
  //  sVTKwriteBianry(phi,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,startz,endz,filename);
  //  IF (DMYID==0) { pvtkwrite(DLOCAL_N0,DNUMPROCS,startz,endz,filename); } 
  
  //output enei  
  //filename="enei";
  //sVTKwriteBianry(denei,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,startz,endz,filename);   
 // if (dmyid==0) { PVTKWRITE(dlocal_n0,dnumprocs,startz,endz,filename); }
  
 // output phiave    
  filename="phiave";          
  SMOOTHP(phi,phiavee,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,filename);	
 // output conave    
  filename="conave";          
  SMOOTHC(con,conavee,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,filename);
 // output eneiave
  filename="eneiavee";
  SMOOTHE(deneif,enei_avee,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,filename); 
  free(phi);
  free(deneif);
  free(phiavee);
  free(enei_avee);
  free(conavee);
  free(con);
               
}
/*	
static void SMOOTHP(float *dphi,float *dphiave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename) {
 int nsum,ix,jy,kz,lspace=2, iteration=30;
// int nsum,ix,jy,kz,lspace=2, iteration=20;
 float psum;
 float *buf,*rbuf,*temp;
 float *phi; 
 phi= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 //dphiave= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 buf = (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+2*lspace));
 MPI_Status stat;
 int left=dmyid-1,right=dmyid+1;
 int head=0,tail= dnumprocs-1;
 if (dmyid>0) {left=dmyid-1;}else{left=MPI_PROC_NULL;}
 if (dmyid<dnumprocs-1){right=dmyid+1;}else{right=MPI_PROC_NULL;} 

//"calculate intial value of dphiave"<<endl;
   for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {
       phi[n+N*(m+M*l)]=dphi[n+N*(m+M*l)];
   dphiave[n+N*(m+M*l)]=dphi[n+N*(m+M*l)]; //abs(phi[n+N*(m+M*l)]-u0)
   }

//"input buf of dphiave"<<endl; 
   for(int l=0; l <dlocal_n0+2*lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=u0;   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=dphiave[n+N*(m+M*(l-lspace))];   }

// "smoothening"
   for(int ll=1;ll<iteration+1;ll++){
     if(dmyid==head) {MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,tail,40,&buf[0],N*M*lspace,MPI_FLOAT,tail,40,MPI_COMM_WORLD,&stat);//l1---l0
     }
     if(dmyid==tail) {MPI_Sendrecv(&buf[(dlocal_n0)*N*M],N*M*lspace,MPI_FLOAT,head,40,&buf[(dlocal_n0+lspace)*N*M],N*M*lspace,MPI_FLOAT,head,40,MPI_COMM_WORLD,&stat);//dlocal_n0
     }
     MPI_Sendrecv(&buf[N*M*(dlocal_n0)],N*M*lspace,MPI_FLOAT,right,20,&buf[0],N*M*lspace,MPI_FLOAT,left,20,MPI_COMM_WORLD,&stat);//lnp---l0
     MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,left,20,&buf[N*M*(dlocal_n0+lspace)],N*M*lspace,MPI_FLOAT,right,20,MPI_COMM_WORLD,&stat); //l1---lnp+1
     for (int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for (int n=0; n <N; n ++){
        psum=0.;nsum=0;
        for(int lx=-lspace;lx<lspace+1;lx++) for(int ly=-lspace;ly<lspace+1;ly++) for (int lz=-lspace;lz<lspace+1;lz++){
          ix=l+lx;
          jy=m+ly; 
          kz=n+lz;
          if (jy > (M-1)) {jy=jy-M+1;}
          if (kz > (N-1)) {kz=kz-N+1;}
          if (jy < 0) {jy=jy+M-1;}
          if (kz < 0) {kz=kz+N-1;}
//		  if(ix<2){ix=ix+L-1;}
//        if(ix>(L+1)){ix=ix-L+1;}
          psum=psum+buf[kz+N*(jy+M*ix)];
          nsum=nsum+1; 
	      }
        buf[n+N*(m+M*l)]=psum/nsum; 
     }
   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
      dphiave[n+N*(m+M*(l-lspace))]=buf[n+N*(m+M*l)];
   }
   
   
  if ( dmyid == 0) {rbuf= (float *) fftw_malloc(sizeof(float)*M*L);} 
  temp= (float *) fftw_malloc(sizeof(float) * M*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
int   ny=N/2;
     for (int m=0; m <M; m ++) {temp[m+M*l]= dphiave[ny+N*(m+M*l)];}
  }
  MPI_Gather(temp,M*dlocal_n0,MPI_FLOAT,rbuf,M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
      char filename4[20];
     sprintf(filename4, "%s%s%d%s", dfilename, "ave",kd,".dat"); 
	ofstream outfile2(filename4,ios_base::out);
    outfile2.setf(ios_base::fixed,ios_base::floatfield);
    outfile2.precision(4);
	outfile2<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outfile2<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",f=point"<<endl;//dat
    for(int m=0; m <M; m ++)for(int l=0; l <L; l ++) {outfile2<<rbuf[m+M*l]<<endl;}
    outfile2.close();   
phi_s=rbuf[M/2+M*L/2]; 
phi_l=rbuf[M/2+M*(L-8)];
cout<<"phi_s= "<<phi_s<<" phi_l= "<<phi_l<<endl;
	free(temp);
    free(rbuf);
  }  
   
   
   
   
   // output dphiave
//  int startz=0, endz=N;
//  sVTKwriteBianry(dphiave,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,startz,endz,dfilename);
//  if (dmyid==0) { PVTKWRITE(dlocal_n0,dnumprocs,startz,endz,dfilename);} 
  	
  //output dphiave and phi amplititude plot
     
   float *phix,*phixy,*phixyz,*phiavex,*phiavexy,*phiavexyz,*phixbuf,*phixybuf,*phixyzbuf,*phiavexbuf,*phiavexybuf,*phiavexyzbuf;
   if (dmyid==0){
   	 phixbuf=(float *) fftw_malloc(sizeof(float)*L);
   	 phiavexbuf=(float *) fftw_malloc(sizeof(float)*L);

   } 
   phix=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));
   phiavex=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));

   int nz=N/2, my=M/2;
   for (int l=0; l <dlocal_n0; l ++)  {
       phix[l]=dphi[nz+N*(my+M*l)];
       phiavex[l]=dphiave[nz+N*(my+M*l)];
   }  

     
   MPI_Gather(phix,dlocal_n0,MPI_FLOAT,phixbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  
   MPI_Gather(phiavex,dlocal_n0,MPI_FLOAT,phiavexbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  

    if (dmyid == 0){
     char filename[20];
     sprintf(filename, "%s%d%s", dfilename, kd,".txt");  
     ofstream outfile(filename,ios_base::out);
     outfile.setf(ios_base::fixed,ios_base::floatfield);
     outfile.precision(4);
     cout<<"test"<<endl;
     for(int l=0; l <L; l ++){
//       cout <<l<<" "<<N<<endl;
       outfile<<l*a_100<<" "<<phixbuf[l]<<" "<<phiavexbuf[l]<<endl;
     }
     outfile.close(); 
    }
 free(phi);
  free(buf);
}	

static void SMOOTHC(float *dcon,float *dconave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename) {
 int nsum,ix,jy,kz,lspace=2, iteration=30;
// int nsum,ix,jy,kz,lspace=2, iteration=20;
 float psum;
 float *buf,*rbuf,*temp;
 float *phi; 
 phi= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 //dconave= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 buf = (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+2*lspace));
 MPI_Status stat;
 int left=dmyid-1,right=dmyid+1;
 int head=0,tail= dnumprocs-1;
 if (dmyid>0) {left=dmyid-1;}else{left=MPI_PROC_NULL;}
 if (dmyid<dnumprocs-1){right=dmyid+1;}else{right=MPI_PROC_NULL;} 

//"calculate intial value of dconave"<<endl;
   for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {
       phi[n+N*(m+M*l)]=dcon[n+N*(m+M*l)];
   dconave[n+N*(m+M*l)]=dcon[n+N*(m+M*l)]; //abs(phi[n+N*(m+M*l)]-u0)
   }

//"input buf of dconave"<<endl; 
   for(int l=0; l <dlocal_n0+2*lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=u0;   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=dconave[n+N*(m+M*(l-lspace))];   }

// "smoothening"
   for(int ll=1;ll<iteration+1;ll++){
     if(dmyid==head) {MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,tail,40,&buf[0],N*M*lspace,MPI_FLOAT,tail,40,MPI_COMM_WORLD,&stat);//l1---l0
     }
     if(dmyid==tail) {MPI_Sendrecv(&buf[(dlocal_n0)*N*M],N*M*lspace,MPI_FLOAT,head,40,&buf[(dlocal_n0+lspace)*N*M],N*M*lspace,MPI_FLOAT,head,40,MPI_COMM_WORLD,&stat);//dlocal_n0
     }
     MPI_Sendrecv(&buf[N*M*(dlocal_n0)],N*M*lspace,MPI_FLOAT,right,20,&buf[0],N*M*lspace,MPI_FLOAT,left,20,MPI_COMM_WORLD,&stat);//lnp---l0
     MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,left,20,&buf[N*M*(dlocal_n0+lspace)],N*M*lspace,MPI_FLOAT,right,20,MPI_COMM_WORLD,&stat); //l1---lnp+1
     for (int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for (int n=0; n <N; n ++){
        psum=0.;nsum=0;
        for(int lx=-lspace;lx<lspace+1;lx++) for(int ly=-lspace;ly<lspace+1;ly++) for (int lz=-lspace;lz<lspace+1;lz++){
          ix=l+lx;
          jy=m+ly; 
          kz=n+lz;
          if (jy > (M-1)) {jy=jy-M+1;}
          if (kz > (N-1)) {kz=kz-N+1;}
          if (jy < 0) {jy=jy+M-1;}
          if (kz < 0) {kz=kz+N-1;}
//		  if(ix<2){ix=ix+L-1;}
//        if(ix>(L+1)){ix=ix-L+1;}
          psum=psum+buf[kz+N*(jy+M*ix)];
          nsum=nsum+1; 
	      }
        buf[n+N*(m+M*l)]=psum/nsum; 
     }
   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
      dconave[n+N*(m+M*(l-lspace))]=buf[n+N*(m+M*l)];
   }
   
   
  if ( dmyid == 0) {rbuf= (float *) fftw_malloc(sizeof(float)*M*L);} 
  temp= (float *) fftw_malloc(sizeof(float) * M*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
int   ny=N/2;
     for (int m=0; m <M; m ++) {temp[m+M*l]= dconave[ny+N*(m+M*l)];}
  }
  MPI_Gather(temp,M*dlocal_n0,MPI_FLOAT,rbuf,M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
      char filename4[20];
     sprintf(filename4, "%s%s%d%s", dfilename, "ave",kd,".dat"); 
	ofstream outfile2(filename4,ios_base::out);
    outfile2.setf(ios_base::fixed,ios_base::floatfield);
    outfile2.precision(4);
	outfile2<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outfile2<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",f=point"<<endl;//dat
    for(int m=0; m <M; m ++)for(int l=0; l <L; l ++) {outfile2<<rbuf[m+M*l]<<endl;}
    outfile2.close();   
con_s=rbuf[M/2+M*L/2]; 
con_l=rbuf[M/2+M*(L-8)];
cout<<"con_s= "<<con_s<<" con_l= "<<con_l<<endl;
	free(temp);
    free(rbuf);
  }  
   
   
   
   
   // output dconave
//  int startz=0, endz=N;
//  sVTKwriteBianry(dconave,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,startz,endz,dfilename);
//  if (dmyid==0) { PVTKWRITE(dlocal_n0,dnumprocs,startz,endz,dfilename);} 
  	
  //output dconave and phi amplititude plot
     
   float *phix,*phixy,*phixyz,*phiavex,*phiavexy,*phiavexyz,*phixbuf,*phixybuf,*phixyzbuf,*phiavexbuf,*phiavexybuf,*phiavexyzbuf;
   if (dmyid==0){
   	 phixbuf=(float *) fftw_malloc(sizeof(float)*L);
   	 phiavexbuf=(float *) fftw_malloc(sizeof(float)*L);

   } 
   phix=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));
   phiavex=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));

   int nz=N/2, my=M/2;
   for (int l=0; l <dlocal_n0; l ++)  {
       phix[l]=dcon[nz+N*(my+M*l)];
       phiavex[l]=dconave[nz+N*(my+M*l)];
   }  

     
   MPI_Gather(phix,dlocal_n0,MPI_FLOAT,phixbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  
   MPI_Gather(phiavex,dlocal_n0,MPI_FLOAT,phiavexbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  

    if (dmyid == 0){
     char filename[20];
     sprintf(filename, "%s%d%s", dfilename, kd,".txt");  
     ofstream outfile(filename,ios_base::out);
     outfile.setf(ios_base::fixed,ios_base::floatfield);
     outfile.precision(4);
     cout<<"test"<<endl;
     for(int l=0; l <L; l ++){
//       cout <<l<<" "<<N<<endl;
       outfile<<l*a_100<<" "<<phixbuf[l]<<" "<<phiavexbuf[l]<<endl;
     }
     outfile.close(); 
    }
 free(phi);
  free(buf);
}	

static void SMOOTHE(float *dphi,float *denei_ave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename) {
 int nsum,ix,jy,kz,lspace=2, iteration=30;
// int nsum,ix,jy,kz,lspace=2, iteration=20;
 float psum;
 float *buf,*rbuf,*temp;
 //phiave= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 buf = (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+2*lspace));
 MPI_Status stat;
 int left=dmyid-1,right=dmyid+1;
 int head=0,tail= dnumprocs-1;
 if (dmyid>0) {left=dmyid-1;}else{left=MPI_PROC_NULL;}
 if (dmyid<dnumprocs-1){right=dmyid+1;}else{right=MPI_PROC_NULL;} 

//"calculate intial value of phiave"<<endl;
   for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {
              denei_ave[n+N*(m+M*l)]=dphi[n+N*(m+M*l)]; //abs(phi[n+N*(m+M*l)]-u0)
   }

//"input buf of denei_ave"<<endl; 
   for(int l=0; l <dlocal_n0+2*lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=u0;   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=denei_ave[n+N*(m+M*(l-lspace))];   }

// "smoothening"
   for(int ll=1;ll<iteration+1;ll++){
     if(dmyid==head) {MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,tail,40,&buf[0],N*M*lspace,MPI_FLOAT,tail,40,MPI_COMM_WORLD,&stat);//l1---l0
     }
     if(dmyid==tail) {MPI_Sendrecv(&buf[(dlocal_n0)*N*M],N*M*lspace,MPI_FLOAT,head,40,&buf[(dlocal_n0+lspace)*N*M],N*M*lspace,MPI_FLOAT,head,40,MPI_COMM_WORLD,&stat);//dlocal_n0
     }
     MPI_Sendrecv(&buf[N*M*(dlocal_n0)],N*M*lspace,MPI_FLOAT,right,20,&buf[0],N*M*lspace,MPI_FLOAT,left,20,MPI_COMM_WORLD,&stat);//lnp---l0
     MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,left,20,&buf[N*M*(dlocal_n0+lspace)],N*M*lspace,MPI_FLOAT,right,20,MPI_COMM_WORLD,&stat); //l1---lnp+1
     for (int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for (int n=0; n <N; n ++){
        psum=0.;nsum=0;
        for(int lx=-lspace;lx<lspace+1;lx++) for(int ly=-lspace;ly<lspace+1;ly++) for (int lz=-lspace;lz<lspace+1;lz++){
          ix=l+lx;
          jy=m+ly; 
          kz=n+lz;
          if (jy > (M-1)) {jy=jy-M+1;}
          if (kz > (N-1)) {kz=kz-N+1;}
          if (jy < 0) {jy=jy+M-1;}
          if (kz < 0) {kz=kz+N-1;}
//		  if(ix<2){ix=ix+L-1;}
//        if(ix>(L+1)){ix=ix-L+1;}
          psum=psum+buf[kz+N*(jy+M*ix)];
          nsum=nsum+1; 
	      }
        buf[n+N*(m+M*l)]=psum/nsum; 
     }
   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
      denei_ave[n+N*(m+M*(l-lspace))]=buf[n+N*(m+M*l)];
   }


   
   
  if ( dmyid == 0) {rbuf= (float *) fftw_malloc(sizeof(float)*M*L);} 
  temp= (float *) fftw_malloc(sizeof(float) * M*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
int   ny=N/2;
     for (int m=0; m <M; m ++) {temp[m+M*l]= denei_ave[ny+N*(m+M*l)];}
  }
  MPI_Gather(temp,M*dlocal_n0,MPI_FLOAT,rbuf,M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
      char filename4[20];
     sprintf(filename4, "%s%s%d%s", dfilename, "ave",kd,".dat"); 
	ofstream outfile2(filename4,ios_base::out);
    outfile2.setf(ios_base::fixed,ios_base::floatfield);
    outfile2.precision(4);
	outfile2<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outfile2<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",f=point"<<endl;//dat
    for(int m=0; m <M; m ++)for(int l=0; l <L; l ++) {outfile2<<rbuf[m+M*l]<<endl;}
    outfile2.close();   
   f_s=rbuf[M/2+M*L/2]; 
   f_l=rbuf[M/2+M*(L-8)];
   cout<<" f_s= "<<f_s<<"f_l= "<<f_l<<endl;
	free(temp);
    free(rbuf);
  }  
   
   
   
   
   // output denei_ave
//  int startz=0, endz=N;
//  sVTKwriteBianry(denei_ave,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,startz,endz,dfilename);
//  if (dmyid==0) { PVTKWRITE(dlocal_n0,dnumprocs,startz,endz,dfilename);} 
  	
  //output denei_ave and phi amplititude plot
     
   float *phix,*phixy,*phixyz,*phiavex,*phiavexy,*phiavexyz,*phixbuf,*phixybuf,*phixyzbuf,*phiavexbuf,*phiavexybuf,*phiavexyzbuf;
   if (dmyid==0){
   	 phixbuf=(float *) fftw_malloc(sizeof(float)*L);
   	 phiavexbuf=(float *) fftw_malloc(sizeof(float)*L);

   } 
   phix=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));
   phiavex=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));

   int nz=N/2, my=M/2;
   for (int l=0; l <dlocal_n0; l ++)  {
       phix[l]=dphi[nz+N*(my+M*l)];
       phiavex[l]=denei_ave[nz+N*(my+M*l)];
   }  

     
   MPI_Gather(phix,dlocal_n0,MPI_FLOAT,phixbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  
   MPI_Gather(phiavex,dlocal_n0,MPI_FLOAT,phiavexbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  

    if (dmyid == 0){
     char filename[20];
     sprintf(filename, "%s%d%s", dfilename, kd,".txt");  
     ofstream outfile(filename,ios_base::out);
     outfile.setf(ios_base::fixed,ios_base::floatfield);
     outfile.precision(4);
     cout<<"test"<<endl;
     for(int l=0; l <L; l ++){
//       cout <<l<<" "<<N<<endl;
       outfile<<l*a_100<<" "<<phixbuf[l]<<" "<<phiavexbuf[l]<<endl;
     }
     outfile.close(); 
    }
//  free(phiave);
  free(buf);
}	
*/

static void SMOOTHP(float *dphi,float *dphiave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename) {
 int nsum,ix,jy,kz,lspace=2, iteration=30;
// int nsum,ix,jy,kz,lspace=2, iteration=20;
 float psum;
 float *buf,*rbuf,*temp;
 float *phi; 
 phi= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 //dphiave= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 buf = (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+2*lspace));
 MPI_Status stat;
 int left=dmyid-1,right=dmyid+1;
 int head=0,tail= dnumprocs-1;
 if (dmyid>0) {left=dmyid-1;}else{left=MPI_PROC_NULL;}
 if (dmyid<dnumprocs-1){right=dmyid+1;}else{right=MPI_PROC_NULL;} 

//"calculate intial value of dphiave"<<endl;
   for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {
       phi[n+N*(m+M*l)]=dphi[n+N*(m+M*l)];
   dphiave[n+N*(m+M*l)]=dphi[n+N*(m+M*l)]; //abs(phi[n+N*(m+M*l)]-u0)
   }

//"input buf of dphiave"<<endl; 
   for(int l=0; l <dlocal_n0+2*lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=u0;   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=dphiave[n+N*(m+M*(l-lspace))];   }

// "smoothening"
   for(int ll=1;ll<iteration+1;ll++){
     if(dmyid==head) {MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,tail,40,&buf[0],N*M*lspace,MPI_FLOAT,tail,40,MPI_COMM_WORLD,&stat);//l1---l0
     }
     if(dmyid==tail) {MPI_Sendrecv(&buf[(dlocal_n0)*N*M],N*M*lspace,MPI_FLOAT,head,40,&buf[(dlocal_n0+lspace)*N*M],N*M*lspace,MPI_FLOAT,head,40,MPI_COMM_WORLD,&stat);//dlocal_n0
     }
     MPI_Sendrecv(&buf[N*M*(dlocal_n0)],N*M*lspace,MPI_FLOAT,right,20,&buf[0],N*M*lspace,MPI_FLOAT,left,20,MPI_COMM_WORLD,&stat);//lnp---l0
     MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,left,20,&buf[N*M*(dlocal_n0+lspace)],N*M*lspace,MPI_FLOAT,right,20,MPI_COMM_WORLD,&stat); //l1---lnp+1
     for (int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for (int n=0; n <N; n ++){
        psum=0.;nsum=0;
        for(int lx=-lspace;lx<lspace+1;lx++) for(int ly=-lspace;ly<lspace+1;ly++) {
          ix=l+lx;
          jy=m+ly; 
          kz=0;
          if (jy > (M-1)) {jy=jy-M+1;}
//          if (kz > (N-1)) {kz=kz-N+1;}
          if (jy < 0) {jy=jy+M-1;}
//          if (kz < 0) {kz=kz+N-1;}
//		  if(ix<2){ix=ix+L-1;}
//        if(ix>(L+1)){ix=ix-L+1;}
          psum=psum+buf[kz+N*(jy+M*ix)];
          nsum=nsum+1; 
	      }
        buf[n+N*(m+M*l)]=psum/nsum; 
     }
   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
      dphiave[n+N*(m+M*(l-lspace))]=buf[n+N*(m+M*l)];
   }
   
   
  if ( dmyid == 0) {rbuf= (float *) fftw_malloc(sizeof(float)*M*L);} 
  temp= (float *) fftw_malloc(sizeof(float) * M*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
int   ny=N/2;
     for (int m=0; m <M; m ++) {temp[m+M*l]= dphiave[ny+N*(m+M*l)];}
  }
  MPI_Gather(temp,M*dlocal_n0,MPI_FLOAT,rbuf,M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
      char filename4[20];
     sprintf(filename4, "%s%s%d%s", dfilename, "ave",kd,".dat"); 
	ofstream outfile2(filename4,ios_base::out);
    outfile2.setf(ios_base::fixed,ios_base::floatfield);
    outfile2.precision(4);
	outfile2<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outfile2<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",f=point"<<endl;//dat
    for(int m=0; m <M; m ++)for(int l=0; l <L; l ++) {outfile2<<rbuf[m+M*l]<<endl;}
    outfile2.close();   
phi_s=rbuf[M/2+M*L/2]; 
phi_l=rbuf[M/2+M*(L-8)];
cout<<"phi_s= "<<phi_s<<" phi_l= "<<phi_l<<endl;
	free(temp);
    free(rbuf);
  }  
   
   
   
   
   // output dphiave
//  int startz=0, endz=N;
//  sVTKwriteBianry(dphiave,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,startz,endz,dfilename);
//  if (dmyid==0) { PVTKWRITE(dlocal_n0,dnumprocs,startz,endz,dfilename);} 
  	
  //output dphiave and phi amplititude plot
     
   float *phix,*phixy,*phixyz,*phiavex,*phiavexy,*phiavexyz,*phixbuf,*phixybuf,*phixyzbuf,*phiavexbuf,*phiavexybuf,*phiavexyzbuf;
   if (dmyid==0){
   	 phixbuf=(float *) fftw_malloc(sizeof(float)*L);
   	 phiavexbuf=(float *) fftw_malloc(sizeof(float)*L);

   } 
   phix=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));
   phiavex=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));

   int nz=N/2, my=M/2;
   for (int l=0; l <dlocal_n0; l ++)  {
       phix[l]=dphi[nz+N*(my+M*l)];
       phiavex[l]=dphiave[nz+N*(my+M*l)];
   }  

     
   MPI_Gather(phix,dlocal_n0,MPI_FLOAT,phixbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  
   MPI_Gather(phiavex,dlocal_n0,MPI_FLOAT,phiavexbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  

    if (dmyid == 0){
     char filename[20];
     sprintf(filename, "%s%d%s", dfilename, kd,".txt");  
     ofstream outfile(filename,ios_base::out);
     outfile.setf(ios_base::fixed,ios_base::floatfield);
     outfile.precision(4);
     cout<<"test"<<endl;
     for(int l=0; l <L; l ++){
//       cout <<l<<" "<<N<<endl;
       outfile<<l*a_100<<" "<<phixbuf[l]<<" "<<phiavexbuf[l]<<endl;
     }
     outfile.close(); 
    }
 free(phi);
  free(buf);
}	

static void SMOOTHC(float *dcon,float *dconave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename) {
 int nsum,ix,jy,kz,lspace=2, iteration=30;
// int nsum,ix,jy,kz,lspace=2, iteration=20;
 float psum;
 float *buf,*rbuf,*temp;
 float *phi; 
 phi= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 //dconave= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 buf = (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+2*lspace));
 MPI_Status stat;
 int left=dmyid-1,right=dmyid+1;
 int head=0,tail= dnumprocs-1;
 if (dmyid>0) {left=dmyid-1;}else{left=MPI_PROC_NULL;}
 if (dmyid<dnumprocs-1){right=dmyid+1;}else{right=MPI_PROC_NULL;} 

//"calculate intial value of dconave"<<endl;
   for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {
       phi[n+N*(m+M*l)]=dcon[n+N*(m+M*l)];
   dconave[n+N*(m+M*l)]=dcon[n+N*(m+M*l)]; //abs(phi[n+N*(m+M*l)]-u0)
   }

//"input buf of dconave"<<endl; 
   for(int l=0; l <dlocal_n0+2*lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=con0;   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=dconave[n+N*(m+M*(l-lspace))];   }

// "smoothening"
   for(int ll=1;ll<iteration+1;ll++){
     if(dmyid==head) {MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,tail,40,&buf[0],N*M*lspace,MPI_FLOAT,tail,40,MPI_COMM_WORLD,&stat);//l1---l0
     }
     if(dmyid==tail) {MPI_Sendrecv(&buf[(dlocal_n0)*N*M],N*M*lspace,MPI_FLOAT,head,40,&buf[(dlocal_n0+lspace)*N*M],N*M*lspace,MPI_FLOAT,head,40,MPI_COMM_WORLD,&stat);//dlocal_n0
     }
     MPI_Sendrecv(&buf[N*M*(dlocal_n0)],N*M*lspace,MPI_FLOAT,right,20,&buf[0],N*M*lspace,MPI_FLOAT,left,20,MPI_COMM_WORLD,&stat);//lnp---l0
     MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,left,20,&buf[N*M*(dlocal_n0+lspace)],N*M*lspace,MPI_FLOAT,right,20,MPI_COMM_WORLD,&stat); //l1---lnp+1
     for (int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++){
        psum=0.;nsum=0;int n=0;
        for(int lx=-lspace;lx<lspace+1;lx++) for(int ly=-lspace;ly<lspace+1;ly++) {
          ix=l+lx;
          jy=m+ly; 
          kz=0;
          if (jy > (M-1)) {jy=jy-M+1;}
   //       if (kz > (N-1)) {kz=kz-N+1;}
          if (jy < 0) {jy=jy+M-1;}
   //       if (kz < 0) {kz=kz+N-1;}
//		  if(ix<2){ix=ix+L-1;}
//        if(ix>(L+1)){ix=ix-L+1;}
          psum=psum+buf[kz+N*(jy+M*ix)];
          nsum=nsum+1; 
	      }
        buf[n+N*(m+M*l)]=psum/nsum; 
     }
   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
      dconave[n+N*(m+M*(l-lspace))]=buf[n+N*(m+M*l)];
   }
   
   
  if ( dmyid == 0) {rbuf= (float *) fftw_malloc(sizeof(float)*M*L);} 
  temp= (float *) fftw_malloc(sizeof(float) * M*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
int   ny=N/2;
     for (int m=0; m <M; m ++) {temp[m+M*l]= dconave[ny+N*(m+M*l)];}
  }
  MPI_Gather(temp,M*dlocal_n0,MPI_FLOAT,rbuf,M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
      char filename4[20];
     sprintf(filename4, "%s%s%d%s", dfilename, "ave",kd,".dat"); 
	ofstream outfile2(filename4,ios_base::out);
    outfile2.setf(ios_base::fixed,ios_base::floatfield);
    outfile2.precision(4);
	outfile2<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outfile2<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",f=point"<<endl;//dat
    for(int m=0; m <M; m ++)for(int l=0; l <L; l ++) {outfile2<<rbuf[m+M*l]<<endl;}
    outfile2.close();   
con_s=rbuf[M/2+M*L/2]; 
con_l=rbuf[M/2+M*(L-8)];
cout<<"con_s= "<<con_s<<" con_l= "<<con_l<<endl;
	free(temp);
    free(rbuf);
  }  
   
   
   
   
   // output dconave
//  int startz=0, endz=N;
//  sVTKwriteBianry(dconave,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,startz,endz,dfilename);
//  if (dmyid==0) { PVTKWRITE(dlocal_n0,dnumprocs,startz,endz,dfilename);} 
  	
  //output dconave and phi amplititude plot
     
   float *phix,*phixy,*phixyz,*phiavex,*phiavexy,*phiavexyz,*phixbuf,*phixybuf,*phixyzbuf,*phiavexbuf,*phiavexybuf,*phiavexyzbuf;
   if (dmyid==0){
   	 phixbuf=(float *) fftw_malloc(sizeof(float)*L);
   	 phiavexbuf=(float *) fftw_malloc(sizeof(float)*L);

   } 
   phix=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));
   phiavex=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));

   int nz=N/2, my=M/2;
   for (int l=0; l <dlocal_n0; l ++)  {
       phix[l]=dcon[nz+N*(my+M*l)];
       phiavex[l]=dconave[nz+N*(my+M*l)];
   }  

     
   MPI_Gather(phix,dlocal_n0,MPI_FLOAT,phixbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  
   MPI_Gather(phiavex,dlocal_n0,MPI_FLOAT,phiavexbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  

    if (dmyid == 0){
     char filename[20];
     sprintf(filename, "%s%d%s", dfilename, kd,".txt");  
     ofstream outfile(filename,ios_base::out);
     outfile.setf(ios_base::fixed,ios_base::floatfield);
     outfile.precision(4);
     cout<<"test"<<endl;
     for(int l=0; l <L; l ++){
//       cout <<l<<" "<<N<<endl;
       outfile<<l*a_100<<" "<<phixbuf[l]<<" "<<phiavexbuf[l]<<endl;
     }
     outfile.close(); 
    }
 free(phi);
  free(buf);
}	

static void SMOOTHE(float *dphi,float *denei_ave,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,char *dfilename) {
 int nsum,ix,jy,kz,lspace=2, iteration=30;
// int nsum,ix,jy,kz,lspace=2, iteration=20;
 float psum;
 float *buf,*rbuf,*temp;
 //phiave= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 buf = (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+2*lspace));
 MPI_Status stat;
 int left=dmyid-1,right=dmyid+1;
 int head=0,tail= dnumprocs-1;
 if (dmyid>0) {left=dmyid-1;}else{left=MPI_PROC_NULL;}
 if (dmyid<dnumprocs-1){right=dmyid+1;}else{right=MPI_PROC_NULL;} 

//"calculate intial value of phiave"<<endl;
   for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {
              denei_ave[n+N*(m+M*l)]=dphi[n+N*(m+M*l)]; //abs(phi[n+N*(m+M*l)]-u0)
   }

//"input buf of denei_ave"<<endl; 
   for(int l=0; l <dlocal_n0+2*lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=u0;   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=denei_ave[n+N*(m+M*(l-lspace))];   }

// "smoothening"
   for(int ll=1;ll<iteration+1;ll++){
     if(dmyid==head) {MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,tail,40,&buf[0],N*M*lspace,MPI_FLOAT,tail,40,MPI_COMM_WORLD,&stat);//l1---l0
     }
     if(dmyid==tail) {MPI_Sendrecv(&buf[(dlocal_n0)*N*M],N*M*lspace,MPI_FLOAT,head,40,&buf[(dlocal_n0+lspace)*N*M],N*M*lspace,MPI_FLOAT,head,40,MPI_COMM_WORLD,&stat);//dlocal_n0
     }
     MPI_Sendrecv(&buf[N*M*(dlocal_n0)],N*M*lspace,MPI_FLOAT,right,20,&buf[0],N*M*lspace,MPI_FLOAT,left,20,MPI_COMM_WORLD,&stat);//lnp---l0
     MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,left,20,&buf[N*M*(dlocal_n0+lspace)],N*M*lspace,MPI_FLOAT,right,20,MPI_COMM_WORLD,&stat); //l1---lnp+1
     for (int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) {
        psum=0.;nsum=0;int n=0;
        for(int lx=-lspace;lx<lspace+1;lx++) for(int ly=-lspace;ly<lspace+1;ly++) {
          ix=l+lx;
          jy=m+ly; 
          kz=0;
          if (jy > (M-1)) {jy=jy-M+1;}
 //         if (kz > (N-1)) {kz=kz-N+1;}
          if (jy < 0) {jy=jy+M-1;}
  //        if (kz < 0) {kz=kz+N-1;}
//		  if(ix<2){ix=ix+L-1;}
//        if(ix>(L+1)){ix=ix-L+1;}
          psum=psum+buf[kz+N*(jy+M*ix)];
          nsum=nsum+1; 
	      }
        buf[n+N*(m+M*l)]=psum/nsum; 
     }
   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
      denei_ave[n+N*(m+M*(l-lspace))]=buf[n+N*(m+M*l)];
   }


   
   
  if ( dmyid == 0) {rbuf= (float *) fftw_malloc(sizeof(float)*M*L);} 
  temp= (float *) fftw_malloc(sizeof(float) * M*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  {
int   ny=N/2;
     for (int m=0; m <M; m ++) {temp[m+M*l]= denei_ave[ny+N*(m+M*l)];}
  }
  MPI_Gather(temp,M*dlocal_n0,MPI_FLOAT,rbuf,M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 
      char filename4[20];
     sprintf(filename4, "%s%s%d%s", dfilename, "ave",kd,".dat"); 
	ofstream outfile2(filename4,ios_base::out);
    outfile2.setf(ios_base::fixed,ios_base::floatfield);
    outfile2.precision(4);
	outfile2<<"variable=\"x\",\"y\",\"function\""<<endl; //dat
	outfile2<<"zone t=\"big zone\""<<",i="<<L<<",j="<<M<<",f=point"<<endl;//dat
    for(int m=0; m <M; m ++)for(int l=0; l <L; l ++) {outfile2<<rbuf[m+M*l]<<endl;}
    outfile2.close();   
   f_s=rbuf[M/2+M*L/2]; 
   f_l=rbuf[M/2+M*(L-8)];
   cout<<" f_s= "<<f_s<<"f_l= "<<f_l<<endl;
	free(temp);
    free(rbuf);
  }  
   
   
   
   
   // output denei_ave
//  int startz=0, endz=N;
//  sVTKwriteBianry(denei_ave,dlocal_n0,dlocal_0_start,dmyid,dnumprocs,startz,endz,dfilename);
//  if (dmyid==0) { PVTKWRITE(dlocal_n0,dnumprocs,startz,endz,dfilename);} 
  	
  //output denei_ave and phi amplititude plot
     
   float *phix,*phixy,*phixyz,*phiavex,*phiavexy,*phiavexyz,*phixbuf,*phixybuf,*phixyzbuf,*phiavexbuf,*phiavexybuf,*phiavexyzbuf;
   if (dmyid==0){
   	 phixbuf=(float *) fftw_malloc(sizeof(float)*L);
   	 phiavexbuf=(float *) fftw_malloc(sizeof(float)*L);

   } 
   phix=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));
   phiavex=(float *) fftw_malloc(sizeof(float)*(dlocal_n0));

   int nz=N/2, my=M/2;
   for (int l=0; l <dlocal_n0; l ++)  {
       phix[l]=dphi[nz+N*(my+M*l)];
       phiavex[l]=denei_ave[nz+N*(my+M*l)];
   }  

     
   MPI_Gather(phix,dlocal_n0,MPI_FLOAT,phixbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  
   MPI_Gather(phiavex,dlocal_n0,MPI_FLOAT,phiavexbuf,dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD);  

    if (dmyid == 0){
     char filename[20];
     sprintf(filename, "%s%d%s", dfilename, kd,".txt");  
     ofstream outfile(filename,ios_base::out);
     outfile.setf(ios_base::fixed,ios_base::floatfield);
     outfile.precision(4);
     cout<<"test"<<endl;
     for(int l=0; l <L; l ++){
//       cout <<l<<" "<<N<<endl;
       outfile<<l*a_100<<" "<<phixbuf[l]<<" "<<phiavexbuf[l]<<endl;
     }
     outfile.close(); 
    }
//  free(phiave);
  free(buf);
}	

static void INTERFACEM(fftw_complex *dphi,int dicd,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs,int &dintertop,int &dinterlow){

 float *dinterlowm,*dintertopm; 
 dinterlowm=(float *) fftw_malloc(sizeof(float) * M*N); 
 dintertopm=(float *) fftw_malloc(sizeof(float) * M*N);

 
 int intertop=0,interlow=0,
	 bufintertop,bufinterlow; 
   MPI_Status stat;
 float *temp2, *rbuf2;
  if ( dmyid == 0) {rbuf2= (float *) fftw_malloc(sizeof(float)*L*M*N);} 
  temp2= (float *) fftw_malloc(sizeof(float) * M*N*(dlocal_n0)); 
  for (int l=0; l <dlocal_n0; l ++)  for (int m=0; m<M; m++) for (int n=0; n <N; n ++) {temp2[n+N*(m+M*l)]=dphi[n+N*(m+M*l)][0];}

  MPI_Gather(temp2,N*M*dlocal_n0,MPI_FLOAT,rbuf2,N*M*dlocal_n0,MPI_FLOAT,0,MPI_COMM_WORLD); 
  if ( dmyid == 0) { 

  float  phic=u0+0.05;  


   int num_low=1;
   int num_top=1;
   
   for (int m=0; m <M; m ++)for (int n=0; n <N; n ++)  {

       for (int l=0; l <L/2; l ++){
           if (rbuf2[n+N*(m+M*l)]>phic){ dinterlowm[n+N*m]=l; num_low++;
             goto endlow;  
          }
        }
       dinterlowm[n+N*m]=0;
       endlow:;
       for (int l=L-1; l >L/2; l --) {
          if (rbuf2[n+N*(m+M*l)]>phic){ dintertopm[n+N*m]=l; num_top++;
            goto endtop;
          }  
        }
        dintertopm[n+N*m]=0;
        endtop:;         
        dintertopm[n+N*m]=dintertopm[n+N*m];
        dinterlowm[n+N*m]=dinterlowm[n+N*m];

   }

   for (int m=0; m <M; m ++)for (int n=0; n <N; n ++)  {
        intertop += int(dintertopm[n+N*m]);
        interlow += int(dinterlowm[n+N*m]);
}
        intertop /=num_top;
        interlow /=num_low;



  dintertop=intertop+dicd;
  dinterlow=interlow-dicd;
  if(dintertop>L-1) {dintertop=L-1;}
  if (dinterlow<0)  {dinterlow=0;}


  
  
    ofstream wafile("face.txt",ios_base::out|ios_base::app);
	wafile << dintertop <<" "<<dinterlow<<endl;



  free(dintertopm);
  free(dinterlowm);
  free(temp2);
  free(rbuf2);
}}


static void interface_energy(fftw_complex *dphi,double *denei,float *dphiave,float *denei_ave,int dmyid,int dnumprocs,int dlocal_n0,int dlocal_0_start){

       float energy_mix,energy_excess,energy_mix1,energy_excess1;
      float  sum_excess=0., sum_excess1=0; 
        int lnew;
	for (int l=0; l <dlocal_n0; l ++){
	   for (int m=0; m <M; m ++){
	       for (int n=0; n <N; n ++){
                    if(dphiave[n+N*(m+M*l)]<phi_s && dphiave[n+N*(m+M*l)]>phi_l)
                     energy_mix=f_s*(dphiave[n+N*(m+M*l)]-phi_l)/(phi_s-phi_l)-f_l*(dphiave[n+N*(m+M*l)]-phi_s)/(phi_s-phi_l);
                     energy_excess=denei_ave[n+N*(m+M*l)]-energy_mix;
                     sum_excess=sum_excess+energy_excess;
                       
                        energy_mix1=f_s*(dphi[n+N*(m+M*l)][0]-phi_l)/(phi_s-phi_l)-f_l*(dphi[n+N*(m+M*l)][0]-phi_s)/(phi_s-phi_l);
                     energy_excess1=denei[n+N*(m+M*l)]-energy_mix;
                     sum_excess1=sum_excess+energy_excess1;





				}
			}
		}


   double * rcounts,* rcounts1;
   rcounts = (double *) malloc(dnumprocs*sizeof(double));
   rcounts1 = (double *) malloc(dnumprocs*sizeof(double));
   MPI_Gather(&sum_excess,1,MPI_FLOAT,rcounts,1,MPI_FLOAT,0,MPI_COMM_WORLD);
   MPI_Gather(&sum_excess1,1,MPI_FLOAT,rcounts1,1,MPI_FLOAT,0,MPI_COMM_WORLD);

if(dmyid==0) { 
  float sum_ex=0,sum_ex1=0;
  for(int i=0;i<dnumprocs;i++){

    sum_ex=rcounts[i]+sum_ex;
    sum_ex1=rcounts1[i]+sum_ex1;}

  // cout<<"total_atoms="<<sum_atom<<"kd"<<kd<<endl;

   ofstream wafile("sum_Inter.txt",ios_base::out|ios_base::app);
   wafile << kd << " " << sum_ex <<" "<<sum_ex1<<endl;}

free(rcounts1);
free(rcounts);
}

static void inirotationzb(float dvtr2[9], float drotmatrix[9]){
	float dvtr1[9]={ 1,0,0, 0,1,0, 0,0,1 };		  
	float dc11,dc12,dc13,dc21,dc22,dc23,dc31,dc32,dc33;
	float dmagn2_x, dmagn2_y, dmagn2_z;
	    dmagn2_x=sqrt(dvtr2[0]*dvtr2[0]+dvtr2[1]*dvtr2[1]+dvtr2[2]*dvtr2[2]);
		dmagn2_y=sqrt(dvtr2[3]*dvtr2[3]+dvtr2[4]*dvtr2[4]+dvtr2[5]*dvtr2[5]);
		dmagn2_z=sqrt(dvtr2[6]*dvtr2[6]+dvtr2[7]*dvtr2[7]+dvtr2[8]*dvtr2[8]);

		dc11=dvtr1[0]*dvtr2[0]+dvtr1[1]*dvtr2[1]+dvtr1[2]*dvtr2[2];
		dc12=dvtr1[0]*dvtr2[3]+dvtr1[1]*dvtr2[4]+dvtr1[2]*dvtr2[5];
		dc13=dvtr1[0]*dvtr2[6]+dvtr1[1]*dvtr2[7]+dvtr1[2]*dvtr2[8];

		dc21=dvtr1[3]*dvtr2[0]+dvtr1[4]*dvtr2[1]+dvtr1[5]*dvtr2[2];
		dc22=dvtr1[3]*dvtr2[3]+dvtr1[4]*dvtr2[4]+dvtr1[5]*dvtr2[5];
		dc23=dvtr1[3]*dvtr2[6]+dvtr1[4]*dvtr2[7]+dvtr1[5]*dvtr2[8];

		dc31=dvtr1[6]*dvtr2[0]+dvtr1[7]*dvtr2[1]+dvtr1[8]*dvtr2[2];
		dc32=dvtr1[6]*dvtr2[3]+dvtr1[7]*dvtr2[4]+dvtr1[8]*dvtr2[5];
		dc33=dvtr1[6]*dvtr2[6]+dvtr1[7]*dvtr2[7]+dvtr1[8]*dvtr2[8];
	drotmatrix[0]=(dc11)/(dmagn2_x);
	drotmatrix[1]=(dc12)/(dmagn2_y);
	drotmatrix[2]=(dc13)/(dmagn2_z);
	drotmatrix[3]=(dc21)/(dmagn2_x);
	drotmatrix[4]=(dc22)/(dmagn2_y);
	drotmatrix[5]=(dc23)/(dmagn2_z);
	drotmatrix[6]=(dc31)/(dmagn2_x);
	drotmatrix[7]=(dc32)/(dmagn2_y);
	drotmatrix[8]=(dc33)/(dmagn2_z);
}

static void ROT_XYZ(float dang_base[9], float dang[9], float *dtheta, char *daxis){    
	float ta=*dtheta; //ȡֵ
//	if(daxis=="x") float drot[9]={cos(ta), 0., -sin(ta), 0., 1., 0., sin(ta), 0., cos(ta)};  
	float drot[9]={cos(ta), 0., -sin(ta), 0., 1., 0., sin(ta), 0., cos(ta)};  
//	if(daxis=="z") float drot[9]={cos(ta), 0., -sin(ta), 0., 1., 0., sin(ta), 0., cos(ta)};  
	dang[0]=dang_base[0]*drot[0]+dang_base[1]*drot[3]+dang_base[2]*drot[6];
	dang[1]=dang_base[0]*drot[1]+dang_base[1]*drot[4]+dang_base[2]*drot[7];
    dang[2]=dang_base[0]*drot[2]+dang_base[1]*drot[5]+dang_base[2]*drot[8];

	dang[3]=dang_base[3]*drot[0]+dang_base[4]*drot[3]+dang_base[5]*drot[6];
	dang[4]=dang_base[3]*drot[1]+dang_base[4]*drot[4]+dang_base[5]*drot[7];
	dang[5]=dang_base[3]*drot[2]+dang_base[4]*drot[5]+dang_base[5]*drot[8];

	dang[6]=dang_base[6]*drot[0]+dang_base[7]*drot[3]+dang_base[8]*drot[6];
	dang[7]=dang_base[6]*drot[1]+dang_base[7]*drot[4]+dang_base[8]*drot[7];
	dang[8]=dang_base[6]*drot[2]+dang_base[7]*drot[5]+dang_base[8]*drot[8];
}


static void READPHIVTIN(fftw_complex *dphi, int dlocal_n0,int dmyid){
     
     char filename[20];
	 sprintf(filename, "%s%d%s%d%s", "phi_",111001,"_",dmyid,".vti");

	 float *treal,*ddphi;
	 treal=(float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+1)); //[N][M][L];
	 ddphi=(float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+1));


     FILE * inf_field;
     inf_field=fopen(filename,"r");
	 if(inf_field==NULL){
	   printf("can't open this file\n");
		   exit(0);
	 }
     char bufc[100]; //�洢�ַ�
     int i = 1;
     while (i<13) {
       fgets(bufc, 100, inf_field) ;
       printf("%d %s", i++, bufc);
     }  
     char ch;
     ch=fgetc(inf_field);
     cout<<"ch="<<ch<<endl;

     int bytes;
     fread(&bytes,sizeof(int),1,inf_field);
//     cout<<"dmyid="<<dmyid<<" "<<bytes<<endl;
	 cout<<bytes<<endl;  
	 int sumi=0;
     if(dmyid==0){           
	                for(int n=0; n<N; n++){
                    for(int m=0; m<M; m++){
                    for(int l=0; l<dlocal_n0; l++){
						 // if(fread(&treal[l+dlocal_n0*(m+M*n)],sizeof(float),1,inf_field)!=1){  // //[n][m][l]
							if(fread(&treal[n+N*(m+M*l)],sizeof(float),1,inf_field)!=1){  // //[n][m][l]
							cout<<"file read error"<<endl;
							}
							else{
								sumi=sumi+1;
                        //  	dphi[n+N*(m+M*l)][0]=treal[l+dlocal_n0*(m+M*n)]; //[n][m][l]; 
								dphi[n+N*(m+M*l)][0]=treal[n+N*(m+M*l)]; //[n][m][l]; 
						    
							}
                        dphi[n+N*(m+M*l)][1]=0; 
						}
                     }
                 }
	 }else{
                    for(int n=0; n<N; n++){
                    for(int m=0; m<M; m++){
                    for(int l=0; l<dlocal_n0+1; l++){
						 // if(fread(&treal[l+dlocal_n0*(m+M*n)],sizeof(float),1,inf_field)!=1){  // //[n][m][l]
							if(fread(&treal[n+N*(m+M*l)],sizeof(float),1,inf_field)!=1){  // //[n][m][l]
							cout<<"file read error"<<endl;
							}
							else{
								sumi=sumi+1;
                       	// ddphi[n+N*(m+M*l)][0]=treal[l+dlocal_n0*(m+M*n)]; //[n][m][l]; 
				   		ddphi[n+N*(m+M*l)]=treal[n+N*(m+M*l)]; //[n][m][l]; 
						    
							}
                        }
                     }
                 }

	                for(int n=0; n<N; n++){
                    for(int m=0; m<M; m++){
                    for(int l=0; l<dlocal_n0; l++){
	                dphi[n+N*(m+M*l)][0]=ddphi[n+N*(m+M*(l+1))];
	                dphi[n+N*(m+M*l)][1]=0;
	 						}
                     }
                 }
	 
	 
	 
	 }
     cout<<"sumi="<<sumi<<endl;
     fclose(inf_field); 
//	 delete [] treal;
	 free(treal);
//free(buf);
}


static void READPHIVTIC(fftw_complex *dcon, int dlocal_n0,int dmyid){
     
     char filename[20];
	 sprintf(filename, "%s%d%s%d%s", "con_",111001,"_",dmyid,".vti");

	 float *treal,*ddphi;
	 treal=(float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+1)); //[N][M][L];
	 ddphi=(float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+1));


     FILE * inf_field;
     inf_field=fopen(filename,"r");
	 if(inf_field==NULL){
	   printf("can't open this file\n");
		   exit(0);
	 }
     char bufc[100]; //�洢�ַ�
     int i = 1;
     while (i<13) {
       fgets(bufc, 100, inf_field) ;
       printf("%d %s", i++, bufc);
     }  
     char ch;
     ch=fgetc(inf_field);
     cout<<"ch="<<ch<<endl;

     int bytes;
     fread(&bytes,sizeof(int),1,inf_field);
//     cout<<"dmyid="<<dmyid<<" "<<bytes<<endl;
	 cout<<bytes<<endl;  
	 int sumi=0;
     if(dmyid==0){           
	                for(int n=0; n<N; n++){
                    for(int m=0; m<M; m++){
                    for(int l=0; l<dlocal_n0; l++){
						 // if(fread(&treal[l+dlocal_n0*(m+M*n)],sizeof(float),1,inf_field)!=1){  // //[n][m][l]
							if(fread(&treal[n+N*(m+M*l)],sizeof(float),1,inf_field)!=1){  // //[n][m][l]
							cout<<"file read error"<<endl;
							}
							else{
								sumi=sumi+1;
                        //  	dcon[n+N*(m+M*l)][0]=treal[l+dlocal_n0*(m+M*n)]; //[n][m][l]; 
								dcon[n+N*(m+M*l)][0]=treal[n+N*(m+M*l)]; //[n][m][l]; 
						    
							}
                        dcon[n+N*(m+M*l)][1]=0.0; 
						}
                     }
                 }
	 }else{
                    for(int n=0; n<N; n++){
                    for(int m=0; m<M; m++){
                    for(int l=0; l<dlocal_n0+1; l++){
						 // if(fread(&treal[l+dlocal_n0*(m+M*n)],sizeof(float),1,inf_field)!=1){  // //[n][m][l]
							if(fread(&treal[n+N*(m+M*l)],sizeof(float),1,inf_field)!=1){  // //[n][m][l]
							cout<<"file read error"<<endl;
							}
							else{
								sumi=sumi+1;
                       	// ddphi[n+N*(m+M*l)][0]=treal[l+dlocal_n0*(m+M*n)]; //[n][m][l]; 
				   		ddphi[n+N*(m+M*l)]=treal[n+N*(m+M*l)]; //[n][m][l]; 
						    
							}
                        }
                     }
                 }

	                for(int n=0; n<N; n++){
                    for(int m=0; m<M; m++){
                    for(int l=0; l<dlocal_n0; l++){
	                dcon[n+N*(m+M*l)][0]=ddphi[n+N*(m+M*(l+1))];
	                dcon[n+N*(m+M*l)][1]=0.0;
	 						}
                     }
                 }
	 
	 
	 
	 }
     cout<<"sumi="<<sumi<<endl;
     fclose(inf_field); 
//	 delete [] treal;
	 free(treal);
//free(buf);
}


static void SMOOTHCON(fftw_complex *dcon,int dlocal_n0,int dlocal_0_start,int dmyid,int dnumprocs) {
 int nsum,ix,jy,kz,lspace=1, iteration=50;
// int nsum,ix,jy,kz,lspace=2, iteration=20;
 float psum;
 float *buf,*rbuf,*temp;
 float *phi, *dconave; 
 phi= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 dconave= (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0));
 buf = (float *) fftw_malloc(sizeof(float) * N*M*(dlocal_n0+2*lspace));
 MPI_Status stat;
 int left=dmyid-1,right=dmyid+1;
 int head=0,tail= dnumprocs-1;
 if (dmyid>0) {left=dmyid-1;}else{left=MPI_PROC_NULL;}
 if (dmyid<dnumprocs-1){right=dmyid+1;}else{right=MPI_PROC_NULL;} 

//"calculate intial value of dconave"<<endl;
   for (int l=0; l <dlocal_n0; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++) {
       phi[n+N*(m+M*l)]=dcon[n+N*(m+M*l)][0];
   dconave[n+N*(m+M*l)]=dcon[n+N*(m+M*l)][0]; //abs(phi[n+N*(m+M*l)]-u0)
   }

//"input buf of dconave"<<endl; 
   for(int l=0; l <dlocal_n0+2*lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=con0;   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
       buf[n+N*(m+M*l)]=dconave[n+N*(m+M*(l-lspace))];   }

// "smoothening"
   for(int ll=1;ll<iteration+1;ll++){
     if(dmyid==head) {MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,tail,40,&buf[0],N*M*lspace,MPI_FLOAT,tail,40,MPI_COMM_WORLD,&stat);//l1---l0
     }
     if(dmyid==tail) {MPI_Sendrecv(&buf[(dlocal_n0)*N*M],N*M*lspace,MPI_FLOAT,head,40,&buf[(dlocal_n0+lspace)*N*M],N*M*lspace,MPI_FLOAT,head,40,MPI_COMM_WORLD,&stat);//dlocal_n0
     }
     MPI_Sendrecv(&buf[N*M*(dlocal_n0)],N*M*lspace,MPI_FLOAT,right,20,&buf[0],N*M*lspace,MPI_FLOAT,left,20,MPI_COMM_WORLD,&stat);//lnp---l0
     MPI_Sendrecv(&buf[N*M*lspace],N*M*lspace,MPI_FLOAT,left,20,&buf[N*M*(dlocal_n0+lspace)],N*M*lspace,MPI_FLOAT,right,20,MPI_COMM_WORLD,&stat); //l1---lnp+1
     for (int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++){
        psum=0.;nsum=0;int n=0;
        for(int lx=-lspace;lx<lspace+1;lx++) for(int ly=-lspace;ly<lspace+1;ly++) {
          ix=l+lx;
          jy=m+ly; 
          kz=0;
          if (jy > (M-1)) {jy=jy-M+1;}
   //       if (kz > (N-1)) {kz=kz-N+1;}
          if (jy < 0) {jy=jy+M-1;}
   //       if (kz < 0) {kz=kz+N-1;}
//		  if(ix<2){ix=ix+L-1;}
//        if(ix>(L+1)){ix=ix-L+1;}
          psum=psum+buf[kz+N*(jy+M*ix)];
          nsum=nsum+1; 
	      }
        buf[n+N*(m+M*l)]=psum/nsum; 
     }
   }
   for(int l=lspace; l <dlocal_n0+lspace; l ++) for (int m=0; m <M; m ++) for(int n=0; n <N; n ++){
      dconave[n+N*(m+M*(l-lspace))]=buf[n+N*(m+M*l)];
   }
   
    for(int l=0; l <dlocal_n0; l ++){
	   for (int m=0; m <M; m ++){
		for (int n=0; n <N; n ++) {
        dcon[n+N*(m+M*l)][0]=dconave[n+N*(m+M*l)];
        dcon[n+N*(m+M*l)][1]=1.0E-10;
	   }}}
   

}	




