static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
(maxarg1) : (maxarg2))
void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
void free_vector(float *v,long nl, long nh);
void free_ivector(int *v,long nl, long nh);
void free_matrix(float **m,long nrl,long nrh, long ncl, long nch);
void powell(float p[],float **xi,int n,float ftol,int *iter,float *fret,
	    float (*func)(float []));
void linmin(float p[], float xi[], int n, float *fret, 
	    float (*func)(float []));
float brent(float ax, float bx, float cx, float (*f)(float), float tol,
              float *xmin);
float f1dim(float x);
void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
              float *fc, float (*func)(float));
char *ggets(char *s);
void moment(float data[], int n, float *ave, float *adev, float *sdev,
	    float *var, float *skew, float *curt);
