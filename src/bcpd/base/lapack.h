extern "C" {
    extern int dsyev_ (char *jobz,char *uplo,int *n,double *a,int *lda,double *w,double *work,int *lwork,int *info);
    extern int dposv_ (char *uplo,int *n,int *nrhs,double *A,int *lda,double *B,int *ldb,int *info);
    extern int dpotrs_(char *uplo,int *n,int *nrhs,double *A,int *lda,double *B,int *ldb,int *info);
    extern int dpotrf_(char *uplo,int *n,double *A,int *lda,int *info);
    extern int dpotri_(char *uplo,int *n,double *A,int *lda,int *info);
    extern int dgetrf_(int *m,int *n,double *A,int *lda,int *ipiv,int *info);
}
