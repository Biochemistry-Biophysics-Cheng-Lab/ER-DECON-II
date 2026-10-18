


int         core       (  double   *Imn,
						  double   *TF,
						  double    *TFc,
						  double   *X,
						  long     *Sz,
						  double   lamratio,
						  double   delfac,
						  double   lamf,
						  double   lampos,
						  double   lampc,
						  double   lamsmooth,
						  double   alpha,
						  long     NIter,
						  long     NIerT,
						  double   tol,
						  double   omega,
			                          int      nthread,
			                          void     (*res_reporter)(long iiters, long iitert, double res, void*),
			                          void*    res_reporter_closure
);
