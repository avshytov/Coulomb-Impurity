#include <Python.h>
#include <math.h>
#include <stdio.h>
#ifdef NEED_MALLOC_H
#include <malloc.h>
#endif 

static PyObject * error_obj;

static
double G(double k, int *mvals, int Nm, double r1, double r2) {
       double s = 0.0;
       double ra = r1, rb = r2; 
   
       if (r2 < r1) {
           ra = r2; 
	   rb = r1; 
       }
   
       double xa = k * ra, xb = k * rb; 
       for (int i_m = 0; i_m < Nm; i_m ++) {
            int m = mvals[i_m]; 
	    ///fprintf (stderr, "m = %d  k = %g xa = %g xb = %g\n", m, k, xa, xb); 
	    double j1a = jn(m,     xa); 
	    double j2a = jn(m + 1, xa); 
	    double j1b = jn(m,     xb); 
	    double j2b = jn(m + 1, xb); 
	    double y1b = yn(m,     xb); 
	    double y2b = yn(m + 1, xb);
	    s += (j1a * j1a + j2a * j2a) * (j1b * y1b + j2b * y2b); 
       }
       s *= k * k; 
       return s; 
}

static double wg[] = {
    0.066671344308688137593568809893332, 
    0.149451349150580593145776339657697,
    0.219086362515982043995534934228163,
    0.269266719309996355091226921569469,
    0.295524224714752870173892994651338
};

static
double xgk[] = {
    0.995657163025808080735527280689003, 
    0.973906528517171720077964012084452,
    0.930157491355708226001207180059508,
    0.865063366688984510732096688423493,
    0.780817726586416897063717578345042,
    0.679409568299024406234327365114874,
    0.562757134668604683339000099272694,
    0.433395394129247190799265943165784,
    0.294392862701460198131126603103866,
    0.148874338981631210884826001129720,
    0.000000000000000000000000000000000
};

static
double wgk[] = {
    0.011694638867371874278064396062192,
    0.032558162307964727478818972459390,
    0.054755896574351996031381300244580,
    0.075039674810919952767043140916190,
    0.093125454583697605535065465083366,
    0.109387158802297641899210590325805,
    0.123491976262065851077958109831074,
    0.134709217311473325928054001771707,
    0.142775938577060080797094273138717,
    0.147739104901338491374841515972068,
    0.149445554002916905664936468389821
};

static 
inline double min(double x, double y) {
       return (x < y) ? x : y; 
}

static 
inline double max(double x, double y) {
       return (x > y) ? x : y; 
}

static
double Q_intra_sub(double kmin, double kmax, int *mvals, int Nm, double r1, double r2, 
		   double *abserr, double *resabs, double *resasc) {
    double kc = 0.5 * (kmin + kmax); 
    double hlength = 0.5 * (kmax - kmin); 
    double dlength = fabs(hlength); 
    double resg = 0.0; 
    double fc = G(kc, mvals, Nm, r1, r2); 
    double resk = wgk[10] * fc; 
    double reskh = 0.0;
    double result = 0.0; 
   
    double fv1[10], fv2[10]; 
    double eps_mach = 1e-10, uflow = 1e-10; 
   
    *resabs = fabs(resk); 
    
    for (int j = 0; j < 5; j ++) {
         int jtw = 2 * j + 1; 
         double absc = hlength * xgk[jtw];
         double fval1 = G(kc - absc, mvals, Nm, r1, r2);
         double fval2 = G(kc + absc, mvals, Nm, r1, r2);
         fv1[jtw] = fval1; 
         fv2[jtw] = fval2; 
         double fsum = fval1 + fval2; 
         resg += wg[j] * fsum; 
         resk += wgk[jtw] * fsum; 
         *resabs += wgk[jtw] * (fabs(fval1) + fabs(fval2)); 
    }
    for (int j = 0; j < 5; j ++) {
         int jtwm1 = 2 * j; 
         double absc = hlength * xgk[jtwm1]; 
         double fval1 = G(kc - absc, mvals, Nm, r1, r2);
         double fval2 = G(kc + absc, mvals, Nm, r1, r2);
         fv1[jtwm1] = fval1;
         fv2[jtwm1] = fval2;
         double fsum = fval1 + fval2; 
         resk += wgk[jtwm1] * fsum; 
         *resabs += wgk[jtwm1] * (fabs(fval1) + fabs(fval2));          
    }
   
    reskh = resk * 0.5; 
    *resasc = wgk[10] * fabs(fc - reskh);
    for (int j = 0; j < 10; j ++)
         *resasc += wgk[j] * (fabs(fv1[j] - reskh) + fabs(fv2[j] - reskh)); 
    
    result = resk * hlength; 
    *resabs *= dlength; 
    *resasc *= dlength; 
    
    *abserr = fabs((resk - resg) * hlength); 
    if (*resasc != 0.0 &&  *abserr != 0.0) {
        *abserr = *resasc * min(1.0, pow(200.0 * (*abserr) / (*resasc), 1.5));
    }
    if (*resabs > uflow / (50 * eps_mach))
        *abserr = max(*abserr, 50 * eps_mach * (*resabs) ); 

    return result; 
}


static
double Q_intra(double kF, int *mvals, int Nm, double r1, double r2) {
       double rmax = (r1 > r2) ? r1 : r2; 
       double dxi0 = 4.0 * M_PI; 
       double dk = dxi0 / rmax; 
       int Nk = floor(kF / dk) + 1; 
       double s = 0.0; 
       dk = kF / Nk;
       for (int i = 0; i < Nk; i ++) {
            double kmin = dk * i; 
	    double kmax = kmin + dk; 
	    double abserr = 0.0, resabs = 0.0, resasc = 0.0; 
	    s += Q_intra_sub(kmin, kmax, mvals, Nm, r1, r2, 
			     &abserr, &resabs, &resasc); 
	    //fprintf (stderr, "%g %g %g\n", abserr, resabs, resasc); 
       }
       return s / 8.0 / M_PI; 
}

static 
double Q_intra_old(double kF, int *mvals, int Nm, double r1, double r2) {
       double rmax = (r1 > r2) ? r1 : r2; 
       double dxi0 = 0.1; 
       double ximax = kF * rmax;  
       int Nk = 50 + floor(ximax / dxi0);
       //fprintf (stderr, "Nk = %d\n", Nk); 
       // We do not include the k=0 point, 
       // due to k^2 * log(k) behaviour at k->0
       double s = 0.0; 
       double dk = kF / Nk; 
       for (int i = 1; i < Nk; i ++) {
            double k_i = i * dk; 
	    s += 4.0 * G(k_i, mvals, Nm, r1, r2); 
	    s += 2.0 * G(k_i + 0.5 * dk, mvals, Nm, r1, r2); 
       }
       s += 1.0 * G(kF, mvals, Nm, r1, r2); 
       s *= dk / 6.0; 
       s *= 1.0 / 8.0 / M_PI; 
       return s; 
}

static
void Q_intra_quad(double kF, int *mvals, int Nm, double ri, 
		  double *xi, double *w, int nw,
		  double r1, double r2, 
		  double *S0, double *S1) 
{
     *S0 = *S1 = 0.0;
     double rc = (r1 + r2) / 2.0; 
     double dr = (r2 - r1); 
     for (int i = 0; i < nw; i ++) {
          double r  = rc + xi[i] * dr / 2.0;
	  double Qr = Q_intra(kF, mvals, Nm, ri, r);
	  *S0 += Qr * r * w[i]; 
	  *S1 += Qr * r * (r - rc) * w[i]; 
     }
     *S0 *= dr / 2.0; 
     *S1 *= dr / 2.0; 
}

static 
void Q_intra_int(double kF, int *mvals, int Nm, double ri,  
		 double r1, double r2, int nj, 
		 double *S0, double *S1) {
       /* if (nj == 1) {
	   double xi[] = {0.0}; 
	   double w[]  = {2.0}; 
	   Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w)/sizeof(*w), r1, r2, S0, S1); 
       } else if (nj == 2) {
           double xi[] = {-0.5773502691896258, 0.5773502691896258}; 
	   double w [] = {1.0,                 1.0}; 
	   Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w)/sizeof(*w), r1, r2, S0, S1);
       } else */if (nj <= 3) {
           double xi[] = {-0.7745966692414834, 0.0,               0.7745966692414834}; 
	   double w [] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
	   Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w)/sizeof(*w), r1, r2, S0, S1);
       } else if (nj == 4) {
           double xi[] = {0.8611363115940526, 0.3399810435848562, -0.3399810435848562, -0.8611363115940526}; 
	   double w [] = {0.347854845137, 0.652145154863, 0.652145154863, 0.347854845137};
	   Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w) /sizeof(*w), r1, r2, S0, S1);
       }  else if (nj <= 7) {
           double xi[] = {0.906179845939, 0.538469310106, 0.0,          -0.538469310106, -0.906179845939}; 
	   double w [] = {0.236926885056, 0.478628670499, 0.56888888888, 0.478628670499, 0.236926885056};
	   Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w) / sizeof(*w), r1, r2, S0, S1);
       } else if (nj <= 13) { // OK
	   double xi[] = {
	        0.978228658146, 0.887062599768, 0.730152005574,  
	        0.519096129207, 0.26954315595, 
	        0.0, 
	        -0.26954315595, -0.519096129207, 
	        -0.730152005574,  -0.887062599768, -0.978228658146
	   }; 
	   double w[] = {
	        0.0556685671162, 0.125580369465, 0.186290210928, 
	        0.23319376459, 0.26280454451, 
	        0.272925086778, 
	        0.26280454451,  0.233193764592, 
	        0.186290210928, 0.125580369465, 0.055668567116
	   }; 
	   Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w)/sizeof(*w), 
			r1, r2, S0, S1);
       } else {
	   double xi[] = {
	        0.98799251802,  0.937273392401, 0.84820658341, 
	        0.72441773136,  0.570972172609, 0.39415134707, 
	        0.201194093997, 
	        0.0, 
	        -0.201194093997, 
	        -0.394151347078, -0.570972172609, -0.72441773136, 
	        -0.84820658341,  -0.937273392401, -0.98799251802
	   };
	   double w[]= {
	        0.0307532419961, 0.0703660474881, 0.107159220467, 
	        0.139570677926, 0.166269205817, 0.186161000016, 
	        0.198431485327, 0.202578241926, 0.198431485327, 
	        0.186161000016, 0.166269205817, 0.139570677926, 
	        0.107159220467, 0.0703660474881, 0.0307532419961
	   }; 
           Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w)/sizeof(*w), 
			r1, r2, S0, S1);                   
       } /*else if (nj <= 21) {
      	   Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w)/sizeof(*w), 
			r1, r2, S0, S1);        
       }*/ /*else if (nj <= 25) {
	   Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w)/sizeof(*w), 
			r1, r2, S0, S1);        
       } else if (nj <= 31) {
	   Q_intra_quad(kF, mvals, Nm, ri, xi, w, sizeof(w)/sizeof(*w), 
			r1, r2, S0, S1);        
       } else { 
	       fprintf (stderr, "normal quad: %d %g %g %g\n", nj, ri, r1, r2); 
               *S0 = 0.0; 
	       *S1 = 0.0;
               double rc = (r1 + r2) / 2.0; 
               double dr = (r2 - r1); 
	       double drk = dr / nj; 
	       //fprintf (stderr, "nj = %d\n", nj); 
	       
	       for (int k = 1; k < nj + 1; k ++) {
		    double rk = r1 + dr / (nj + 1) * k; 
		    double Qk = Q_intra(kF, mvals, Nm, ri, rk); 
	            *S0 += Qk * rk; 
		    *S1 += Qk * rk * (rk - rc); 
	       }
	       *S0 *= drk; 
	       *S1 *= drk;       
       }*/
}


static
void do_kernel_m_intra(double **Q, double *r, int Nr, 
		       int *mvals, int Nm, 
		       double kF) {
   
     double dr0 = 0.1;
     double S0, S1; 
   
     if (kF > 1.0)
         dr0 = dr0 / kF ; 
   
     for (int i = 0; i < Nr; i ++)
          for (int j = 0; j < Nr; j ++) 
               Q[i][j] = 0.0; 
     
     for (int i = 0; i < Nr; i ++) {
	  for (int j = 0; j < Nr - 1; j ++) {
	       //fprintf (stderr, "i = %d j = %d\n", i, j); 
	       double r1 = r[j]; 
	       double r2 = r[j + 1]; 
	       //double rc = (r1 + r2) / 2.0; 
	       double dr = r2 - r1; 
	       int nj = floor(dr / dr0)  + 1; 
	       S0 = 0.0; 
	       S1 = 0.0;
	       
	       Q_intra_int(kF, mvals, Nm, r[i], r1, r2, nj, &S0, &S1); 
	       //double drk = dr / nj; 
	       //fprintf (stderr, "nj = %d\n", nj); 
	       
	       //for (int k = 1; k < nj + 1; k ++) {
	       //	    double rk = r1 + (r2 - r1) / (nj + 1) * k; 
	       //	    double Qk = Q_intra(kF, mvals, Nm, r[i], rk); 
	       //     S0 += Qk * rk; 
	       //	    S1 += Qk * rk * (rk - rc); 
	       //}
	       //S0 *= drk; 
	       //S1 *= drk;
	       /**/
	       //Q[i][j] += Q_intra(kF, mvals, Nm, r[i], r[j]) * rc * dr; 
	       Q[i][j]     += 0.5 * S0 - S1 / dr; 
	       Q[i][j + 1] += 0.5 * S0 + S1 / dr; 
	  }
	  //Q[i][0] += Q_intra(kF, mvals, Nm, r[i], 0.0) * r[0] * r[0] / 2.0; 
	  S0 = 0.0; 
	  S1 = 0.0; 
	  Q_intra_int(kF, mvals, Nm, r[i], 0.0, r[0], 5, &S0, &S1);
	  Q[i][0] += S0; 
	  double sq = 0.0; 
	  for (int j = 0; j < Nr; j ++)
	       sq += Q[i][j]; 
	  sq *= 8.0 * M_PI * M_PI; 
	  fprintf (stderr, "intra: %d %g\n", i, sq); 
     }
     for (int i = 0; i < Nr; i ++)
          for (int j = 0; j < Nr; j ++)
               Q[i][j] *= 4.0 * M_PI; 
}

static 
PyObject* kernel_m_intra(PyObject * self, PyObject * args) {
          double kF; 
          double *r; 
          int *mvals; 
          double **Q; 
          PyObject *Q_obj; 
          int Nr, Nm; 
   
          PyObject *mvals_obj, *r_obj; 
          if (!PyArg_ParseTuple(args, "OOd", &r_obj, &mvals_obj, &kF)) {
	       PyErr_SetString(error_obj, "invalid arguments to kernel_m_intra"); 
	       return NULL; 
	  }
          if (!PyList_Check(mvals_obj)) {
	       PyErr_SetString(error_obj, "mvals is not a list"); 
	       return NULL; 
	  }
          if (!PyList_Check(r_obj)) {
	       PyErr_SetString(error_obj, "r is not a list"); 
	       return NULL; 
	  }
          Nr = PyList_Size(r_obj); 
          Nm = PyList_Size(mvals_obj); 
          r = PyMem_New(double, Nr); 
          mvals = PyMem_New(int, Nm); 
          Q = PyMem_New(double *, Nr); 
          for (int i = 0; i < Nr; i ++) {
	       PyObject *val = PyList_GetItem(r_obj, i); 
	       double r_i; 
	       PyArg_Parse(val, "d", &r_i);
	       r[i] = r_i; 
	       Q[i] = PyMem_New(double, Nr); 
	  }
          for (int i = 0; i < Nm; i ++) {
	       PyObject *val = PyList_GetItem(mvals_obj, i); 
	       double m_i; 
	       PyArg_Parse(val, "d", &m_i); 
	       mvals[i] = (int)floor(m_i);
   	       if (fabs(m_i - floor(m_i)) > 1e-2) {
	           fprintf (stderr, "Warning: only integer values of m are now supported! %g %d\n", m_i, mvals[i]); 
	       }
	  }
   
          do_kernel_m_intra(Q, r, Nr, mvals, Nm, kF); 
   
          Q_obj = PyList_New(Nr); 
          for (int i = 0; i < Nr; i ++) {
	       PyObject *Q_i = PyList_New(Nr); 
	       for (int j = 0; j < Nr; j ++) {
	           PyList_SetItem(Q_i, j, Py_BuildValue("d", Q[i][j])); 
	       }
	       PyList_SetItem(Q_obj, i, Q_i); 
	  }
   
          for (int i = 0; i < Nr; i ++) {
	       PyMem_Free(Q[i]); 
	  }
          PyMem_Free(Q); 
          PyMem_Free(mvals); 
          PyMem_Free(r); 
          return Py_BuildValue("N", Q_obj);    
}

static 
PyObject* qm_intra(PyObject * self, PyObject * args) {
          double kF; 
          double r1, r2; 
          int *mvals; 
          double Qvalue; 
          //PyObject *Q_obj; 
          int Nm; 
   
          PyObject *mvals_obj, *r_obj; 
          if (!PyArg_ParseTuple(args, "dOdd", &kF, &mvals_obj, &r1, &r2)) {
	       PyErr_SetString(error_obj, "invalid arguments to qm_intra"); 
	       return NULL; 
	  }
          if (!PyList_Check(mvals_obj)) {
	       PyErr_SetString(error_obj, "mvals is not a list"); 
	       return NULL; 
	  }
          Nm = PyList_Size(mvals_obj); 
          mvals = PyMem_New(int, Nm); 
          for (int i = 0; i < Nm; i ++) {
	       PyObject *val = PyList_GetItem(mvals_obj, i); 
	       double m_i; 
	       PyArg_Parse(val, "d", &m_i); 
	       mvals[i] = (int)floor(m_i);
   	       if (fabs(m_i - floor(m_i)) > 1e-2) {
	           fprintf (stderr, "Warning: only integer values of m are now supported! %g %d\n", m_i, mvals[i]); 
	       }
	  }
   
          Qvalue = Q_intra(kF, mvals, Nm, r1, r2);    
          PyMem_Free(mvals); 
          return Py_BuildValue("d", Qvalue);    
}

static struct PyMethodDef module_methods[] = {
   {"kernel_m_intra", kernel_m_intra, METH_VARARGS, "Calculates RPA intraband contribution"}, 
   {"qm_intra",       qm_intra,       METH_VARARGS, "Calculates RPA intraband contribution"}, 
   {NULL, NULL, 0, NULL},
};

void init_rpam() {
     PyObject *module, *module_dict; 
     PyMethodDef *method_def; 
     module = Py_InitModule ("_rpam", module_methods); 
     module_dict = PyModule_GetDict (module); 
     error_obj = Py_BuildValue("s", "_rpam.error"); 
     if (PyErr_Occurred()) {
	            Py_FatalError ("cannot initialize the module _rpam"); 
    }
}

/*
int main(int argc, char **argv) {
    double kF = 3.0; 
    double rmin = 0.01, rmax = 50.0; 
    int Nr = 500; 
    int Mmax = 31; 
    int *mvals = (int *) malloc(sizeof(*mvals) * Mmax); 
    double *r = (double *) malloc(sizeof(*r) * Nr);
    double **Q = (double **) malloc(sizeof(*Q) * Nr); 
    for (int i = 0; i < Nr; i ++) {
         Q[i] = (double *) malloc(sizeof(*Q[i]) * Nr); 
    }
    
    for (int i = 0; i < Nr; i ++) {
        r[i] = rmin + (rmax - rmin) / (Nr - 1) * i; 
    }
   
    for (int m = 0; m < Mmax; m ++) 
         mvals[m] = m; 
   kernel_m_intra(Q, r, Nr, mvals, Mmax, kF); 
}
*/
