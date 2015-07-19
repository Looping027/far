/* The comments in this file can be compiled into a
   documentation file using doxygen (http://www.doxygen.org) */
/*! \mainpage far library C source code

 \section license license

This software is distributed under the terms of the 
LESSER GNU GENERAL PUBLIC LICENSE version 2.1. 
The terms of this license can be obtained via WWW at
http://www.gnu.org/copyleft/lesser.html

``Share and Enjoy.''
 
  \section whatis what is it ?

Those programs are part of R library named far. The far library
need some C code to work faster.

*/

/*!
  \file   Far.c
  \brief  Estimation of FAR time series

  Nonparametric estimation of time series
  Those programs implement different
  nonparametric methods of estimation for
  time series: scalar and functional kernel.

  \author Julien Damon <julien.damon@free.fr>
  \date   2001-02-27
*/
#include <R.h>

/*! Multidimentionnal gaussian kernel */
/*!
  Multidimentionnal gaussian kernel K(u) without the
  normalization constant.

  \param serie_x       vector u
  \param n             length of serie_x
  \param h             width of the window (h)
  \param result        result
  \return              the exponential of minus squared norm of u 
                       divided by squared h
  \version 0.5-1
  \date 2001-02-27
*/
/* Multidimentionnal gaussian kernel */
void kersca (double *serie_x, int *n, double *h, double *result)
{
  int i,nn=*n;
  result[0] = 0.0;
  for (i=0; i<nn; i++)
  {
    result[0] += serie_x[i] * serie_x[i] ;
  }
  result[0] = exp(- result[0] / 2 / *h / *h );
}

/*! Kernel estimation */
/*!
  Kernel estimation of a functional autoregressive process using 
  a given window width

  \param serie_x        past values of the serie (1 --> n-1)
  \param serie_y        past values of the serie (2 --> n)
  \param serie_x0       actual value of the functional serie 
                        (the starting point)
  \param h              window width
  \param n              number of observations in the functional serie serie_x
  \param p              number of point for each observation of serie_x 
                        (length of *serie_x0)
  \param result         vector containing the estimation K(serie_x0)
  \return               the functional kernel estimation of *serie_x1 obtained 
                        by the formula serie_x1=K(serie_x0) where K depends on 
                        the past values of the serie
*/
/* Kernel estimation */
void prevkerfon (double *serie_x, double *serie_y, double *serie_x0,
                 double *h, int *n, int *p, double *result)
{
  int i, j, pp=*p;
  double *tmpkval, ttmpkval, kval;
  double *vide = Calloc((size_t) pp, double);
  ttmpkval = 0.0 ;
  tmpkval = &ttmpkval ;
  kval = 0.0;
  for (j= 0; j< pp; j++)
  {
    result[j] = 0.0;
  }
  for (i = 0; i < (*n - 1); i++)
  {
    for (j= 0; j< pp; j++)
    {
      vide[j] = serie_x[(i * pp) + j] - serie_x0[j];
    }
    kersca(vide, p, h, tmpkval);
    kval += *tmpkval;
    for (j= 0; j< pp; j++)
    {
      result[j] += serie_y[(i * pp) + j] * *tmpkval;
    }
    ttmpkval = 0.0 ;
  }
  for (j= 0; j< pp; j++)
  {
    result[j] /= kval;
  }
  Free(vide);
}

/*! Cross validation of Functional Kernel */
/*!
  Cross validation of the window width used in the functional kernel 
  estimation of an autoregressive process

  \param serie_x        past values of the serie (1 --> n-1)
  \param serie_y        past values of the serie (2 --> n)
  \param h              window width
  \param n              number of observations in the functional serie serie_x
  \param p              number of point for each observation of serie_x 
                        (length of serie_x0)
  \param r              length of the cross-validation period
  \param result         the criteria of the cross validation
  \result               the sum of the square errors obtain with a width of h
*/
/* Cross validation of Functional Kernel */
void CVkerfon (double *serie_x, double *serie_y, double *h, 
               int *n, int *p, int *r, double *result)
{
  int i, tj, *j, k, pp=*p, nn=*n, rr=*r;
  double *vide_court = Calloc((size_t) pp, double),
         *res        = Calloc((size_t) pp, double),
         *vide_long  = Calloc((size_t) ((nn - rr) * pp), double),
         *vide_long2 = Calloc((size_t) ((nn - rr) * pp), double);
  j = &tj;
  tj = nn - rr ;
  *result = 0.0 ;
  for (k = 0;k < ((nn - rr) * pp);k++)
  {
    vide_long[k] = serie_x[k];
    vide_long2[k] = serie_y[k];
  }
  for (k = 0;k < rr;k++)
  {
    for (i = 0;i < pp; i++)
    {
      vide_court[i] = serie_x[((k + nn - rr) * pp) + i];
    }
    prevkerfon(vide_long, vide_long2, vide_court, h, j, p, res);
    for (i = 0;i < pp; i++)
    {
      *result += (res[i] - serie_y[((k + nn - rr) * pp)
              + i]) * (res[i] - serie_y[((k + nn - rr) * pp) + i]);
    }
  }
  Free(vide_court);Free(res);Free(vide_long);Free(vide_long2);
}
