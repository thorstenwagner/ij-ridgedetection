/*
 * #%L
 * Ridge Detection plugin for ImageJ
 * %%
 * Copyright (C) 2014 - 2015 Thorsten Wagner (ImageJ java plugin), 1996-1998 Carsten Steger (original C code), 1999 R. Balasubramanian (detect lines code to incorporate within GRASP)
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-2.0.html>.
 * #L%
 */
package de.biomedical_imaging.ij.steger;

import org.apache.commons.lang3.mutable.MutableLong;


public class Convol {
	/* 1/sqrt(2*PI) */
	private final double  SQRT_2_PI_INV = 0.398942280401432677939946059935;
	
	/* Functions to compute the integral, and the 0th and 1st derivative of the
	   Gaussian function 1/(sqrt(2*PI)*sigma)*exp(-0.5*x^2/sigma^2) */

	/* Integral of the Gaussian function */
	private double phi0(double x, double sigma)
	{
	  return Normal.getNormal(x/sigma);
	}
	
	/* The Gaussian function */
	private double phi1(double x, double sigma)
	{
	  double t;
	  t = x/sigma;
	  return SQRT_2_PI_INV/sigma*Math.exp(-0.5*t*t);
	}
	
	/* First derivative of the Gaussian function */
	public double phi2(double x, double sigma)
	{
	  double t;
	  t = x/sigma;
	  return -x*SQRT_2_PI_INV/Math.pow(sigma,3.0)*Math.exp(-0.5*t*t);
	}
	
	/* Functions to compute the one-dimensional convolution masks of the 0th, 1st,
	   and 2nd derivative of the Gaussian kernel for a certain smoothing level
	   given by sigma.  The mask is allocated by the function and given as the
	   return value.  The caller must ensure that this memory is freed.  The
	   output is intended to be used as an array with range [-num:num].  Therefore,
	   the caller should add num to the return value.  Examples for the calling
	   sequence can be found in convolve_gauss.  Examples for the usage of the
	   masks are given in convolve_rows_gauss and convolve_cols_gauss. */

	/* Gaussian smoothing mask */
	
	/* num ist eigentlich pointer - aufrufende Funkion nimmt an, dass num geändert wird. Übergebe es deswegen 
	 * als MutableDouble aus CommonsLang
	 *  */
	public double[] compute_gauss_mask_0(MutableLong num, double sigma)
	{
	  
	  int   i, n;
	  double limit;
	  double[] h;

	  limit = LinesUtil.MASK_SIZE(LinesUtil.MAX_SIZE_MASK_0,sigma); /* Error < 0.001 on each side */
	  n = (int)limit;
	  h = new double[2*n+1];
	  for (i=-n+1;i<=n-1;i++)
	    h[n+i] = phi0(-i+0.5,sigma) - phi0(-i-0.5,sigma);
	  h[0] = 1.0 - phi0(n-0.5,sigma);
	  h[2*n] = phi0(-n+0.5,sigma);
	  num.setValue(n);
	  return h;
	}
	
	/* First derivative of Gaussian smoothing mask */
	/* num ist eigentlich pointer - aufrufende Funkion nimmt an, dass num geändert wird. Übergebe es deswegen 
	 * als MutableDouble aus CommonsLang
	 *  */
	public double[] compute_gauss_mask_1(MutableLong num, double sigma)
	{
	  int   i, n;
	  double limit;
	  double[] h;

	  limit = LinesUtil.MASK_SIZE(LinesUtil.MAX_SIZE_MASK_1,sigma); /* Error < 0.001 on each side */
	  n = (int)limit;
	  h = new double[2*n+1];
	  
	  for (i=-n+1;i<=n-1;i++)
	    h[n+i] = phi1(-i+0.5,sigma) - phi1(-i-0.5,sigma);
	  h[0] = -phi1(n-0.5,sigma);
	  h[2*n] = phi1(-n+0.5,sigma);
	  num.setValue(n);
	  return h;
	}
	
	/* Second derivative of Gaussian smoothing mask */
	/* num ist eigentlich pointer - aufrufende Funkion nimmt an, dass num geändert wird. Übergebe es deswegen 
	 * als MutableDouble aus CommonsLang
	 *  */
	public double[] compute_gauss_mask_2(MutableLong num, double sigma)
	{
	  int   i, n;
	  double limit;
	  double[] h;

	  limit = LinesUtil.MASK_SIZE(LinesUtil.MAX_SIZE_MASK_2,sigma); /* Error < 0.001 on each side */
	  n = (int)limit;
	  h = new double[2*n+1];
	  
	  for (i=-n+1;i<=n-1;i++)
	    h[n+i] = phi2(-i+0.5,sigma) - phi2(-i-0.5,sigma);
	  h[0] = -phi2(n-0.5,sigma);
	  h[2*n] = phi2(-n+0.5,sigma);
	  num.setValue(n);
	  return h;
	}
	
	/* Convolve an image with the derivatives of a Gaussian smoothing kernel.
	   Since all of the masks are separable, this is done in two steps in the
	   function convolve_gauss.  Firstly, the rows of the image are convolved by
	   an appropriate one-dimensional mask in convolve_rows_gauss, yielding an
	   intermediate float-image h.  Then the columns of this image are convolved
	   by another appropriate mask in convolve_cols_gauss to yield the final
	   result k.  At the border of the image the gray values are mirrored. */

	/* Convolve the rows of an image with the derivatives of a Gaussian. */
	private void convolve_rows_gauss(float[] image, double[] mask, int n, float[] h,int width,int height)
	{
	  int      j, r, c, l;
	  double    sum;

	  /* Inner region */
	  for (r=n; r<height-n; r++) {
	    for (c=0; c<width; c++) {
	      l = LinesUtil.LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += (double)(image[(l+j*width)])*mask[(j+n)];
	      h[ l] = (float) sum;
	    }
	  }
	  /* Border regions */
	  for (r=0; r<n; r++) {
	    for (c=0; c<width; c++) {
	      l = LinesUtil.LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += (double)(image[ LinesUtil.LINCOOR(LinesUtil.BR(r+j,height),c,width)])*mask[(j+n)];
	      h[ l] = (float) sum;
	    }
	  }
	  for (r=height-n; r<height; r++) {
	    for (c=0; c<width; c++) {
	      l = LinesUtil.LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += (double)(image[ LinesUtil.LINCOOR(LinesUtil.BR(r+j,height),c,width)])*mask[(j+n)];
	      h[ l] = (float) sum;
	    }
	  }
	}
	
	/* Convolve the columns of an image with the derivatives of a Gaussian. */
	private void convolve_cols_gauss(float[] h, double[] mask,int n, float[] k, int width, int height)
	{
	  int      j, r, c, l;
	  double    sum;

	  /* Inner region */
	  for (r=0; r<height; r++) {
	    for (c=n; c<width-n; c++) {
	      l = LinesUtil.LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += h[ (l+j)]*mask[(j+n)];
	      k[ l] = (float)sum;
	    }
	  }
	  /* Border regions */
	  for (r=0; r<height; r++) {
	    for (c=0; c<n; c++) {
	      l = LinesUtil.LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += h[ LinesUtil.LINCOOR(r,LinesUtil.BC(c+j,width),width)]*mask[(j+n)];
	      k[ l] = (float)sum;
	    }
	  }
	  for (r=0; r<height; r++) {
	    for (c=width-n; c<width; c++) {
	      l = LinesUtil.LINCOOR(r,c,width);
	      sum = 0.0;
	      for (j=-n;j<=n;j++)
	        sum += h[ LinesUtil.LINCOOR(r,LinesUtil.BC(c+j,width),width)]*mask[(j+n)];
	      k[ l] = (float)sum;
	    }
	  }
	}
	
	/* Convolve an image with a derivative of the Gaussian. */
	public void convolve_gauss(float[] image,float[] k,int width,int height,double sigma,int deriv_type)
	{
	  double[]  hr = null, hc = null;
	  double[]  maskr, maskc;
	  MutableLong    nr = new MutableLong(), nc= new MutableLong();
	  float[]   h;

	  h = new float[ (width*height)];

	  switch (deriv_type) {
	    case LinesUtil.DERIV_R:
	      hr = compute_gauss_mask_1(nr,sigma);
	      hc = compute_gauss_mask_0(nc,sigma);
	      break;
	    case LinesUtil.DERIV_C:
	      hr = compute_gauss_mask_0(nr,sigma);
	      hc = compute_gauss_mask_1(nc,sigma);
	      break;
	    case LinesUtil.DERIV_RR:
	      hr = compute_gauss_mask_2(nr,sigma);
	      hc = compute_gauss_mask_0(nc,sigma);
	      break;
	    case LinesUtil.DERIV_RC:
	      hr = compute_gauss_mask_1(nr,sigma);
	      hc = compute_gauss_mask_1(nc,sigma);
	      break;
	    case LinesUtil.DERIV_CC:
	      hr = compute_gauss_mask_0(nr,sigma);
	      hc = compute_gauss_mask_2(nc,sigma);
	      break;
	  }

	  maskr = hr;// + nr; Wird ersetzt in den eigentlichen Funktionen, indem ich z.B. in convolve_rows_gauss immer beim Zugriff auf mask n dazuaddiere
	  maskc = hc;// + nc;

	  convolve_rows_gauss(image,maskr,nr.intValue(),h,width,height);
	  convolve_cols_gauss(h,maskc,nc.intValue(),k,width,height);

	}

}
