/*  detect-lines, extract lines and their width from images.
    Copyright (C) 1996-1998 Carsten Steger

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/* 	Changes Made by R. Balasubramanian for incorporating the the detect lines code to incorporate
   	within GRASP (May 10th 1999) */

/*	Translated into an ImageJ java plugin by Thorsten Wagner (Dez. 2014) */

package de.biomedical_imaging.ij.steger;

public class Correction {
	  double w_est;       /* Total line width extracted from the image */
	  double r_est;       /* Gradient ratio extracted from the image */
	  double w;           /* True line width */
	  double h;           /* True asymmetry */
	  double correction;  /* Line position correction */
	  double w_strong;    /* True width on the side with the stronger gradient */
	  double w_weak;      /* True width on the side with the weaker gradient */
	  boolean   is_valid;    /* Is this table entry valid? */
	  
	  public Correction(double w_est,double r_rest, double w, double h, double correction, double w_strong, double w_weak, boolean is_valid) {
		  this.w_est = w_est;
		  this.r_est = r_rest;
		  this.w = w;
		  this.h = h;
		  this.correction = correction;
		  this.w_strong = w_strong;
		  this.w_weak = w_weak;
		  this.is_valid = is_valid;
	}
}
