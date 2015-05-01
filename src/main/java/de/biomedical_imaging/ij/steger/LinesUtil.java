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
    aint with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. */

/* 	Changes Made by R. Balasubramanian for incorporating the the detect lines code to incorporate
   	within GRASP (May 10th 1999) */

/*	Translated into an ImageJ java plugin by Thorsten Wagner (Dez. 2014) */

package de.biomedical_imaging.ij.steger;

public class LinesUtil {
	
	public static final int DERIV_R  = 1;  /* Derivative in row direction */
	public static final int DERIV_C  = 2;  /* Derivative in column direction */
	public static final int DERIV_RR = 3;  /* Second derivative in row direction */
	public static final int DERIV_RC = 4;  /* Second derivative in row and column direction */
	public static final int DERIV_CC = 5;  /* Second derivative in column direction */

	public static final int MODE_LIGHT = 1;  /* Extract bright lines */
	public static final int MODE_DARK  = 2 ; /* Extract dark lines */

	public static final int INITIAL_SIZE  = 100;
	public static final int REALLOC_FACTOR = 2; 

	public static final double MAX_SIZE_MASK_0 =  3.09023230616781;    /* Size for Gaussian mask */
	public static final double MAX_SIZE_MASK_1 = 3.46087178201605;   /* Size for 1st derivative mask */
	public static final double MAX_SIZE_MASK_2 = 3.82922419517181;    /* Size for 2nd derivative mask */

	public static int MASK_SIZE(double MAX, double sigma) {
		return (int)Math.ceil(MAX*sigma); /* Maximum mask index */
	}
	
	public static final String ERR_NOMEM = "Out of memory";
	public static final String ERR_FNF   = "Could not open file";
	public static final String ERR_NOPGM = "Not a valid pgm file:";
	public static final String ERR_SNS   = "Sigma not specified";
	public static final String ERR_SOR   = "Sigma out of range:";
	public static final String ERR_LNS   = "Low not specified";
	public static final String ERR_LOR   = "Low out of range:";
	public static final String ERR_HNS   = "High not specified";
	public static final String ERR_HOR   = "High out of range:";
	public static final String ERR_LGH   = "Low > High";
	public static final String ERR_CNW   = "Line position correction impossible without line width";
	public static final String ERR_INP   = "Include-image option requires PostScript output";
	public static final String ERR_UKO   = "Unknown option:";
	public static final String ERR_TMF   = "Too many files specified:";
	
	
	/* Translate row and column coordinates of an image into an index into its
	   one-dimensional array. */
	public static int LINCOOR(int row, int col, int width) {
		return row*width+col;
	}
	
	/* Mirror the row coordinate at the borders of the image; height must be a
	   defined variable in the calling function containing the image height. */
	public static int BR(int row, int height) {
		return ((row) < 0 ? -(row) : (row) >= height ? height - (row) + height - 2 : (row));
	}
	
	/* Mirror the column coordinate at the borders of the image; width must be a
	   defined variable in the calling function containing the image width. */
	public static int BC(int col, int width) { 
		return ((col) < 0 ? -(col) :  (col) >= width ? width - (col) + width - 2 : (col));
	}
	               
	public enum contour_class {
		  cont_no_junc,    /* no end point is a junction */
		  cont_start_junc, /* only the start point of the line is a junction */
		  cont_end_junc,   /* only the end point of the line is a junction */
		  cont_both_junc,  /* both end points of the line are junctions */
		  cont_closed;      /* the contour is closed */
	}
	               

}
