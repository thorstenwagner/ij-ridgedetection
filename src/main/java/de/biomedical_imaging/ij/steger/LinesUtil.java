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

// TODO: Auto-generated Javadoc
/**
 * The Class LinesUtil.
 */
public class LinesUtil {

	/** The Constant DERIV_R. */
	public static final int DERIV_R = 1; /* Derivative in row direction */

	/** The Constant DERIV_C. */
	public static final int DERIV_C = 2; /* Derivative in column direction */

	/** The Constant DERIV_RR. */
	public static final int DERIV_RR = 3; /* Second derivative in row direction */

	/** The Constant DERIV_RC. */
	public static final int DERIV_RC = 4; /* Second derivative in row and column direction */

	/** The Constant DERIV_CC. */
	public static final int DERIV_CC = 5; /* Second derivative in column direction */

	/** The Constant MODE_LIGHT. */
	public static final int MODE_LIGHT = 1; /* Extract bright lines */

	/** The Constant MODE_DARK. */
	public static final int MODE_DARK = 2; /* Extract dark lines */

	/** The Constant INITIAL_SIZE. */
	public static final int INITIAL_SIZE = 100;

	/** The Constant REALLOC_FACTOR. */
	public static final int REALLOC_FACTOR = 2;

	/** The Constant MAX_SIZE_MASK_0. */
	public static final double MAX_SIZE_MASK_0 = 3.09023230616781; /* Size for Gaussian mask */

	/** The Constant MAX_SIZE_MASK_1. */
	public static final double MAX_SIZE_MASK_1 = 3.46087178201605; /* Size for 1st derivative mask */

	/** The Constant MAX_SIZE_MASK_2. */
	public static final double MAX_SIZE_MASK_2 = 3.82922419517181; /* Size for 2nd derivative mask */

	/**
	 * Mask size.
	 *
	 * @param MAX
	 *            the max
	 * @param sigma
	 *            the sigma
	 * @return the int
	 */
	public static int MASK_SIZE(double MAX, double sigma) {
		return (int) Math.ceil(MAX * sigma); /* Maximum mask index */
	}

	/** The Constant ERR_NOMEM. */
	public static final String ERR_NOMEM = "Out of memory";

	/** The Constant ERR_FNF. */
	public static final String ERR_FNF = "Could not open file";

	/** The Constant ERR_NOPGM. */
	public static final String ERR_NOPGM = "Not a valid pgm file:";

	/** The Constant ERR_SNS. */
	public static final String ERR_SNS = "Sigma not specified";

	/** The Constant ERR_SOR. */
	public static final String ERR_SOR = "Sigma out of range:";

	/** The Constant ERR_LNS. */
	public static final String ERR_LNS = "Low not specified";

	/** The Constant ERR_LOR. */
	public static final String ERR_LOR = "Low out of range:";

	/** The Constant ERR_HNS. */
	public static final String ERR_HNS = "High not specified";

	/** The Constant ERR_HOR. */
	public static final String ERR_HOR = "High out of range:";

	/** The Constant ERR_LGH. */
	public static final String ERR_LGH = "Low > High";

	/** The Constant ERR_CNW. */
	public static final String ERR_CNW = "Line position correction impossible without line width";

	/** The Constant ERR_INP. */
	public static final String ERR_INP = "Include-image option requires PostScript output";

	/** The Constant ERR_UKO. */
	public static final String ERR_UKO = "Unknown option:";

	/** The Constant ERR_TMF. */
	public static final String ERR_TMF = "Too many files specified:";

	/**
	 * Lincoor.
	 *
	 * @param row
	 *            the row
	 * @param col
	 *            the col
	 * @param width
	 *            the width
	 * @return the int
	 */
	/*
	 * Translate row and column coordinates of an image into an index into its
	 * one-dimensional array.
	 */
	public static int LINCOOR(int row, int col, int width) {
		return row * width + col;
	}

	/**
	 * Br.
	 *
	 * @param row
	 *            the row
	 * @param height
	 *            the height
	 * @return the int
	 */
	/*
	 * Mirror the row coordinate at the borders of the image; height must be a
	 * defined variable in the calling function containing the image height.
	 */
	public static int BR(int row, int height) {
		return ((row) < 0 ? -(row) : (row) >= height ? height - (row) + height - 2 : (row));
	}

	/**
	 * Bc.
	 *
	 * @param col
	 *            the col
	 * @param width
	 *            the width
	 * @return the int
	 */
	/*
	 * Mirror the column coordinate at the borders of the image; width must be a
	 * defined variable in the calling function containing the image width.
	 */
	public static int BC(int col, int width) {
		return ((col) < 0 ? -(col) : (col) >= width ? width - (col) + width - 2 : (col));
	}

	/**
	 * The Enum contour_class.
	 */
	public enum contour_class {

		/** The cont no junc. */
		cont_no_junc,
		/** The cont start junc. */
		/* no end point is a junction */
		cont_start_junc,
		/** The cont end junc. */
		/* only the start point of the line is a junction */
		cont_end_junc,
		/** The cont both junc. */
		/* only the end point of the line is a junction */
		cont_both_junc,
		/** The cont closed. */
		/* both end points of the line are junctions */
		cont_closed; /* the contour is closed */
	}

}
