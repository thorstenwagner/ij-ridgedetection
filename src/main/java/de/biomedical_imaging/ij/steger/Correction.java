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
 * The Class Correction.
 */
public class Correction {

	/** The w est. */
	double w_est; /* Total line width extracted from the image */

	/** The r est. */
	double r_est; /* Gradient ratio extracted from the image */

	/** The w. */
	double w; /* True line width */

	/** The h. */
	double h; /* True asymmetry */

	/** The correction. */
	double correction; /* Line position correction */

	/** The w strong. */
	double w_strong; /* True width on the side with the stronger gradient */

	/** The w weak. */
	double w_weak; /* True width on the side with the weaker gradient */

	/** The is valid. */
	boolean is_valid; /* Is this table entry valid? */

	/**
	 * Instantiates a new correction.
	 *
	 * @param w_est
	 *            the w est
	 * @param r_rest
	 *            the r rest
	 * @param w
	 *            the w
	 * @param h
	 *            the h
	 * @param correction
	 *            the correction
	 * @param w_strong
	 *            the w strong
	 * @param w_weak
	 *            the w weak
	 * @param is_valid
	 *            the is valid
	 */
	public Correction(double w_est, double r_rest, double w, double h, double correction, double w_strong,
			double w_weak, boolean is_valid) {
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
