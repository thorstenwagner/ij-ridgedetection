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
 * The Class Options.
 */
/* Command line options for the program */
public class Options {

	/** The sigma. */
	double sigma;

	/** The low. */
	double low;

	/** The high. */
	double high;

	/** The mode. */
	int mode;

	/** The min length. */
	double minLength;

	/** The max length. */
	double maxLength;

	/** The correct. */
	boolean correct;

	/** The width. */
	boolean width;

	/** The extend. */
	boolean extend;

	/** The postscript. */
	boolean postscript;

	/** The encapsulated. */
	boolean encapsulated;

	/** The image. */
	boolean image;

	/** The overlap. */
	OverlapOption overlap;

	/**
	 * Instantiates a new options.
	 *
	 * @param sigma
	 *            the sigma
	 * @param low
	 *            the low
	 * @param high
	 *            the high
	 * @param mode
	 *            the mode
	 * @param minLength
	 *            the min length
	 * @param maxLength
	 *            the max length
	 * @param correct
	 *            the correct
	 * @param width
	 *            the width
	 * @param extend
	 *            the extend
	 * @param postscript
	 *            the postscript
	 * @param encapsulated
	 *            the encapsulated
	 * @param image
	 *            the image
	 * @param overlap
	 *            the overlap
	 */
	public Options(double sigma, double low, double high, int mode, double minLength, double maxLength, boolean correct,
			boolean width, boolean extend, boolean postscript, boolean encapsulated, boolean image,
			OverlapOption overlap) {
		// TODO Auto-generated constructor stub
		this.sigma = sigma;
		this.low = low;
		this.high = high;
		this.mode = mode;
		this.correct = correct;
		this.width = width;
		this.extend = extend;
		this.postscript = postscript;
		this.encapsulated = encapsulated;
		this.image = image;
		this.overlap = overlap;
	}

	/**
	 * Gets the sigma.
	 *
	 * @return the sigma
	 */
	public double getSigma() {
		return sigma;
	}

	/**
	 * Sets the sigma.
	 *
	 * @param sigma
	 *            the new sigma
	 */
	public void setSigma(double sigma) {
		this.sigma = sigma;
	}

	/**
	 * Gets the low.
	 *
	 * @return the low
	 */
	public double getLow() {
		return low;
	}

	/**
	 * Sets the low.
	 *
	 * @param low
	 *            the new low
	 */
	public void setLow(double low) {
		this.low = low;
	}

	/**
	 * Gets the high.
	 *
	 * @return the high
	 */
	public double getHigh() {
		return high;
	}

	/**
	 * Sets the high.
	 *
	 * @param high
	 *            the new high
	 */
	public void setHigh(double high) {
		this.high = high;
	}

	/**
	 * Gets the mode.
	 *
	 * @return the mode
	 */
	public int getMode() {
		return mode;
	}

	/**
	 * Gets the min length.
	 *
	 * @return the min length
	 */
	public double getminLength() {
		return minLength;
	}

	/**
	 * Sets the min length.
	 *
	 * @param minLength
	 *            the new min length
	 */
	public void setminLength(double minLength) {
		this.minLength = minLength;
	}

	/**
	 * Gets the max length.
	 *
	 * @return the max length
	 */
	public double getmaxLength() {
		return maxLength;
	}

	/**
	 * Sets the max length.
	 *
	 * @param maxLength
	 *            the new max length
	 */
	public void setmaxLength(double maxLength) {
		this.maxLength = maxLength;
	}

	/**
	 * Sets the mode.
	 *
	 * @param mode
	 *            the new mode
	 */
	public void setMode(int mode) {
		this.mode = mode;
	}

	/**
	 * Checks if is correct.
	 *
	 * @return true, if is correct
	 */
	public boolean isCorrect() {
		return correct;
	}

	/**
	 * Sets the correct.
	 *
	 * @param correct
	 *            the new correct
	 */
	public void setCorrect(boolean correct) {
		this.correct = correct;
	}

	/**
	 * Checks if is width.
	 *
	 * @return true, if is width
	 */
	public boolean isWidth() {
		return width;
	}

	/**
	 * Sets the width.
	 *
	 * @param width
	 *            the new width
	 */
	public void setWidth(boolean width) {
		this.width = width;
	}

	/**
	 * Checks if is extend.
	 *
	 * @return true, if is extend
	 */
	public boolean isExtend() {
		return extend;
	}

	/**
	 * Sets the extend.
	 *
	 * @param extend
	 *            the new extend
	 */
	public void setExtend(boolean extend) {
		this.extend = extend;
	}

	/**
	 * Checks if is postscript.
	 *
	 * @return true, if is postscript
	 */
	public boolean isPostscript() {
		return postscript;
	}

	/**
	 * Sets the postscript.
	 *
	 * @param postscript
	 *            the new postscript
	 */
	public void setPostscript(boolean postscript) {
		this.postscript = postscript;
	}

	/**
	 * Checks if is encapsulated.
	 *
	 * @return true, if is encapsulated
	 */
	public boolean isEncapsulated() {
		return encapsulated;
	}

	/**
	 * Sets the encapsulated.
	 *
	 * @param encapsulated
	 *            the new encapsulated
	 */
	public void setEncapsulated(boolean encapsulated) {
		this.encapsulated = encapsulated;
	}

	/**
	 * Checks if is image.
	 *
	 * @return true, if is image
	 */
	public boolean isImage() {
		return image;
	}

	/**
	 * Sets the image.
	 *
	 * @param image
	 *            the new image
	 */
	public void setImage(boolean image) {
		this.image = image;
	}

	/**
	 * Gets the overlap resolution.
	 *
	 * @return the overlap resolution
	 */
	public OverlapOption getOverlapResolution() {
		return overlap;
	}

	/**
	 * Sets the overlap resolution.
	 *
	 * @param overlap
	 *            the new overlap resolution
	 */
	public void setOverlapResolution(OverlapOption overlap) {
		this.overlap = overlap;
	}
}
