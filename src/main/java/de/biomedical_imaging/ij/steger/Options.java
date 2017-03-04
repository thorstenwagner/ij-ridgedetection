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
/* Command line options for the program */
public class Options {
	  double sigma;
	  double low;
	  double high;
	  int   mode;
      double minLength;
      double maxLength;
	  boolean   correct;
	  boolean   width;
	  boolean   extend;
	  boolean   postscript;
	  boolean   encapsulated;
	  boolean   image;
		OverlapOption overlap;
	  
	  public Options(double sigma, double low, double high, int mode, double minLength, double maxLength, boolean correct, boolean width, boolean extend, boolean postscript, boolean encapsulated, boolean image, OverlapOption overlap) {
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

	public double getSigma() {
		return sigma;
	}

	public void setSigma(double sigma) {
		this.sigma = sigma;
	}

	public double getLow() {
		return low;
	}

	public void setLow(double low) {
		this.low = low;
	}

	public double getHigh() {
		return high;
	}

	public void setHigh(double high) {
		this.high = high;
	}

	public int getMode() {
		return mode;
	}
    
    public double getminLength() {
        return minLength;
    }
    
    public void setminLength(double minLength) {
        this.minLength = minLength;
    }
    
    public double getmaxLength() {
        return maxLength;
    }
    
    public void setmaxLength(double maxLength) {
        this.maxLength = maxLength;
    }

	public void setMode(int mode) {
		this.mode = mode;
	}

	public boolean isCorrect() {
		return correct;
	}

	public void setCorrect(boolean correct) {
		this.correct = correct;
	}

	public boolean isWidth() {
		return width;
	}

	public void setWidth(boolean width) {
		this.width = width;
	}

	public boolean isExtend() {
		return extend;
	}

	public void setExtend(boolean extend) {
		this.extend = extend;
	}

	public boolean isPostscript() {
		return postscript;
	}

	public void setPostscript(boolean postscript) {
		this.postscript = postscript;
	}

	public boolean isEncapsulated() {
		return encapsulated;
	}

	public void setEncapsulated(boolean encapsulated) {
		this.encapsulated = encapsulated;
	}

	public boolean isImage() {
		return image;
	}

	public void setImage(boolean image) {
		this.image = image;
	}

	public OverlapOption getOverlapResolution() {
		return overlap;
	}

	public void setOverlapResolution(OverlapOption overlap) {
		this.overlap = overlap;
	}
}
