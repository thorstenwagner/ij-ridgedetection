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
/* Command line options for the program */
public class Options {
	  double sigma;
	  double low;
	  double high;
	  long   mode;
	  boolean   correct;
	  boolean   width;
	  boolean   extend;
	  boolean   postscript;
	  boolean   encapsulated;
	  boolean   image;
	  
	  public Options(double sigma, double low, double high, long mode, boolean correct, boolean width, boolean extend, boolean postscript, boolean encapsulated, boolean image) {
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

	public long getMode() {
		return mode;
	}

	public void setMode(long mode) {
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
}
