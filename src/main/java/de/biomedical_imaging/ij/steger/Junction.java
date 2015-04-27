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
/** This data structure is used to accumulate junction information.  It is
   needed to split lines at junction points. */
public class Junction implements Comparable<Junction> {
	 /** Index of line that is already processed */
	 long  cont1; 
	 /** Index of line tnat runs into cont1 */
	 long  cont2; 
	 /** Index of the junction point in cont1 */
	 long  pos; 
	 /** x-(row-)coordinate of the junction point */
	 float x;     
	 /** y-(col-)coordinate of the junction point */
	 float y;     
	 /** line that is already processed */
	 Line lineCont1; 
	 /** line that runs into idCont1 */
	 Line lineCont2; 
	 

	@Override
	public int compareTo(Junction o) {
		return (int) (((this.cont1 - o.cont1) != 0) ? this.cont1 - o.cont1
				: this.pos - o.pos);
	}
	
	/**
	 * @return  x-coordinate of the junction point
	 */
	public float getX(){
		return x;
	}
	
	/**
	 * @return  y-coordinate of the junction point
	 */
	public float getY(){
		return y;
	}
	
	/**
	 * @return  The line that is already processed
	 */
	public Line getLine1(){
		return lineCont1;
	}
	
	/**
	 * @return  The line that runs into line1
	 */
	public Line getLine2(){
		return lineCont2;
	}
	
	
	 
}
