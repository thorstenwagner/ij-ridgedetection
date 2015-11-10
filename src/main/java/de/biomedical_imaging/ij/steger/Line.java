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

import ij.IJ;

/** This class holds one extracted line.  The field num contains the number of
   points in the line.  The coordinates of the line points are given in the
   arrays row and col.  The array angle contains the direction of the normal
   to each line point, as measured from the row-axis.  Some people like to
   call the col-axis the x-axis and the row-axis the y-axis, and measure the
   angle from the x-axis.  To convert the angle into this convention, subtract
   PI/2 from the angle and normalize it to be in the interval [0,2*PI).  The
   array response contains the response of the operator, i.e., the second
   directional derivative in the direction of angle, at each line point.  The
   arrays width_l and width_r contain the width information for each line point
   if the algorithm was requested to extract it; otherwise they are NULL.  If
   the line position and width correction was applied the contents of width_l
   and width_r will be identical.  The arrays asymmetry and contrast contain
   the true asymmetry and contrast of each line point if the algorithm was
   instructed to apply the width and position correction.  Otherwise, they are
   set to NULL.  If the asymmetry, i.e., the weaker gradient, is on the right
   side of the line, the asymmetry is set to a positive value, while if it is
   on the left side it is set to a negative value. */
public class Line {
	/** number of points */
	int  num;         
	
	/** row coordinates of the line points */
	float[] row;   
	
	/** column coordinates of the line points */
	float[] col; 
	
	/** angle of normal (measured from the row axis) */
	float[] angle;    
	
	/** response of line point (second derivative) */
	float[] response;    
	
	/** width to the left of the line */
	float[] width_l;    
	
	/** width to the right of the line */
	float[] width_r;    
	
	/** asymmetry of the line point */
	float[] asymmetry;      
	
	/** intensity of the line point */
	float[] intensity;          
	
	/** contour class (e.g., closed, no_junc) */
	private LinesUtil.contour_class cont_class; 
	
	static int idCounter = 0;
	private int id;
	private int frame;
	public Line() {
		// TODO Auto-generated constructor stub
		assignID();
		
	}
	
	public LinesUtil.contour_class getContourClass(){
		
		return cont_class;
	}
	
	public void setContourClass(LinesUtil.contour_class cont_class){
		this.cont_class = cont_class;
	}
	
	public void setFrame(int frame){
		this.frame = frame;
	}
	
	/**
	 * @return the frame index where the line was detected
	 */
	public int getFrame(){
		return frame;
	}
	
	/**
	 * @return x coordinates of the line points
	 */
	public float[] getXCoordinates(){
		return col;
	}
	
	/**
	 * @return y coordinates of the line points
	 */
	public float[] getYCoordinates(){
		return row;
	}
	/**
	 * @return response of line point (second derivative)
	 */
	public float[] getResponse(){
		return response;
	}
	/**
	 * 
	 * @return intensity of the line points
	 */
	public float[] getIntensity(){
		return intensity;
	}
	/**
	 * @return angle of normal (measured from the y-axis)
	 */
	public float[] getAngle(){
		return angle;
	}
	/**
	 * @return asymmetry of the line point
	 */
	public float[] getAsymmetry(){
		return asymmetry;
	}
	/**
	 * @return width to the left of the line
	 */
	public float[] getLineWidthL(){
		return width_l;
	}
	
	/**
	 * @return width to the right of the line
	 */
	public float[] getLineWidthR(){
		return width_r;
	}
	
	/**
	 * 
	 * @return Return the number of points
	 */
	public int getNumber(){
		return  num;
	}
	
	/**
	 * @return unique ID of the line
	 */
	public int getID(){
		return id;
	}
	
	
	public LinesUtil.contour_class getLineClass(){
		return cont_class;
	}

	public int getStartOrdEndPosition(float x, float y){
		double distStart = Math.sqrt( Math.pow(col[0]-x, 2)+Math.pow(row[0]-y, 2));
		double distEnd = Math.sqrt( Math.pow(col[ (this.num-1)]-x, 2)+Math.pow(row[ (this.num-1)]-y, 2));
		return  (distStart < distEnd ? 0:this.num-1);
	}
	
	/**
	 * 
	 * @return the estimated length of the line
	 */
	public double estimateLength(){
		double length = 0;
		for( int i = 1 ; i < num; i++){
			length += Math.sqrt(Math.pow(col[i]-col[i-1], 2)+Math.pow(row[i]-row[i-1], 2));
		}
		return length;
	}
	
	private synchronized  void assignID(){
		this.id = idCounter;
		idCounter++;
	}
	
	static void resetCounter(){
		idCounter = 0;
	}
}
