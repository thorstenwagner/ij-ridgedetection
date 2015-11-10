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
/* This data structure facilitates the quick search for the next possible
   starting point of a line.  An array of crossrefs will be accumulated and
   sorted according to its value.  x and y are the coordinates of a point in
   the image.  When this point has been processed it will be marked as done. */
public class Crossref implements Comparable<Crossref>{
	 short  x;
	 short  y;
	 double value;
	 boolean   done;
	 /*
 int compare_crossrefs(crossref p1, crossref p2) {
		if (p1.value > p2.value)
			return -1;
		if (p1.value < p2.value)
			return 1;
		return 0;
	}
	*/
	
	/*
	 * This function compares two crossrefs according to their value. It is
	 * called by qsort.
	 */
	 @Override
	public int compareTo(Crossref arg0) {
		if (this.value > arg0.value)
			return -1;
		if (this.value < arg0.value)
			return 1;
		return 0;
	}
}
