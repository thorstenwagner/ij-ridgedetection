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
 * This data structure is used to accumulate junction information. It is needed
 * to split lines at junction points.
 */
public class Junction implements Comparable<Junction> {

	/** Index of line that is already processed. */
	int cont1;

	/** Index of line tnat runs into cont1. */
	int cont2;

	/** Index of the junction point in cont1. */
	int pos;

	/** x-(row-)coordinate of the junction point. */
	float x;

	/** y-(col-)coordinate of the junction point. */
	float y;

	/** line that is already processed. */
	Line lineCont1;

	/** line that runs into idCont1. */
	Line lineCont2;

	/** True if this junction sits on a start/end of at least one line. */
	boolean isNonTerminal = false;

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(Junction o) {
		 // Total order. Comparing only (cont1, pos) left junctions sharing both unordered, so the
		// line-splitting loop (grouping by cont1, splitting at pos) could split differently between
		// runs. Tie-break on cont2, then x, then y. Also use Integer/Float.compare instead of
		// subtraction to avoid integer overflow.
		if (this.cont1 != o.cont1) return Integer.compare(this.cont1, o.cont1);
		if (this.pos   != o.pos)   return Integer.compare(this.pos,   o.pos);
		if (this.cont2 != o.cont2) return Integer.compare(this.cont2, o.cont2);
		if (this.x     != o.x)     return Float.compare(this.x, o.x);
		return Float.compare(this.y, o.y);
	}

	/**
	 * Gets the x.
	 *
	 * @return x-coordinate of the junction point
	 */
	public float getX() {
		return x;
	}

	/**
	 * Gets the y.
	 *
	 * @return y-coordinate of the junction point
	 */
	public float getY() {
		return y;
	}

	/**
	 * Gets the line 1.
	 *
	 * @return The line that is already processed
	 */
	public Line getLine1() {
		return lineCont1;
	}

	/**
	 * Gets the line 2.
	 *
	 * @return The line that runs into line1
	 */
	public Line getLine2() {
		return lineCont2;
	}

	/**
	 * Checks if is non terminal.
	 *
	 * @return True iff this junction point does not sit on either line's start/end
	 */
	public boolean isNonTerminal() {
		return isNonTerminal;
	}

}
