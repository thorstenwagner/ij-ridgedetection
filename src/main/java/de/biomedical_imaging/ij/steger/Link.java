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

import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;

public class Link {

	public static final double MAX_ANGLE_DIFFERENCE = Math.PI / 6.0;
	// public static final double MAX_LINE_EXTENSION = 2.5*sigma;

	/*
	 * This table contains the three appropriate neighbor pixels that the
	 * linking algorithm must examine. It is indexed by the octant the current
	 * line angle lies in, e.g., 0 if the angle in degrees lies within
	 * [-22.5,22.5].
	 */
	final static int[][][] dirtab = new int[][][] {
			{ { 1, 0 }, { 1, -1 }, { 1, 1 } },
			{ { 1, 1 }, { 1, 0 }, { 0, 1 } },
			{ { 0, 1 }, { 1, 1 }, { -1, 1 } },
			{ { -1, 1 }, { 0, 1 }, { -1, 0 } },
			{ { -1, 0 }, { -1, 1 }, { -1, -1 } },
			{ { -1, -1 }, { -1, 0 }, { 0, -1 } },
			{ { 0, -1 }, { -1, -1 }, { 1, -1 } },
			{ { 1, -1 }, { 0, -1 }, { 1, 0 } } };

	/*
	 * This table contains the two neighbor pixels that the linking algorithm
	 * should examine and mark as processed in case there are double responses.
	 */
	final static int[][][] cleartab = new int[][][] { { { 0, 1 }, { 0, -1 } },
			{ { -1, 1 }, { 1, -1 } }, { { -1, 0 }, { 1, 0 } },
			{ { -1, -1 }, { 1, 1 } }, { { 0, -1 }, { 0, 1 } },
			{ { 1, -1 }, { -1, 1 } }, { { 1, 0 }, { -1, 0 } },
			{ { 1, 1 }, { -1, -1 } } };

	/*
	 * Compute the response of the operator with sub-pixel accuracy by using the
	 * facet model to interpolate the pixel accurate responses.
	 */
	private double interpolate_response(float[] resp, int x, int y,
			double px, double py, int width, int height) {
		double i1, i2, i3, i4, i5, i6, i7, i8, i9;
		double t1, t2, t3, t4, t5, t6;
		double d, dr, dc, drr, drc, dcc;
		double xx, yy;

		i1 = resp[ LinesUtil.LINCOOR(LinesUtil.BR(x - 1, height),
				LinesUtil.BC(y - 1, width), width)];
		i2 = resp[ LinesUtil
				.LINCOOR(LinesUtil.BR(x - 1, height), y, width)];
		i3 = resp[ LinesUtil.LINCOOR(LinesUtil.BR(x - 1, height),
				LinesUtil.BC(y + 1, width), width)];
		i4 = resp[ LinesUtil.LINCOOR(x, LinesUtil.BC(y - 1, width), width)];
		i5 = resp[ LinesUtil.LINCOOR(x, y, width)];
		i6 = resp[ LinesUtil.LINCOOR(x, LinesUtil.BC(y + 1, width), width)];
		i7 = resp[ LinesUtil.LINCOOR(LinesUtil.BR(x + 1, height),
				LinesUtil.BC(y - 1, width), width)];
		i8 = resp[ LinesUtil
				.LINCOOR(LinesUtil.BR(x + 1, height), y, width)];
		i9 = resp[ LinesUtil.LINCOOR(LinesUtil.BR(x + 1, height),
				LinesUtil.BC(y + 1, width), width)];
		t1 = i1 + i2 + i3;
		t2 = i4 + i5 + i6;
		t3 = i7 + i8 + i9;
		t4 = i1 + i4 + i7;
		t5 = i2 + i5 + i8;
		t6 = i3 + i6 + i9;
		d = (-i1 + 2 * i2 - i3 + 2 * i4 + 5 * i5 + 2 * i6 - i7 + 2 * i8 - i9) / 9;
		dr = (t3 - t1) / 6;
		dc = (t6 - t4) / 6;
		drr = (t1 - 2 * t2 + t3) / 6;
		dcc = (t4 - 2 * t5 + t6) / 6;
		drc = (i1 - i3 - i7 + i9) / 4;
		xx = px - x;
		yy = py - y;
		return d + xx * dr + yy * dc + xx * xx * drr + xx * yy * drc + yy * yy
				* dcc;
	}

	/*
	 * Calculate the closest point to (px,py) on the line (lx,ly) + t*(dx,dy)
	 * and return the result in (cx,cy), plus the parameter in t.
	 */
	private void closest_point(double lx, double ly, double dx, double dy,
			double px, double py, MutableDouble cx, MutableDouble cy,
			MutableDouble t) {
		double mx, my, den, nom, tt;
		mx = px - lx;
		my = py - ly;
		den = dx * dx + dy * dy;
		nom = mx * dx + my * dy;
		if (den != 0)
			tt = nom / den;
		else
			tt = 0;
		cx.setValue(lx + tt * dx);
		cy.setValue(ly + tt * dy);
		t.setValue(tt);
	}

	/*
	 * Interpolate the gradient of the gradient images gradx and grady with
	 * width width at the point (px,py) using linear interpolation, and return
	 * the result in (gx,gy).
	 */
	private void interpolate_gradient(float[] gradx, float[] grady, double px,
			double py, int width, MutableDouble gx, MutableDouble gy) {
		int gix, giy, gpos;
		double gfx, gfy, gx1, gy1, gx2, gy2, gx3, gy3, gx4, gy4;

		gix =  (int)Math.floor(px);
		giy =  (int)Math.floor(py);

		gfx = px % 1.0;
		;
		gfy = py % 1.0;
		gpos = LinesUtil.LINCOOR(gix, giy, width);
		gx1 = gradx[ gpos];
		gy1 = grady[ gpos];
		gpos = LinesUtil.LINCOOR(gix + 1, giy, width);
		gx2 = gradx[ gpos];
		gy2 = grady[ gpos];
		gpos = LinesUtil.LINCOOR(gix, giy + 1, width);
		gx3 = gradx[ gpos];
		gy3 = grady[ gpos];
		gpos = LinesUtil.LINCOOR(gix + 1, giy + 1, width);
		gx4 = gradx[ gpos];
		gy4 = grady[ gpos];
		gx.setValue((1 - gfy) * ((1 - gfx) * gx1 + gfx * gx2) + gfy
				* ((1 - gfx) * gx3 + gfx * gx4));
		gy.setValue((1 - gfy) * ((1 - gfx) * gy1 + gfx * gy2) + gfy
				* ((1 - gfx) * gy3 + gfx * gy4));
	}

	/*
	 * This function links the line points into lines. The input to this
	 * function are the response of the filter, i.e., the second directional
	 * derivative along (nx[l],ny[l]), contained in eigval[l], and the sub-pixel
	 * position of each line point, contained in (px[l],py[l]). The parameters
	 * low and high are the hysteresis thresholds for the linking, while width
	 * and height are the dimensions of the five float-images. The linked lines
	 * are returned in result, and the number of lines detected is returned in
	 * num_result.
	 */
	public void compute_contours(byte[] ismax, float[] eigval, float[] normx,
			float[] normy, float[] posx, float[] posy, float[] gradx,
			float[] grady, Lines contours, MutableInt num_result,
			double sigma, boolean extend_lines, int mode, double low,
			double high, int width, int height, Junctions junctions) {
		int i = 0, j = 0, k, l, it, pos, nextpos, nexti;
		int begin, end;
		int x, y;
		int octant, last_octant;
		int[] label;
		int num_cont, num_pnt;
		int size_cont, size_pnt;
		float[] row, col, trow, tcol;
		float[] angle, tangle;
		float[] resp, tresp;
		Junction[] junc;
		int num_junc, size_junc;
		Line[] cont;
		Line tmp_cont;
		LinesUtil.contour_class cls;
		double max;
		int maxx, maxy;
		int nextx, nexty;
		double nx, ny;
		double alpha, nextalpha, diff, mindiff, dist, mindist;
		double beta, last_beta, diff1, diff2;
		double px, py, nextpx = 0, nextpy = 0;
		double dx, dy;
		float tmp;
		int area;
		int[] indx;
		int indx_max;
		boolean nextismax;
		Crossref[] cross;
		int m = 0, max_line, num_add;
		MutableInt num_line = new MutableInt();
		double length, response;
		Offset[] line;
		double mx, my, gx, gy, s, end_angle = 0, end_resp = 0;
		MutableDouble t = new MutableDouble();
		float[] extx, exty;
		boolean add_ext;
		Region seg = new Region();
		Chord[] rl;
		Width w = new Width();

		/*
		 * The image label contains information on the pixels that have been
		 * processed by the linking algorithm.
		 */
		label = new int[ (width * height)];
		java.util.Arrays.fill(label, 0); // Vermutlich nicht notwendig, da
											// standardmäßig 0.

		/*
		 * The image indx is an index into the table of all pixels that possibly
		 * could be starting points for new lines. It is used to quickly
		 * determine the next starting point of a line.
		 */
		indx = new int[ (width * height)];
		java.util.Arrays.fill(indx, 0);

		num_cont = 0;
		num_junc = 0;
		size_cont = LinesUtil.INITIAL_SIZE;
		size_pnt = LinesUtil.INITIAL_SIZE;
		size_junc = LinesUtil.INITIAL_SIZE;
		cont = new Line[ size_cont];
		for (int o = 0; o < cont.length; o++) {
			cont[o] = new Line();
		}
		row = new float[ size_pnt];
		col = new float[ size_pnt];
		angle = new float[ size_pnt];
		resp = new float[ size_pnt];
		junc = new Junction[ size_junc];
		for (int o = 0; o < junc.length; o++) {
			junc[o] = new Junction();
		}

		/* Select all pixels that can be starting points for lines. */
		Threshold.threshold(ismax, 2, width, height, seg);

		/* Count the number of possible starting points. */
		area = 0;
		for (i = 0; i < seg.num; i++)
			area += seg.rl[ i].ce - seg.rl[ i].cb + 1;

		/* Create the index of possible starting points. */
		cross = new Crossref[ area];
		for (int o = 0; o < cross.length; o++) {
			cross[o] = new Crossref();
		}
		k = 0;
		rl = seg.rl;
		for (i = 0; i < seg.num; i++) {
			x = rl[ i].r;
			for (y = rl[ i].cb; y <= rl[ i].ce; y++) {
				pos = LinesUtil.LINCOOR(x, y, width);
				cross[ k].x = (short) x;
				cross[ k].y = (short) y;
				cross[ k].value = eigval[ pos];
				cross[ k].done = false;
				k++;
			}
		}
		
		java.util.Arrays.sort(cross);
		// qsort(cross,area,sizeof(*cross),compare_crossrefs);
		for (i = 0; i < area; i++)
			indx[ LinesUtil.LINCOOR(cross[ i].x, cross[ i].y,
					width)] = i + 1;

		/* Link lines points. */
		indx_max = 0;
		for (;;) {
			/*
			 * Contour class unknown at this point; therefore assume both ends
			 * free.
			 */
			cls = LinesUtil.contour_class.cont_no_junc;
			/* Search for next starting point. */
			while (indx_max < area && cross[ indx_max].done)
				indx_max++;
			/* Stop if no feasible starting point exists. */
			if (indx_max == area)
				break;
			max = cross[ indx_max].value;
			maxx = cross[ indx_max].x;
			maxy = cross[ indx_max].y;
			if (max == 0.0)
				break;

			/* Add starting point to the line. */
			num_pnt = 0;
			pos = LinesUtil.LINCOOR(maxx, maxy, width);
			label[ pos] =  (num_cont + 1);
			if (!(indx[ pos] == 0))
				cross[ (indx[ pos] - 1)].done = true;
			row[ num_pnt] = posx[ pos];
			col[ num_pnt] = posy[ pos];
			/* Select line direction. */
			nx = -normy[ pos];
			ny = normx[ pos];
			alpha = Math.atan2(ny, nx);
			if (alpha < 0.0)
				alpha += 2.0 * Math.PI;
			if (alpha >= Math.PI)
				alpha -= Math.PI;
			octant =  (int)(Math.floor(4.0 / Math.PI * alpha + 0.5)) % 4;
			/*
			 * Select normal to the line. The normal points to the right of the
			 * line as the line is traversed from 0 to num-1. Since the points
			 * are sorted in reverse order before the second iteration, the
			 * first beta actually has to point to the left of the line!
			 */
			beta = alpha + Math.PI / 2.0;
			if (beta >= 2.0 * Math.PI)
				beta -= 2.0 * Math.PI;
			angle[ num_pnt] = (float) beta;
			resp[ num_pnt] = (float) interpolate_response(eigval, maxx,
					maxy, posx[ pos], posy[ pos], width, height);
			num_pnt++;
			/* Mark double responses as processed. */
			for (i = 0; i < 2; i++) {
				nextx = maxx + cleartab[ octant][ i][0];
				nexty = maxy + cleartab[ octant][ i][1];
				if (nextx < 0 || nextx >= height || nexty < 0 || nexty >= width)
					continue;
				nextpos = LinesUtil.LINCOOR(nextx, nexty, width);
				if (ismax[ nextpos] > 0) {
					nx = -normy[ nextpos];
					ny = normx[ nextpos];
					nextalpha = Math.atan2(ny, nx);
					if (nextalpha < 0.0)
						nextalpha += 2.0 * Math.PI;
					if (nextalpha >= Math.PI)
						nextalpha -= Math.PI;
					diff = Math.abs(alpha - nextalpha);
					if (diff >= Math.PI / 2.0)
						diff = Math.PI - diff;
					if (diff < MAX_ANGLE_DIFFERENCE) {
						label[ nextpos] =  (num_cont + 1);
						if (!(indx[ nextpos] == 0))
							cross[ (indx[ nextpos] - 1)].done = true;
					}
				}
			}

			for (it = 1; it <= 2; it++) {
				if (it == 1) {
					/*
					 * Search along the initial line direction in the first
					 * iteration.
					 */
					x = maxx;
					y = maxy;
					pos = LinesUtil.LINCOOR(x, y, width);
					nx = -normy[ pos];
					ny = normx[ pos];
					alpha = Math.atan2(ny, nx);
					if (alpha < 0.0)
						alpha += 2.0 * Math.PI;
					if (alpha >= Math.PI)
						alpha -= Math.PI;
					last_octant =  (int)(Math.floor(4.0 / Math.PI * alpha
							+ 0.5)) % 4;
					last_beta = alpha + Math.PI / 2.0;
					if (last_beta >= 2.0 * Math.PI)
						last_beta -= 2.0 * Math.PI;
				} else {
					/* Search in the opposite direction in the second iteration. */
					x = maxx;
					y = maxy;
					pos = LinesUtil.LINCOOR(x, y, width);
					nx = -normy[ pos];
					ny = normx[ pos];
					alpha = Math.atan2(ny, nx);
					if (alpha < 0.0)
						alpha += 2.0 * Math.PI;
					if (alpha >= Math.PI)
						alpha -= Math.PI;
					last_octant =  (int)(Math.floor(4.0 / Math.PI * alpha
							+ 0.5)) % 4 + 4;
					last_beta = alpha + Math.PI / 2.0;
					if (last_beta >= 2.0 * Math.PI)
						last_beta -= 2.0 * Math.PI;
				}
				if (it == 2) {
					/* Sort the points found in the first iteration in reverse. */
					for (i = 0; i < num_pnt / 2; i++) {
						tmp = row[ i];
						row[ i] = row[ (num_pnt - 1 - i)];
						row[ (num_pnt - 1 - i)] = tmp;
						tmp = col[ i];
						col[ i] = col[ (num_pnt - 1 - i)];
						col[ (num_pnt - 1 - i)] = tmp;
						tmp = angle[ i];
						angle[ i] = angle[ (num_pnt - 1 - i)];
						angle[ (num_pnt - 1 - i)] = tmp;
						tmp = resp[ i];
						resp[ i] = resp[ (num_pnt - 1 - i)];
						resp[ (num_pnt - 1 - i)] = tmp;
					}
				}

				/* Now start adding appropriate neighbors to the line. */
				for (;;) {
					pos = LinesUtil.LINCOOR(x, y, width);
					nx = -normy[ pos];
					ny = normx[ pos];
					px = posx[ pos];
					py = posy[ pos];
					/* Orient line direction w.r.t. the last line direction. */
					alpha = Math.atan2(ny, nx);
					if (alpha < 0.0)
						alpha += 2.0 * Math.PI;
					if (alpha >= Math.PI)
						alpha -= Math.PI;
					octant =  (int)(Math.floor(4.0 / Math.PI * alpha + 0.5)) % 4;
					switch ( octant) {
					case 0:
						if (last_octant >= 3 && last_octant <= 5)
							octant = 4;
						break;
					case 1:
						if (last_octant >= 4 && last_octant <= 6)
							octant = 5;
						break;
					case 2:
						if (last_octant >= 4 && last_octant <= 7)
							octant = 6;
						break;
					case 3:
						if (last_octant == 0 || last_octant >= 6)
							octant = 7;
						break;
					}
					last_octant = octant;

					/* Determine appropriate neighbor. */
					nextismax = false;
					nexti = 1;
					mindiff = Double.MAX_VALUE;
					for (i = 0; i < 3; i++) {
						nextx = x + dirtab[ octant][ i][0];
						nexty = y + dirtab[ octant][ i][1];
						if (nextx < 0 || nextx >= height || nexty < 0
								|| nexty >= width)
							continue;
						nextpos = LinesUtil.LINCOOR(nextx, nexty, width);
						if (ismax[ nextpos] == 0)
							continue;
						nextpx = posx[ nextpos];
						nextpy = posy[ nextpos];
						dx = nextpx - px;
						dy = nextpy - py;
						dist = Math.sqrt(dx * dx + dy * dy);
						nx = -normy[ nextpos];
						ny = normx[ nextpos];
						nextalpha = Math.atan2(ny, nx);
						if (nextalpha < 0.0)
							nextalpha += 2.0 * Math.PI;
						if (nextalpha >= Math.PI)
							nextalpha -= Math.PI;
						diff = Math.abs(alpha - nextalpha);
						if (diff >= Math.PI / 2.0)
							diff = Math.PI - diff;
						diff = dist + diff;
						if (diff < mindiff) {
							mindiff = diff;
							nexti = i;
						}
						if (!(ismax[ nextpos] == 0))
							nextismax = true;
					}

					/* Mark double responses as processed. */
					for (i = 0; i < 2; i++) {
						nextx = x + cleartab[ octant][ i][0];
						nexty = y + cleartab[ octant][ i][1];
						if (nextx < 0 || nextx >= height || nexty < 0
								|| nexty >= width)
							continue;
						nextpos = LinesUtil.LINCOOR(nextx, nexty, width);
						if (ismax[ nextpos] > 0) {
							nx = -normy[ nextpos];
							ny = normx[ nextpos];
							nextalpha = Math.atan2(ny, nx);
							if (nextalpha < 0.0)
								nextalpha += 2.0 * Math.PI;
							if (nextalpha >= Math.PI)
								nextalpha -= Math.PI;
							diff = Math.abs(alpha - nextalpha);
							if (diff >= Math.PI / 2.0)
								diff = Math.PI - diff;
							if (diff < MAX_ANGLE_DIFFERENCE) {
								label[ nextpos] =  (num_cont + 1);
								if (!(indx[ nextpos] == 0))
									cross[ (indx[ nextpos] - 1)].done = true;
							}
						}
					}

					/* Have we found the end of the line? */
					if (!nextismax)
						break;
					/* If not, add the neighbor to the line. */
					x += dirtab[ octant][ nexti][0];
					y += dirtab[ octant][ nexti][1];
					if (num_pnt >= size_pnt) {
						size_pnt =  (int)Math
								.floor((double) (size_pnt * LinesUtil.REALLOC_FACTOR));
						float[] newArr = new float[ size_pnt];
						for (int o = 0; o < row.length; o++) {
							newArr[o] = row[o];
						}
						row = newArr;
						
						newArr = new float[ size_pnt];
						for (int o = 0; o < col.length; o++) {
							newArr[o] = col[o];
						}
						col = newArr;
						
						newArr = new float[ size_pnt];
						for (int o = 0; o < angle.length; o++) {
							newArr[o] = angle[o];
						}
						angle = newArr;
						
						newArr = new float[ size_pnt];
						for (int o = 0; o < resp.length; o++) {
							resp[o] = resp[o];
						}
						resp = newArr;
					}
					pos = LinesUtil.LINCOOR(x, y, width);
					row[ num_pnt] = posx[ pos];
					col[ num_pnt] = posy[ pos];

					/*
					 * Orient normal to the line direction w.r.t. the last
					 * normal.
					 */
					nx = normx[ pos];
					ny = normy[ pos];
					beta = Math.atan2(ny, nx);
					if (beta < 0.0)
						beta += 2.0 * Math.PI;
					if (beta >= Math.PI)
						beta -= Math.PI;
					diff1 = Math.abs(beta - last_beta);
					if (diff1 >= Math.PI)
						diff1 = 2.0 * Math.PI - diff1;
					diff2 = Math.abs(beta + Math.PI - last_beta);
					if (diff2 >= Math.PI)
						diff2 = 2.0 * Math.PI - diff2;
					if (diff1 < diff2) {
						angle[ num_pnt] = (float) beta;
						last_beta = beta;
					} else {
						angle[ num_pnt] = (float) (beta + Math.PI);
						last_beta = beta + Math.PI;
					}

					resp[ num_pnt] = (float) interpolate_response(eigval,
							x, y, posx[ pos], posy[ pos], width,
							height);
					num_pnt++;

					/*
					 * If the appropriate neighbor is already processed a
					 * junction point is found.
					 */
					if (label[ pos] > 0) {
						if (num_junc >= size_junc) {
							size_junc =  (int)Math
									.floor((double) (size_junc * LinesUtil.REALLOC_FACTOR));
							Junction[] junch = new Junction[ size_junc];
							for (int o = 0; o < junch.length; o++) {
								if (o < junc.length)
									junch[o] = junc[o];
								else
									junch[o] = new Junction();
							}
							junc = junch;
						}
						/* Look for the junction point in the other line. */
						k = label[ pos] - 1;
						if (k == num_cont) {
							/* Line intersects itself. */
							for (j = 0; j < num_pnt - 1; j++) {
								if (row[ j] == posx[ pos]
										&& col[ j] == posy[ pos]) {
									if (j == 0) {
										/* Contour is closed. */
										cls = LinesUtil.contour_class.cont_closed;
										for (i = 0; i < num_pnt / 2; i++) {
											tmp = row[ i];
											row[ i] = row[ (num_pnt - 1 - i)];
											row[ (num_pnt - 1 - i)] = tmp;
											tmp = col[ i];
											col[ i] = col[ (num_pnt - 1 - i)];
											col[ (num_pnt - 1 - i)] = tmp;
											tmp = angle[ i];
											angle[ i] = angle[ (num_pnt - 1 - i)];
											angle[ (num_pnt - 1 - i)] = tmp;
											tmp = resp[ i];
											resp[ i] = resp[ (num_pnt - 1 - i)];
											resp[ (num_pnt - 1 - i)] = tmp;
										}
										it = 2;
									} else {
										if (it == 2) {
											/* Determine contour class. */
											if (cls == LinesUtil.contour_class.cont_start_junc)
												cls = LinesUtil.contour_class.cont_both_junc;
											else
												cls = LinesUtil.contour_class.cont_end_junc;
											/* Index j is the correct index. */
											junc[ num_junc].cont1 = num_cont;
											junc[ num_junc].cont2 = num_cont;
											junc[ num_junc].pos = j;
											junc[ num_junc].x = posx[ pos];
											junc[ num_junc].y = posy[ pos];
											num_junc++;
										} else {
											/* Determine contour class. */
											cls = LinesUtil.contour_class.cont_start_junc;
											/*
											 * Index num_pnt-1-j is the correct
											 * index since the line is going to
											 * be sorted in reverse.
											 */
											junc[ num_junc].cont1 = num_cont;
											junc[ num_junc].cont2 = num_cont;
											junc[ num_junc].pos = num_pnt
													- 1 - j;
											junc[ num_junc].x = posx[ pos];
											junc[ num_junc].y = posy[ pos];
											num_junc++;
										}
									}
									break;
								}
							}
							/*
							 * Mark this case as being processed for the
							 * algorithm below.
							 */
							j = -1;
						} else {

							for (j = 0; j < cont[ k].num; j++) {
								if (cont[ k].row[ j] == posx[ pos]
										&& cont[ k].col[ j] == posy[ pos])
									break;
							}
							/*
							 * If no point can be found on the other line a
							 * double response must have occured. In this case,
							 * find the nearest point on the other line and add
							 * it to the current line.
							 */
							if (j ==  cont[ k].num) {
								mindist = Double.MAX_VALUE;
								j = -1;
								for (l = 0; l < cont[ k].num; l++) {
									dx = posx[ pos]
											- cont[ k].row[ l];
									dy = posy[ pos]
											- cont[ k].col[ l];
									dist = Math.sqrt(dx * dx + dy * dy);
									if (dist < mindist) {
										mindist = dist;
										j = l;
									}
								}
								/*
								 * Add the point with index j to the current
								 * line.
								 */
								if (num_pnt >= size_pnt) {
									size_pnt =  (int)Math
											.floor((double) (size_pnt * LinesUtil.REALLOC_FACTOR));
									float[] newArr = new float[ size_pnt];
									for (int o = 0; o < row.length; o++) {
										newArr[o] = row[o];
									}
									row = newArr;
									newArr = new float[ size_pnt];
									for (int o = 0; o < col.length; o++) {
										newArr[o] = col[o];
									}
									col = newArr;
									newArr = new float[ size_pnt];
									for (int o = 0; o < angle.length; o++) {
										newArr[o] = angle[o];
									}
									angle = newArr;
									newArr = new float[ size_pnt];
									for (int o = 0; o < resp.length; o++) {
										resp[o] = resp[o];
									}
									resp = newArr;
								}

								row[ num_pnt] = cont[ k].row[ j];
								col[ num_pnt] = cont[ k].col[ j];
								beta = cont[ k].angle[ j];
								if (beta >= Math.PI)
									beta -= Math.PI;
								diff1 = Math.abs(beta - last_beta);
								if (diff1 >= Math.PI)
									diff1 = 2.0 * Math.PI - diff1;
								diff2 = Math.abs(beta + Math.PI - last_beta);
								if (diff2 >= Math.PI)
									diff2 = 2.0 * Math.PI - diff2;
								if (diff1 < diff2)
									angle[ num_pnt] = (float) beta;
								else
									angle[ num_pnt] = (float) (beta + Math.PI);
								resp[ num_pnt] =  cont[ k].response[ j];
								num_pnt++;
							}
						}
						/*
						 * Add the junction point only if it is not one of the
						 * other line's endpoints.
						 */
						if (j > 0 && j < cont[ k].num - 1) {
							/* Determine contour class. */
							if (it == 1)
								cls = LinesUtil.contour_class.cont_start_junc;
							else if (cls == LinesUtil.contour_class.cont_start_junc)
								cls = LinesUtil.contour_class.cont_both_junc;
							else
								cls = LinesUtil.contour_class.cont_end_junc;
							/* Add the new junction. */

							junc[ num_junc].cont1 = k;
							junc[ num_junc].cont2 = num_cont;
							junc[ num_junc].pos = j;
							junc[ num_junc].x = row[ (num_pnt - 1)];
							junc[ num_junc].y = col[ (num_pnt - 1)];
							num_junc++;
						}
						break;
					}
					label[ pos] =  (num_cont + 1);
					if (!(indx[ pos] == 0))
						cross[ (indx[ pos] - 1)].done = true;
				}
			}

			if (num_pnt > 1) {
				/* Only add lines with at least two points. */
				if (num_cont >= size_cont) {
					size_cont =  (int)Math
							.floor((double) (size_cont * LinesUtil.REALLOC_FACTOR));
					Line[] conth = new Line[ size_cont];
					for (int o = 0; o < conth.length; o++) {
						// true ? (conth[o] = cont[0]) : (conth[o] = new
						// contour());
						if (o < cont.length)
							conth[o] = cont[o];
						else
							conth[o] = new Line();
					}
					cont = conth;
				}
				cont[ num_cont] = new Line();
				
				cont[ num_cont].row = java.util.Arrays.copyOf(row,
						 num_pnt);
				cont[ num_cont].col = java.util.Arrays.copyOf(col,
						 num_pnt);
				cont[ num_cont].angle = java.util.Arrays.copyOf(angle,
						 num_pnt);
				cont[ num_cont].response = java.util.Arrays.copyOf(resp,
						 num_pnt);

				cont[ num_cont].width_r = null;
				cont[ num_cont].width_l = null;
				cont[ num_cont].asymmetry = null;
				cont[ num_cont].intensity = null;
				cont[ num_cont].num = num_pnt;
				cont[ num_cont].setContourClass(cls);
				num_cont++;
			} else {
				/*
				 * Delete the point from the label image; we can use maxx and
				 * maxy as the coordinates in the label image in this case.
				 */
				for (i = -1; i <= 1; i++) {
					for (j = -1; j <= 1; j++) {
						pos = LinesUtil.LINCOOR(LinesUtil.BR(maxx + i, height),
								LinesUtil.BC(maxy + j, width), width);
						if (label[ pos] == num_cont + 1)
							label[ pos] = 0;
					}
				}
			}
		}

		/*
		 * Now try to extend the lines at their ends to find additional
		 * junctions.
		 */
		if (extend_lines) {
			/* Sign by which the gradient has to be multiplied below. */
			if (mode == LinesUtil.MODE_LIGHT)
				s = 1;
			else
				s = -1;
			double MAX_LINE_EXTENSION = 2.5 * sigma;
			length = MAX_LINE_EXTENSION;
			max_line =  (int)Math.ceil(length * 3);
			line = new Offset[ max_line];
			for (int o = 0; o < line.length; o++) {
				line[o] = new Offset();
			}
			extx = new float[ max_line];
			exty = new float[ max_line];
			for (i = 0; i < num_cont; i++) {
				tmp_cont = cont[ i];
				num_pnt = tmp_cont.num;
				if (num_pnt == 1)
					continue;
				if (tmp_cont.getContourClass() == LinesUtil.contour_class.cont_closed)
					continue;
				trow = tmp_cont.row;
				tcol = tmp_cont.col;
				tangle = tmp_cont.angle;
				tresp = tmp_cont.response;
				/* Check both ends of the line (it==-1: start, it==1: end). */
				for (it = -1; it <= 1; it += 2) {
					/*
					 * Determine the direction of the search line. This is done
					 * by using the normal to the line (angle). Since this
					 * normal may point to the left of the line (see below) we
					 * have to check for this case by comparing the normal to
					 * the direction of the line at its respective end point.
					 */
					if (it == -1) {
						/* Start point of the line. */
						if (tmp_cont.getContourClass() == LinesUtil.contour_class.cont_start_junc
								|| tmp_cont.getContourClass() == LinesUtil.contour_class.cont_both_junc)
							continue;
						dx = trow[1] - trow[0];
						dy = tcol[1] - tcol[0];
						alpha = tangle[0];
						nx = Math.cos(alpha);
						ny = Math.sin(alpha);
						if (nx * dy - ny * dx < 0) {
							/* Turn the normal by +90 degrees. */
							mx = -ny;
							my = nx;
						} else {
							/* Turn the normal by -90 degrees. */
							mx = ny;
							my = -nx;
						}
						px = trow[0];
						py = tcol[0];
						response = tresp[0];
					} else {
						/* End point of the line. */
						if (tmp_cont.getContourClass() == LinesUtil.contour_class.cont_end_junc
								|| tmp_cont.getContourClass() == LinesUtil.contour_class.cont_both_junc)
							continue;
						dx = trow[ (num_pnt - 1)]
								- trow[ (num_pnt - 2)];
						dy = tcol[ (num_pnt - 1)]
								- tcol[ (num_pnt - 2)];
						alpha = tangle[ (num_pnt - 1)];
						nx = Math.cos(alpha);
						ny = Math.sin(alpha);
						if (nx * dy - ny * dx < 0) {
							/* Turn the normal by -90 degrees. */
							mx = ny;
							my = -nx;
						} else {
							/* Turn the normal by +90 degrees. */
							mx = -ny;
							my = nx;
						}
						px = trow[ (num_pnt - 1)];
						py = tcol[ (num_pnt - 1)];
						response = tresp[ (num_pnt - 1)];
					}
					/*
					 * Determine the current pixel and calculate the pixels on
					 * the search line.
					 */
					x =  (int)Math.floor(px + 0.5);
					y =  (int)Math.floor(py + 0.5);
					dx = px - x;
					dy = py - y;
					w.bresenham(mx, my, dx, dy, length, line, num_line);
					/*
					 * Now determine whether we can go only uphill (bright
					 * lines) or downhill (dark lines) until we hit another
					 * line.
					 */
					num_add = 0;
					add_ext = false;
					for (k = 0; k < num_line.intValue(); k++) {
						nextx = x + line[ k].x;
						nexty = y + line[ k].y;
						MutableDouble hnextpx = new MutableDouble(nextpx);
						MutableDouble hnextpy = new MutableDouble(nextpy);
						closest_point(px, py, mx, my, (double) nextx,
								(double) nexty, hnextpx, hnextpy, t);
						nextpx = hnextpx.getValue();
						nextpy = hnextpy.getValue();
						/*
						 * Ignore points before or less than half a pixel away
						 * from the true end point of the line.
						 */
						if (t.getValue() <= 0.5)
							continue;
						/*
						 * Stop if the gradient can't be interpolated any more
						 * or if the next point lies outside the image.
						 */
						if (nextpx < 0 || nextpy < 0 || nextpx >= height - 1
								|| nextpy >= width - 1 || nextx < 0
								|| nexty < 0 || nextx >= height
								|| nexty >= width)
							break;
						MutableDouble hgx = new MutableDouble();
						MutableDouble hgy = new MutableDouble();
						interpolate_gradient(gradx, grady, nextpx, nextpy,
								 width, hgx, hgy);
						gx = hgx.getValue();
						gy = hgy.getValue();
						/*
						 * Stop if we can't go uphill anymore. This is
						 * determined by the dot product of the line direction
						 * and the gradient. If it is smaller than 0 we go
						 * downhill (reverse for dark lines).
						 */
						nextpos = LinesUtil.LINCOOR(nextx, nexty, width);
						if (s * (mx * gx + my * gy) < 0
								&& label[ nextpos] == 0)
							break;
						/* Have we hit another line? */
						if (label[ nextpos] > 0) {
							m = label[ nextpos] - 1;
							/* Search for the junction point on the other line. */
							mindist = Double.MAX_VALUE;
							j = -1;
							for (l = 0; l < cont[ m].num; l++) {
								dx = nextpx - cont[ m].row[ l];
								dy = nextpy - cont[ m].col[ l];
								dist = Math.sqrt(dx * dx + dy * dy);
								if (dist < mindist) {
									mindist = dist;
									j = l;
								}
							}
							/*
							 * This should not happen... But better safe than
							 * sorry...
							 */
							if (mindist > 3.0){
								break;
							}
							extx[ num_add] = cont[ m].row[ j];
							exty[ num_add] = cont[ m].col[ j];
							end_resp = cont[ m].response[ j];
							end_angle = cont[ m].angle[ j];
							beta = end_angle;
							if (beta >= Math.PI)
								beta -= Math.PI;
							diff1 = Math.abs(beta - alpha);
							if (diff1 >= Math.PI)
								diff1 = 2.0 * Math.PI - diff1;
							diff2 = Math.abs(beta + Math.PI - alpha);
							if (diff2 >= Math.PI)
								diff2 = 2.0 * Math.PI - diff2;
							if (diff1 < diff2)
								end_angle = beta;
							else
								end_angle = beta + Math.PI;
							num_add++;
							add_ext = true;
							break;
						} else {
							extx[ num_add] = (float) nextpx;
							exty[ num_add] = (float) nextpy;
							num_add++;
						}
					}
					if (add_ext) {
						/* Make room for the new points. */
						num_pnt += num_add;
						float[] newArr = new float[ num_pnt];
						for (int o = 0; o < trow.length; o++) {
							newArr[o] = trow[o];
						}
						trow = newArr;

						newArr = new float[ num_pnt];
						for (int o = 0; o < tcol.length; o++) {
							newArr[o] = tcol[o];
						}
						tcol = newArr;

						newArr = new float[ num_pnt];
						for (int o = 0; o < tangle.length; o++) {
							newArr[o] = tangle[o];
						}
						tangle = newArr;

						newArr = new float[ num_pnt];
						for (int o = 0; o < tresp.length; o++) {
							newArr[o] = tresp[o];
						}
						tresp = newArr;

						tmp_cont.row = trow;
						tmp_cont.col = tcol;
						tmp_cont.angle = tangle;
						tmp_cont.response = tresp;
						tmp_cont.num = num_pnt;
						if (it == -1) {
							/* Move points on the line up num_add places. */
							for (k = num_pnt - 1 - num_add; k >= 0; k--) {
								trow[ (k + num_add)] = trow[ k];
								tcol[ (k + num_add)] = tcol[ k];
								tangle[ (k + num_add)] = tangle[ k];
								tresp[ (k + num_add)] = tresp[ k];
							}
							/* Insert points at the beginning of the line. */
							for (k = 0; k < num_add; k++) {
								trow[ k] = extx[ (num_add - 1 - k)];
								tcol[ k] = exty[ (num_add - 1 - k)];
								tangle[ k] = (float) alpha;
								tresp[ k] = (float) response;
							}
							tangle[0] = (float) end_angle;
							tresp[0] = (float) end_resp;
							/* Adapt indices of the previously found junctions. */
							for (k = 0; k < num_junc; k++) {
								if (junc[ k].cont1 == i)
									junc[ k].pos += num_add;
							}
						} else {
							/* Insert points at the end of the line. */
							for (k = 0; k < num_add; k++) {
								trow[ (num_pnt - num_add + k)] = extx[ k];
								tcol[ (num_pnt - num_add + k)] = exty[ k];
								tangle[ (num_pnt - num_add + k)] = (float) alpha;
								tresp[ (num_pnt - num_add + k)] = (float) response;
							}
							tangle[ (num_pnt - 1)] = (float) end_angle;
							tresp[ (num_pnt - 1)] = (float) end_resp;
						}
						/* If necessary, make room for the new junction. */
						if (num_junc >= size_junc) {
							size_junc =  (int)Math
									.floor((double) (size_junc * LinesUtil.REALLOC_FACTOR));
							Junction[] junch = new Junction[ size_junc];
							for (int o = 0; o < junch.length; o++) {
								if (o < junc.length)
									junch[o] = junc[o];
								else
									junch[o] = new Junction();
							}
							junc = junch;

						}
						/*
						 * Add the junction point only if it is not one of the
						 * other line's endpoints.
						 */
						if (j > 0 && j < cont[ m].num - 1) {
							if (it == -1) {
								if (tmp_cont.getContourClass() == LinesUtil.contour_class.cont_end_junc)
									tmp_cont.setContourClass(LinesUtil.contour_class.cont_both_junc);
								else
									tmp_cont.setContourClass(LinesUtil.contour_class.cont_start_junc);
							} else {
								if (tmp_cont.getContourClass() == LinesUtil.contour_class.cont_start_junc)
									tmp_cont.setContourClass(LinesUtil.contour_class.cont_both_junc);
								else
									tmp_cont.setContourClass(LinesUtil.contour_class.cont_end_junc);
							}
							junc[ num_junc].cont1 = m;
							junc[ num_junc].cont2 = i;
							junc[ num_junc].pos = j;
							if (it == -1) {
								junc[ num_junc].x = trow[0];
								junc[ num_junc].y = tcol[0];
							} else {
								junc[ num_junc].x = trow[ (num_pnt - 1)];
								junc[ num_junc].y = tcol[ (num_pnt - 1)];
							}
							num_junc++;
						}
					}
				}
				
			}
		}

		/* Done with linking. Now split the lines at the junction points. */
		java.util.Arrays.sort(junc);
		for (i = 0; i < num_junc; i += k) {
			j = junc[ i].cont1;
			tmp_cont = cont[ j];
			num_pnt = tmp_cont.num;
			/* Count how often line j needs to be split. */
			for (k = 0; junc[ (i + k)].cont1 == j && i + k < num_junc; k++);
				
			if (k == 1 && tmp_cont.row[0] == tmp_cont.row[ (num_pnt - 1)]
					&& tmp_cont.col[0] == tmp_cont.col[ (num_pnt - 1)]) {
				/*
				 * If only one junction point is found and the line is closed it
				 * only needs to be rearranged cyclically, but not split.
				 */
				begin = junc[ i].pos;
				trow = tmp_cont.row;
				tcol = tmp_cont.col;
				tangle = tmp_cont.angle;
				tresp = tmp_cont.response;
				tmp_cont.row = new float[ num_pnt];
				tmp_cont.col = new float[ num_pnt];
				tmp_cont.angle = new float[ num_pnt];
				tmp_cont.response = new float[ num_pnt];
				for (l = 0; l < num_pnt; l++) {
					pos = begin + l;
					/* Skip starting point so that it is not added twice. */
					if (pos >= num_pnt)
						pos = begin + l - num_pnt + 1;
					tmp_cont.row[ l] = trow[ pos];
					tmp_cont.col[ l] = tcol[ pos];
					tmp_cont.angle[ l] = tangle[ pos];
					tmp_cont.response[ l] = tresp[ pos];
				}
				/* Modify contour class. */
				tmp_cont.setContourClass(LinesUtil.contour_class.cont_both_junc);

			} else {
				/* Otherwise the line has to be split. */
				for (l = 0; l <= k; l++) {
					if (l == 0)
						begin = 0;
					else
						begin = junc[ (i + l - 1)].pos;
					if (l == k)
						end = tmp_cont.num - 1;
					else
						end = junc[ (i + l)].pos;
					num_pnt = end - begin + 1;
					if (num_pnt == 1 && k > 1) {
						/* Do not add one point segments. */
						continue;
					}
					if (num_cont >= size_cont) {
						size_cont =  (int)Math
								.floor((double) (size_cont * LinesUtil.REALLOC_FACTOR));
						Line[] conth = new Line[ size_cont];
						for (int o = 0; o < cont.length; o++) {
								conth[o] = cont[o];
						}
						cont = conth;
					}
					cont[ num_cont] = new Line();
					
					cont[ num_cont].row = new float[ num_pnt];
					cont[ num_cont].col = new float[ num_pnt];
					cont[ num_cont].angle = new float[ num_pnt];
					cont[ num_cont].response = new float[ num_pnt];

					System.arraycopy(tmp_cont.row,  begin,
							cont[ num_cont].row, 0,  num_pnt);

					System.arraycopy(tmp_cont.col,  begin,
							cont[ num_cont].col, 0,  num_pnt);

					System.arraycopy(tmp_cont.angle,  begin,
							cont[ num_cont].angle, 0,  num_pnt);

					System.arraycopy(tmp_cont.response,  begin,
							cont[ num_cont].response, 0,  num_pnt);

					cont[ num_cont].width_r = null;
					cont[ num_cont].width_l = null;
					cont[ num_cont].asymmetry = null;
					cont[ num_cont].intensity = null;
					cont[ num_cont].num = num_pnt;
					/* Modify contour class. */
					if (l == 0) {
						if (tmp_cont.getContourClass() == LinesUtil.contour_class.cont_start_junc
								|| tmp_cont.getContourClass() == LinesUtil.contour_class.cont_both_junc)
							cont[ num_cont].setContourClass(LinesUtil.contour_class.cont_both_junc);
						else
							cont[ num_cont].setContourClass(LinesUtil.contour_class.cont_end_junc);
					} else if (l == k) {
						if (tmp_cont.getContourClass() == LinesUtil.contour_class.cont_end_junc
								|| tmp_cont.getContourClass()== LinesUtil.contour_class.cont_both_junc)
							cont[ num_cont].setContourClass(LinesUtil.contour_class.cont_both_junc);
						else
							cont[ num_cont].setContourClass(LinesUtil.contour_class.cont_start_junc);
					} else {
						cont[ num_cont].setContourClass(LinesUtil.contour_class.cont_both_junc);
					}
					num_cont++;
				}
				cont[ j] = cont[ --num_cont];
			}
		}

		/* Finally, check whether all angles point to the right of the line. */
		for (i = 0; i < num_cont; i++) {
			tmp_cont = cont[ i];
			num_pnt = tmp_cont.num;
			if (num_pnt > 1) {
				trow = tmp_cont.row;
				tcol = tmp_cont.col;
				tangle = tmp_cont.angle;
				/*
				 * One point of the contour is enough to determine the
				 * orientation.
				 */
				
				k = (num_pnt - 1) / 2;
				/*
				 * The next few lines are ok because lines have at least two
				 * points.
				 */
				dx = trow[ (k + 1)] - trow[ k];
				dy = tcol[ (k + 1)] - tcol[ k];
				nx = Math.cos(tangle[ k]);
				ny = Math.sin(tangle[ k]);
				/*
				 * If the angles point to the left of the line they have to be
				 * adapted. The orientation is determined by looking at the
				 * z-component of the cross-product of (dx,dy,0) and (nx,ny,0).
				 */
				if (nx * dy - ny * dx < 0) {
					for (j = 0; j < num_pnt; j++) {
						tangle[ j] += Math.PI;
						if (tangle[ j] >= 2 * Math.PI)
							tangle[ j] -= 2 * Math.PI;
					}
				}
			}
		}

		for (Line c : cont) {
			if (c!=null && c.getContourClass()!= null) {
				c.setFrame(contours.getFrame());
				contours.add(c);

			}
		}
	
		for (Junction jun : junc) {
			if (jun != null && !(jun.cont1 == 0 && jun.cont2 == 0)) {
				junctions.add(jun);
			}
		}
		num_result.setValue(num_cont);
	}

}
