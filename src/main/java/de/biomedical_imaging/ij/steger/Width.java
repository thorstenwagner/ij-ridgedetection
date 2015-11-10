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

import java.util.ArrayList;

import org.apache.commons.lang3.mutable.MutableDouble;
import org.apache.commons.lang3.mutable.MutableInt;


public class Width {

	/* Maximum search line length. */
	//#define MAX_LINE_WIDTH (2.5*sigma)

	/* This constant is introduced because for very narrow lines the facet model
	   width detection scheme sometimes extracts the line width too narrow.  Since
	   the correction function has a very steep slope in that area, this will lead
	   to lines of almost zero width, especially since the bilinear interpolation
	   in correct.c will tend to overcorrect.  Therefore it is wise to make the
	   extracted line width slightly larger before correction.  */
	public static final double LINE_WIDTH_COMPENSATION = 1.05;

	/* Minimum line width allowed (used for outlier check in fix_locations()) */
	public static final double MIN_LINE_WIDTH =  0.1;

	/* Maximum contrast allowed (used for outlier check in fix_locations()) */
	public static final double MAX_CONTRAST = 275.0;
	
	
	/* Modified Bresenham algorithm.  It returns in line all pixels that are
	   intersected by a half line less than length away from the point (px,py)
	   aint the direction (nx,ny).  The point (px,py) must lie within the pixel
	   of the origin, i.e., fabs(px) <= 0.5 and fabs(py) <= 0.5. */
	public void bresenham(double nx, double ny, double px, double py, double length, Offset[] line, MutableInt num_points)
	{
	  int i, n, x, y, s1, s2, xchg, maxit;
	  double e, dx, dy, t;

	  x = 0;
	  y = 0;
	  dx = Math.abs(nx);
	  dy = Math.abs(ny);
	  s1 =  (int)Math.signum(nx);
	  s2 =  (int)Math.signum(ny);
	  px *= s1;
	  py *= s2;
	  if (dy > dx) {
	    t = dx;
	    dx = dy;
	    dy = t;
	    t = px;
	    px = py;
	    py = t;
	    xchg = 1;
	  } else {
	    xchg = 0;
	  }
	  maxit =  (int)Math.ceil(length*dx);
	  e = (0.5-px)*dy/dx-(0.5-py);
	  n = 0;
	  for (i=0; i<=maxit; i++) {
	    line[n].x = x;
	    line[n].y = y;
	    n++;
	    while (e >= -1e-8) {
	      if (!(xchg==0)) x += s1;
	      else y += s2;
	      e--;
	      if (e > -1) {
	        line[n].x = x;
	        line[n].y = y;
	        n++;
	      }
	    }
	    if (!(xchg==0)) y += s2;
	    else x += s1;
	    e += dy/dx;
	  }
	  num_points.setValue(n);
	}
	
	
	/* Fill gaps in the arrays master, slave1, and slave2, i.e., points where
	   master=0, by interpolation (interior points) or extrapolation (end points).
	   The array master will usually be the width of the line, while slave1 and
	   slave2 will be values that depend on master[i] being 0, e.g., the gradient
	   at each line point.  The arrays slave1 and slave2 can be NULL. */
	private void fill_gaps(double[] master, double[] slave1, double[] slave2, Line cont)
	{
	  int    i, j, k, s, e;
	  int    num_points;
	  double  m_s, m_e, s1_s, s1_e, s2_s, s2_e, d_r, d_c, arc_len, len;

	  num_points = cont.num;
	  for (i=0; i<num_points; i++) {
	    if (master[ i] == 0) {
	      for (j=i+1; j<num_points; j++) {
	        if (master[ j] > 0)
	          break;
	      }
	      m_s = 0;
	      m_e = 0;
	      s1_s = 0;
	      s1_e = 0;
	      s2_s = 0;
	      s2_e = 0;
	      if (i > 0 && j < num_points-1) {
	        s = i;
	        e = j-1;
	        m_s = master[ (s-1)];
	        m_e = master[ (e+1)];
	        if (slave1 != null) {
	          s1_s = slave1[ (s-1)];
	          s1_e = slave1[ (e+1)];
	        }
	        if (slave2 != null) {
	          s2_s = slave2[ (s-1)];
	          s2_e = slave2[ (e+1)];
	        }
	      } else if (i > 0) {
	        s = i;
	        e = num_points-2;
	        m_s = master[ (s-1)];
	        m_e = master[ (s-1)];
	        master[ (e+1)] = m_e;
	        if (slave1 != null) {
	          s1_s = slave1[ (s-1)];
	          s1_e = slave1[ (s-1)];
	          slave1[ (e+1)] = s1_e;
	        }
	        if (slave2 != null) {
	          s2_s = slave2[ (s-1)];
	          s2_e = slave2[ (s-1)];
	          slave2[ (e+1)] = s2_e;
	        }
	      } else if (j < num_points-1) {
	        s = 1;
	        e = j-1;
	        m_s = master[ (e+1)];
	        m_e = master[ (e+1)];
	        master[ (s-1)] = m_s;
	        if (slave1 != null) {
	          s1_s = slave1[ (e+1)];
	          s1_e = slave1[ (e+1)];
	          slave1[ (s-1)] = s1_s;
	        }
	        if (slave2 != null) {
	          s2_s = slave2[ (e+1)];
	          s2_e = slave2[ (e+1)];
	          slave2[ (s-1)] = s2_s;
	        }
	      } else {
	        s = 1;
	        e = num_points-2;
	        m_s = master[ (s-1)];
	        m_e = master[ (e+1)];
	        if (slave1 != null) {
	          s1_s = slave1[ (s-1)];
	          s1_e = slave1[ (e+1)];
	        }
	        if (slave2 != null) {
	          s2_s = slave2[ (s-1)];
	          s2_e = slave2[ (e+1)];
	        }
	      }
	      arc_len = 0;
	      for (k=s; k<=e+1; k++) {
	        d_r = cont.row[ k]-cont.row[ (k-1)];
	        d_c = cont.col[ k]-cont.col[ (k-1)];
	        arc_len += Math.sqrt(d_r*d_r+d_c*d_c);
	      }
	      len = 0;
	      for (k=s; k<=e; k++) {
	        d_r = cont.row[ k]-cont.row[ (k-1)];
	        d_c = cont.col[ k]-cont.col[ (k-1)];
	        len += Math.sqrt(d_r*d_r+d_c*d_c);
	        master[ k] = (arc_len-len)/arc_len*m_s+len/arc_len*m_e;
	        if (slave1 != null)
	          slave1[ k] = (arc_len-len)/arc_len*s1_s+len/arc_len*s1_e;
	        if (slave2 != null)
	          slave2[ k] = (arc_len-len)/arc_len*s2_s+len/arc_len*s2_e;
	      }
	      i = j;
	    }
	  }
	}
	
	
	/* Correct the extracted line positions and widths.  The algorithm first closes
	   gaps in the extracted data width_l, width_r, grad_l, and grad_r to provide
	   meaningful input over the whole line.  Then the correction is calculated.
	   After this, gaps that have been introduced by the width correction are again
	   closed.  Finally, the position correction is applied if correct_pos is set.
	   The results are returned in width_l, width_r, and cont. */
	private void fix_locations(double[] width_l,double[] width_r, double[] grad_l,double[] grad_r, double[] pos_x, double[] pos_y,
	                          double[] correction,double[] contr, double[] asymm, double sigma, int mode, boolean correct_pos, Line cont)
	{
	  int    i;
	  int    num_points;
	  double  px, py;
	  double  nx, ny;
	  double  w_est, r_est;
	  MutableDouble w_real, h_real, corr;
	  w_real = new MutableDouble();
	  h_real = new MutableDouble();
	  corr = new MutableDouble();
	  MutableDouble w_strong = new MutableDouble();
	  MutableDouble w_weak = new MutableDouble();
	  double  correct, asymmetry, response, width, contrast;
	  boolean    weak_is_r;
	  boolean    correct_start, correct_end;
	  Convol convol = new Convol();
	  fill_gaps(width_l,grad_l,null,cont);
	  fill_gaps(width_r,grad_r,null,cont);

	  num_points = cont.num;

	  /* Calculate true line width, asymmetry, and position correction. */
	  if (correct_pos) {
	    /* Do not correct the position of a junction point if its width is found
	       by interpolation, i.e., if the position could be corrected differently
	       for each junction point, thereby destroying the junction. */
	    correct_start = ((cont.getContourClass() == LinesUtil.contour_class.cont_no_junc ||
	                      cont.getContourClass() == LinesUtil.contour_class.cont_end_junc ||
	                      cont.getContourClass() == LinesUtil.contour_class.cont_closed) &&
	                     (width_r[0] > 0 && width_l[0] > 0));
	    correct_end = ((cont.getContourClass() == LinesUtil.contour_class.cont_no_junc ||
	                    cont.getContourClass() == LinesUtil.contour_class.cont_start_junc ||
	                    cont.getContourClass() == LinesUtil.contour_class.cont_closed) &&
	                   (width_r[ (num_points-1)] > 0 && width_l[ (num_points-1)] > 0));
	    /* Calculate the true width and assymetry, and its corresponding
	       correction for each line point. */
	    for (i=0; i<num_points; i++) {
	      if (width_r[ i] > 0 && width_l[ i] > 0) {
	        w_est = (width_r[ i]+width_l[ i])*LINE_WIDTH_COMPENSATION;
	        if (grad_r[ i] <= grad_l[ i]) {
	          r_est = grad_r[ i]/grad_l[ i];
	          weak_is_r = true;
	        } else {
	          r_est = grad_l[ i]/grad_r[ i];
	          weak_is_r = false;
	        }
	        Correct.line_corrections(sigma,w_est,r_est,w_real,h_real,corr,
	                         w_strong,w_weak);
	        w_real.setValue(w_real.getValue()/LINE_WIDTH_COMPENSATION);
	        corr.setValue(corr.getValue()/LINE_WIDTH_COMPENSATION);
	        width_r[ i] = w_real.getValue();
	        width_l[ i] = w_real.getValue();
	        if (weak_is_r) {
	          asymm[ i] = h_real.getValue();
	          correction[ i] = -corr.getValue();
	        } else {
	          asymm[ i] = -h_real.getValue();
	          correction[ i] = corr.getValue();
	        }
	      }
	    }

	    fill_gaps(width_l,correction,asymm,cont);
	    for (i=0; i<num_points; i++)
	      width_r[ i] = width_l[ i];

	    /* Adapt the correction for junction points if necessary. */
	    if (!correct_start)
	      correction[0] = 0;
	    if (!correct_end)
	      correction[ (num_points-1)] = 0;

	    for (i=0; i<num_points; i++) {
	      px = pos_x[ i];
	      py = pos_y[ i];
	      nx = Math.cos(cont.angle[ i]);
	      ny = Math.sin(cont.angle[ i]);
	      px = px+correction[ i]*nx;
	      py = py+correction[ i]*ny;
	      pos_x[ i] = px;
	      pos_y[ i] = py;
	    }
	  }

	  /* Update the position of a line and add the extracted width. */
	  cont.width_l = new float[num_points];
	  cont.width_r = new float[num_points];
	  for (i=0; i<num_points; i++) {
	    cont.width_l[ i] = (float) width_l[ i];
	    cont.width_r[ i] = (float) width_r[ i];
	    cont.row[ i] = (float) pos_x[ i];
	    cont.col[ i] = (float) pos_y[ i];
	  }

	  /* Now calculate the true contrast. */
	  if (correct_pos) {
	    cont.asymmetry = new float[num_points];
	    cont.intensity = new float[num_points];
	    for (i=0; i<num_points; i++) {
	      response = cont.response[ i];
	      asymmetry = Math.abs(asymm[ i]);
	      correct = Math.abs(correction[ i]);
	      width = cont.width_l[ i];
	      if (width < MIN_LINE_WIDTH)
	        contrast = 0;
	      else
	        contrast = 
	          (response/Math.abs(convol.phi2(correct+width,sigma)+
	                         (asymmetry-1)*convol.phi2(correct-width,sigma)));
	      
	      if (contrast > MAX_CONTRAST)
	        contrast = 0;
	      contr[ i] = contrast;
	    }
	    fill_gaps(contr,null,null,cont);
	    for (i=0; i<num_points; i++) {
	      cont.asymmetry[ i] = (float) asymm[ i];
	      if (mode == LinesUtil.MODE_LIGHT)
	        cont.intensity[ i] = (float) contr[ i];
	      else
	        cont.intensity[ i] = (float) -contr[ i];
	    }
	  }
	}
	
	/* Extract the line width by using a facet model line detector on an image of
	   the absolute value of the gradient. */
	public void compute_line_width(float[] dx, float[] dy, int width, int height, double sigma,int mode,boolean correct_pos, ArrayList<Line> contours,
	                        MutableInt num_contours)
	{
	  float[] grad;
	  int    i, j, k;
	  int    r, c, l;
	  int    x, y, dir;
	  Offset[]  line;
	  int    max_line, num_line=0;
	  double  length;
	  Line cont;
	  int    num_points, max_num_points;
	  double[]  width_r, width_l;
	  double[]  grad_r, grad_l;
	  double[]  pos_x, pos_y, correct, asymm, contrast;
	  double  d, dr, dc, drr, drc, dcc;
	  double  i1, i2, i3, i4, i5, i6, i7, i8, i9;
	  double  t1, t2, t3, t4, t5, t6;
	  double[]  eigval = new double[2];
	  double[][]  eigvec = new double[2][2];
	  double  a, b, t=0;
	  int    num = 0;
	  double  nx, ny;
	  double  n1, n2;
	  double  p1, p2;
	  double  val;
	  double  px, py;
	  Position p = new Position();
	  max_num_points = 0;
	  for (i=0; i<num_contours.getValue(); i++) {
	    num_points = contours.get(i).num;
	    if (num_points > max_num_points)
	      max_num_points = num_points;
	  }

	  width_l = new double[max_num_points];
	  width_r = new double[ max_num_points];
	  grad_l = new double[ max_num_points];
	  grad_r = new double[ max_num_points];
	  pos_x = new double[ max_num_points];
	  pos_y = new double[ max_num_points];
	  correct = new double[ max_num_points];
	  contrast = new double[ max_num_points];
	  asymm = new double[ max_num_points];

	  grad = new float[ (width*height)];

	  length = 2.5*sigma;
	  max_line =  (int)Math.ceil(length*3);
	  line = new Offset[ max_line];
	  for(int o = 0; o < line.length; o++){
		  line[o] = new Offset();
	  }

	  /* Compute the gradient image. */
	  for (r=0; r<height; r++) {
	    for (c=0; c<width; c++) {
	      l = LinesUtil.LINCOOR(r,c,width);
	      grad[ l] = (float) Math.sqrt(dx[ l]*dx[ l]+dy[ l]*dy[ l]);
	    }
	  }

	  for (i=0; i<num_contours.getValue(); i++) {
	    cont = contours.get(i);
	    num_points = cont.num;
	    
	    for (j=0; j<num_points; j++) {
	      px = cont.row[ j];
	      py = cont.col[ j];
	      pos_x[ j] = px;
	      pos_y[ j] = py;
	      r =  (int)Math.floor(px+0.5);
	      c =  (int)Math.floor(py+0.5);
	      nx = Math.cos(cont.angle[ j]);
	      ny = Math.sin(cont.angle[ j]);
	      /* Compute the search line. */
	      MutableInt num_lineh = new MutableInt(num_line);
	      bresenham(nx,ny,0.0,0.0,length,line,num_lineh);
	      num_line = num_lineh.intValue();
	      width_r[ j] = width_l[ j] = 0;
	      /* Look on both sides of the line. */
	      for (dir=-1; dir<=1; dir+=2) {
	        for (k=0; k<num_line; k++) {
	          x = LinesUtil.BR(r+dir*line[ k].x,height);
	          y = LinesUtil.BC(c+dir*line[ k].y,width);
	          i1 = grad[ LinesUtil.LINCOOR(LinesUtil.BR(x-1,height),LinesUtil.BC(y-1,width),width)];
	          i2 = grad[ LinesUtil.LINCOOR(LinesUtil.BR(x-1,height),y,width)];
	          i3 = grad[ LinesUtil.LINCOOR(LinesUtil.BR(x-1,height),LinesUtil.BC(y+1,width),width)];
	          i4 = grad[ LinesUtil.LINCOOR(x,LinesUtil.BC(y-1,width),width)];
	          i5 = grad[ LinesUtil.LINCOOR(x,y,width)];
	          i6 = grad[ LinesUtil.LINCOOR(x,LinesUtil.BC(y+1,width),width)];
	          i7 = grad[ LinesUtil.LINCOOR(LinesUtil.BR(x+1,height),LinesUtil.BC(y-1,width),width)];
	          i8 = grad[ LinesUtil.LINCOOR(LinesUtil.BR(x+1,height),y,width)];
	          i9 = grad[ LinesUtil.LINCOOR(LinesUtil.BR(x+1,height),LinesUtil.BC(y+1,width),width)];
	          t1 = i1+i2+i3;
	          t2 = i4+i5+i6;
	          t3 = i7+i8+i9;
	          t4 = i1+i4+i7;
	          t5 = i2+i5+i8;
	          t6 = i3+i6+i9;
	          dr = (t3-t1)/6;
	          dc = (t6-t4)/6;
	          drr = (t1-2*t2+t3)/6;
	          dcc = (t4-2*t5+t6)/6;
	          drc = (i1-i3-i7+i9)/4;
	          p.compute_eigenvals(2*drr,drc,2*dcc,eigval,eigvec);
	          val = -eigval[0];
	          if (val > 0.0) {
	            n1 = eigvec[0][0];
	            n2 = eigvec[0][1];
	            a = 2.0*(drr*n1*n1+drc*n1*n2+dcc*n2*n2);
	            b = dr*n1+dc*n2;
	            MutableDouble th = new MutableDouble(t);
	            MutableInt numh = new MutableInt(num);
	            p.solve_linear(a,b,th,numh);
	            t = th.getValue();
	            num = numh.getValue();
	            if (num != 0) {
	              p1 = t*n1;
	              p2 = t*n2;
	              if (Math.abs(p1) <= 0.5 && Math.abs(p2) <= 0.5) {
	                /* Project the maximum point position perpendicularly onto the
	                   search line. */
	                a = 1;
	                b = nx*(px-(r+dir*line[ k].x+p1))+ny*(py-(c+dir*line[ k].y+p2));
	                th = new MutableDouble(t);
		            numh = new MutableInt(num);
		            p.solve_linear(a,b,th,numh);
		            t = th.getValue();
		            num = numh.getValue();
	                d = (-i1+2*i2-i3+2*i4+5*i5+2*i6-i7+2*i8-i9)/9;
	                if (dir == 1) {
	                  grad_r[ j] = d+p1*dr+p2*dc+p1*p1*drr+p1*p2*drc+p2*p2*dcc;
	                  width_r[ j] = Math.abs(t);
	                } else {
	                  grad_l[ j] = d+p1*dr+p2*dc+p1*p1*drr+p1*p2*drc+p2*p2*dcc;
	                  width_l[ j] = Math.abs(t);
	                }
	                break;
	              }
	            }
	          }
	        }
	      }
	    }

	    fix_locations(width_l,width_r,grad_l,grad_r,pos_x,pos_y,correct,contrast,
	                  asymm,sigma,mode,correct_pos,cont);
	  }
	}

}
