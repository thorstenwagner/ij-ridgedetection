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

public class Threshold {
	
	static void threshold(byte[] image, int min,int width,int height,Region out)
	{
	  int   grey;
	  int   r,c,l,num,num_max;
	  boolean   inside;
	  Chord[] rl;

	  inside = false;
	  num = 0;
	  num_max = LinesUtil.INITIAL_SIZE;
	  rl =new Chord[ num_max]; 
	  for(int i = 0; i < num_max; i++){
		  rl[i] = new Chord();
	  }
	  out.rl = null;
	  out.num = 0;

	  for (r=0; r<height; r++) {
	    for (c=0; c<width; c++) {
	      l = LinesUtil.LINCOOR(r,c,width);
	      grey = image[ l];
	      if (grey >= min) {
	        if (!inside) {
	          inside = true;
	          rl[ num].r = (short) r;
	          rl[ num].cb = (short) c;
	        }
	      } else {
	        if (inside) {
	          inside = false;
	          rl[ num].ce = (short) (c - 1);
	          num++;
	          if (num >= num_max) {
	            num_max =  (int)Math.floor((double)(num_max*LinesUtil.REALLOC_FACTOR));
	            Chord[] rlh = new Chord[ num_max];
	            for(int i = 0; i < rlh.length; i++){
	            	if(i < rl.length)
		        		rlh[i] = rl[i];
		        	else
		        		rlh[i] = new Chord();
	      	  	}
	            rl=rlh;
	          }
	        }
	      }
	    }
	    if (inside) {
	      inside = false;
	      rl[ num].ce = (short) (width-1);
	      num++;
	      if (num >= num_max) {
	        num_max =  (int)Math.floor((double)(num_max*LinesUtil.REALLOC_FACTOR));
	        Chord[] rlh = new Chord[ num_max];
	        for(int i = 0; i < rlh.length; i++){
	        	if(i < rl.length)
	        		rlh[i] = rl[i];
	        	else
	        		rlh[i] = new Chord();
	        }
	        rl = rlh;
	      }
	    }
	  }
	  out.rl = new Chord[ num];
	  for(int i = 0; i < num; i++){
		  out.rl[i] = rl[i];
	  }
	  out.num = num;
	}

}
